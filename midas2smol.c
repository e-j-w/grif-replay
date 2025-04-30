/* read midas data files
 set sources=( midas2smol.c midas-format.c grif-format.c config.c reorder.c user_sort.c default_sort.c test_config.c )
 gcc -g     -o midas2smol $sources -rdynamic -ldl -lm -lpthread
 gcc -g -O3 -o midas2smol $sources -rdynamic -ldl -lm -lpthread
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include "config.h"
#include "grif-format.h"
#include "midas-format.h"

int sort_next_file(Config *cfg, Sort_status *sort);

int coinc_events_cutoff = 64;
char midas_runtitle[SYS_PATH_LENGTH];
Sort_metrics diagnostics;
static Sort_status sort_status;
FILE *output_tree;
int main(int argc, char *argv[])
{
   if(argc != 3){
      fprintf(stdout,"  midas2smol midas_file output_SMOL_tree\n");
      fprintf(stdout,"    MIDAS files are expected to be named using the\n");
      fprintf(stdout,"    standard run and subrun numbering scheme (eg.\n");
      fprintf(stdout,"    run29623_000.mid). Provide the path of the first \n");
      fprintf(stdout,"    subrun, and all available subruns for that run\n");
      fprintf(stdout,"    number will then be sorted.\n\n");
      exit(0);
   }
   srand(28719747); //seed random number generator (use fixed seed so that results are consistent between sorts)
   Sort_status *sort = &sort_status;
   int web_arg=1;  Config *cfg;
   init_config();
   copy_config(configs[0], configs[1]); // copy config0 to cfg1 for sorting
   cfg = configs[1];
   sort->subrun = 0;
   sort->single_thread = 0;

   strncpy(cfg->out_file,argv[2],SYS_PATH_LENGTH); //setup output filename
   add_sortfile(argv[1]);
   fprintf(stdout,"MIDAS file: %s (%i subrun(s)), output file: %s\n",argv[1],sort->num_subruns,cfg->out_file);
   
   if( open_next_sortfiles(sort) == 0 ){
      sort_next_file(cfg, sort);
      fprintf(stdout,"DONE\n");
      close_sortfiles(sort);
   }
   
   exit(0);
}

static int presort_window_start, sort_window_start;
static int done_events;
extern void grif_main(Sort_status *arg);
extern void reorder_main(Sort_status *arg);
extern void reorder_out(Sort_status *arg);
extern uint64_t sort_main(Sort_status *arg, FILE *out);
static pthread_t midas_thread, grif_thread, ordthrd, ordthr2;
static int reorder_save, singlethread_save, sortthread_save;
int sort_next_file(Config *cfg, Sort_status *sort)
{
   char tmp[64];
   FILE *smolfp;
   time_t end, start=time(NULL);
   done_events = 0;
   presort_window_start = sort_window_start = 0;
   memset(&diagnostics, 0, sizeof(Sort_metrics) );
   diagnostics.run_sort_start = start;
   sort->shutdown_midas = sort->end_of_data = 0;
   sort->reorder_in_done = sort->reorder_out_done = 0;
   sort->grif_sort_done = sort->odb_done = sort->odb_ready = 0;
   reorder_save = sort->reorder;
   singlethread_save = sort->single_thread;
   sortthread_save = sort->sort_thread;
   if( singlethread_save == 1 ){
      midas_main(sort);
   } else {
      printf("creating midas thread\n");
      pthread_create(&midas_thread, NULL, (void* (*)(void*))midas_main, sort);
      printf("creating reorder threads\n");
      pthread_create(&ordthrd,NULL,(void* (*)(void*))reorder_main, sort);
      pthread_create(&ordthr2,NULL,(void* (*)(void*))reorder_out, sort);

      while( !sort->odb_ready ){     // wait for midas thread to read odb event
         usleep(1);
      }
      
      init_default_sort(configs[1], sort);     // depend on odb in datafile
      pthread_create(&grif_thread, NULL,(void* (*)(void*))grif_main, sort);

      sprintf(tmp,"%s", cfg->out_file);
      if( (smolfp=fopen(tmp,"wb")) != NULL ){
         printf("Opened output file: %s\n",tmp);
         uint64_t numSortedEvts = 0U; //placeholder
         fwrite(&numSortedEvts,sizeof(uint64_t),1,smolfp);
         numSortedEvts += sort_main(sort,smolfp); // this exits when sort is done
         //number of sorted events can only be 48 bits
         if(numSortedEvts > 0xFFFFFFFFFFFF){
            printf("WARNING: number of output events (%lu) truncated to %lu.\n",numSortedEvts,(uint64_t)(numSortedEvts & 0xFFFFFFFFFFFF));
         }
         numSortedEvts &= 0xFFFFFFFFFFFF;
         fseek(smolfp,0,SEEK_SET);
         fwrite(&numSortedEvts,sizeof(uint64_t),1,smolfp);
         printf("Wrote %lu separated events to output file.\n",numSortedEvts);
         fclose(smolfp);
      } else {
         printf("Can't open SMOL tree: %s to write\n", tmp);
         return(0);
      }

      sort->shutdown_midas = 1;
      pthread_join(midas_thread, NULL);
      pthread_join(ordthrd, NULL);
      pthread_join(ordthr2, NULL);
      pthread_join(grif_thread, NULL);
   }
   end=time(NULL);
   cfg->midas_start_time = diagnostics.midas_run_start;
   cfg->midas_runtime    = diagnostics.midas_last_timestamp+1;
   cfg->midas_runtime   -= cfg->midas_start_time;
   memcpy(cfg->midas_title, midas_runtitle, SYS_PATH_LENGTH);
   printf("File took %ld seconds\n", end-start);
   return(0);
}

extern volatile long grifevent_wrpos;
volatile long grifevent_rdpos;
extern Grif_event grif_event[MAX_COINC_EVENTS];

extern volatile unsigned long bankbuf_wrpos;
extern volatile unsigned long bankbuf_rdpos;
extern volatile long tsevents_in;
extern long tsevents_out;
extern volatile unsigned long eventbuf_rdpos;
extern volatile unsigned long eventbuf_wrpos;
extern volatile long grif_evcount;
void show_sort_state()
{
   int val = sort_status.midas_bytes/1000;
   int v2 = bankbuf_wrpos/1000000, v3 = bankbuf_rdpos/1000000;
   int v4 = tsevents_in/1000, v5 = tsevents_out/1000;
   int v6 = eventbuf_wrpos/1000000, v7 = eventbuf_rdpos/1000000;
   int v8 = grifevent_wrpos/1000, v9 = grifevent_rdpos/1000;
   int v10 = grifevent_wrpos - grifevent_rdpos;
   printf("MIDAS:read %d Mbytes [~%d Kevents]\n", val/1000, val/50);
   printf("      BUF in:%dMbytes out:%dMbytes  [Cap:%5.1f%%]\n",
          4*v2, 4*v3, (100000000.0*(v2-v3))/BANK_BUFSIZE );
   printf("REORDER: In:%dKevents Out%dKevents\n", v4, v5 );
   printf("     BUF In:%dMbytes  Out:%dMbytes  [Cap:%5.1f%%]\n",
          4*v6, 4*v7, (100000000.0*(v6-v7))/EVENTBUFSIZE );
   printf("GRIF: Unpacked:%dKevents  ", (int)(grif_evcount/1000) );
   printf("      Sorted  :%dKevents\n", done_events/1000 );
   printf("     BUF In:%dKevents Out%dKevents [=%d][Cap:%5.1f%%]\n\n",
          v8, v9, v10, (100.0*v10)/MAX_COINC_EVENTS);
}
Sort_status *get_sort_status(){ return( &sort_status ); }

uint64_t sort_main(Sort_status *arg, FILE *out)
{
   int i, len, nxtpos, rd_avail;
   static long grifevent_nxtpos;
   unsigned int usecs=10000;
   uint64_t numSortedEvts = 0U;

   printf("starting sort_main ...\n");
   grifevent_rdpos = grifevent_nxtpos = nxtpos = 0;
   while(1){
      // if( arg->shutdown_midas != 0 ){  break; }
      rd_avail = grifevent_wrpos - grifevent_nxtpos;
      if( arg->grif_sort_done && rd_avail < 1 ){ break; }
      if( rd_avail < 1 ){ usleep(usecs); continue; }
      numSortedEvts += process_event(&grif_event[nxtpos], nxtpos, out);
      //printf("sorted events: %lu\n",numSortedEvts);
      nxtpos = ++grifevent_nxtpos % MAX_COINC_EVENTS;
   }
   printf("sort_main finished\n");
   return numSortedEvts;
}

static int proc_calls, sorted, skipped, prefull, sortfull, completed_events;
// called when each new event read into ptr -> list[slot]
uint64_t process_event(Grif_event *ptr, int slot, FILE *out)
{
   time_t cur_time = time(NULL);
   static int calls, prv_call;
   static time_t prv_time;
   int dt = cur_time - prv_time, de = calls - prv_call;
   if( prv_time == 0 ){ prv_time = cur_time; }
   prv_call = ++calls;

   if( cur_time-prv_time >= 10 ){
      printf("----------------------------------------------------------------\n");
      midas_status(cur_time); reorder_status(cur_time);  grif_status(cur_time);
      printf("ProcEvt: %10d[Good:%3d%% Skip:%3d%% WinFull:%3d%%] %6.3f Mevt/s\n",
             calls, (int)(100.0*(calls-skipped-prefull)/calls),
             (int)(100.0*skipped/calls), (int)(100.0*prefull/calls), (de/(1000000.0*dt))
      );
      prv_time = cur_time;
   }
   if( singlethread_save == 1 && sort_status.odb_done == 0 ){ calls = 0;
      init_default_sort(configs[1], &sort_status); sort_status.odb_done = 1;
   }
   apply_gains(ptr);
   return insert_presort_win(ptr, slot, out);
}

char *debug_show_ts(long ts)
{
   static char tmp[32];
   int deci_ms = ts/10000;
   int remain = ts - (deci_ms*10000);

   sprintf(tmp,"%5.1fms+%04d", 00000.1*deci_ms, remain);
   return(tmp);
}
char *debug_show_chan(Grif_event *ptr)
{
   static char tmp[32];
   if( ptr->chan != -1 ){
      sprintf(tmp," %3d", ptr->chan);
   } else {
      sprintf(tmp,"%04x", ptr->address);
   }
   return(tmp);
}
extern char chan_name[MAX_DAQSIZE][CHAN_NAMELEN];

// add event to presort window (event has just been read in)
//    recalculate coincwin (sorting any events leaving window)
// => final win of run won't be sorted, as these events will not leave window
uint64_t insert_presort_win(Grif_event *ptr, int slot, FILE *out)
{
   int window_width = 500; // 5us to capture all pileup events - MAXIMUM (indiv. gates can be smaller)
   int win_count, win_end;
   Grif_event *alt;
   int64_t dt;
   uint64_t numSort = 0;

   /*
   static int prv_evt[65535], count;
   ++count;
   i = (slot-1-window_start+2*MAX_COINC_EVENTS) % MAX_COINC_EVENTS;
   printf("PRST:Chan[%4s:prev:%5d][%s] [win:%5d[%05d-%05d:%s]          ",
          debug_show_chan(ptr),
          prv_evt[ptr->address] == 0 ? -1 : count-prv_evt[ptr->address],
          debug_show_ts(ptr->ts),
          i, window_start, slot-1,
          debug_show_ts(grif_event[window_start].ts) );
   prv_evt[ptr->address] = count;
   if( ptr->chan != -1 ){
      printf("%s E=%6d[cal:%8.1f] id:0x%08x\n",
             chan_name[ptr->chan],ptr->energy,ptr->esum, ptr->master_id);
   } else {
      printf("---------- E=%6d[cal:%8.1f] id:0x%08x\n",
             ptr->energy, ptr->esum, ptr->master_id );
   }
   */

   ///////////////// Presort window (used for suppression/addback)
   while( presort_window_start != slot ){ 
      alt = &grif_event[presort_window_start];
      win_count = (slot - presort_window_start+2*MAX_COINC_EVENTS) % MAX_COINC_EVENTS;
      dt = ptr->ts - alt->ts; if( dt < 0 ){ dt *= -1; }

      // should exit while-loop when no more events outside window
      //    *BUT* add error recovery - if window too full, dump events
      if( dt < window_width ){
         //if( win_count < 0.45*MAX_COINC_EVENTS  ){ break; }
         if( win_count < coinc_events_cutoff ){ break; } // LIMIT to ?? events
         else { ++prefull; }
      }

      // event[win_start] is leaving window
      //    ( either because dt > coincwidth OR due to error recovery)
      // NOTE event[slot] is out of window - use slot-1 as window-end
      if( (win_end = slot-1) < 0 ){ win_end = MAX_COINC_EVENTS-1; } // WRAP
      pre_sort(presort_window_start, win_end);
      numSort += insert_sort_win(alt, presort_window_start, out); // add event to next window
      if( ++presort_window_start >= MAX_COINC_EVENTS ){ presort_window_start=0; } // WRAP
   }
   return numSort;
}

// add event to main sort window (event has just left presort window)
uint64_t insert_sort_win(Grif_event *ptr, int slot, FILE *out)
{
   int window_width = 200; // 2us - MAXIMUM (indiv. gates can be smaller)
   int win_count, win_end;
   Grif_event *alt;
   int64_t dt;
   uint64_t numSort = 0;

   /* i = (slot-1-window_start+2*MAX_COINC_EVENTS) % MAX_COINC_EVENTS;
   printf("MAIN:Chan[%4s:    :     ][%s] [win:%5d[%05d-%05d:%s]\n",
          debug_show_chan(ptr),

          debug_show_ts(ptr->ts),
          i, window_start, slot-1,
          debug_show_ts(grif_event[window_start].ts) );
   */
   while( sort_window_start != slot ){
      alt = &grif_event[sort_window_start];
      win_count = (slot - sort_window_start+2*MAX_COINC_EVENTS) % MAX_COINC_EVENTS;
      dt = ptr->ts - alt->ts; if( dt < 0 ){ dt *= -1; }

      // should exit while-loop when no more events outside window
      //    *BUT* add error recovery - if window too full, dump events
      if( dt < window_width ){
       //if( win_count >= 0.45*MAX_COINC_EVENTS ){ ++sortfull; } else {
         if( win_count > coinc_events_cutoff ){ ++sortfull; } else {
            break;
         }
      }

      // event[win_start] is leaving window
      //    ( either because dt > coincwidth OR due to error recovery)
      // NOTE event[slot] is out of window - use slot-1 as window-end
      if( (win_end = slot-1) < 0 ){ win_end = MAX_COINC_EVENTS-1; } // WRAP
      if( alt->chan != -1 ){ ++sorted;
         numSort += (uint64_t)fill_smol_entry(out, sort_window_start, win_end);
      } else { ++skipped; }
      if( ++sort_window_start >= MAX_COINC_EVENTS ){ sort_window_start=0; } // WRAP
      ++grifevent_rdpos;  ++completed_events;
   }
   //printf("insert_sort_win - sorted %u event(s)\n",numSort);
   return numSort;
}

/////////////////////////////////////////////////////////////////////////////////
///////// Alternate sort method [build complete events, then sort them] /////////
/////////   [faster?] but slightly worse at finding all coincidences    /////////
/////////////////////////////////////////////////////////////////////////////////


// add event to window (event has just been read in)
//    (previously would loop over events already in window to see if they are now OUT)
// NEW method - if first event is OUT ...
//    build an event with all fragments before current fragment (which starts a new event)
/*uint64_t build_event(Grif_event *ptr, int slot, FILE *out)
{
   int window_width = 200; // 2us - MAXIMUM (indiv. gates can be smaller)
   int win_count, win_end;
   static int window_start;
   Grif_event *alt;
   int64_t dt;
   uint64_t numSort = 0;

   if( window_start == slot ){ return(0); }

   alt = &grif_event[window_start];
   win_count = (slot - window_start+2*MAX_COINC_EVENTS) % MAX_COINC_EVENTS;
   dt = ptr->ts - alt->ts; if( dt < 0 ){ dt *= -1; }
   if( dt < window_width ){ return(0); }

   if( (win_end = slot-1) < 0 ){ win_end = MAX_COINC_EVENTS-1; } // WRAP
   numSort += sort_built_event(window_start, win_end, out);

   window_start = slot;
   return numSort;
}

uint64_t sort_built_event(int window_start, int win_end, FILE *out)
{
   uint64_t numSort = 0;
   pre_sort(window_start, win_end); // only fold of first fragment is set
   numSort += fill_smol_entry(out, window_start, win_end);
   return numSort;
}*/
