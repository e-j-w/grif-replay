// 0) make sure user sort runs
//
// 1) reset bytes sorted
//     ** Runs are entries in filelist[+may contain multiple subruns] **
//     so handle subruns differently
//     status - subrun is x/y or 1/1 ****
//
// DONE 2) strings replaced with %20 - globals etc
//
// 3) getHistofileList - Histogram files need run title info [+date/time]
//                     [Already have odbstuff in config file]
//    ALSO remember first and last midas-timestamps
//      add config section: SortedFileDetails
//
// DONE 4) Add radware matrix script [optionally all histos]
//
// 5) Add 180-deg-coinc matrix - for summing calcns
//       56Co - 4400keV   1keV/chan
//
//    Add 2-photon decay 1dhisto of sum of pairs of crystals [2photonGe]
//       Add 2-photon Variable as well
//
// *** 6) sort gzipped files
//
// DONE 7) check root script with config file entry
//
// DONE 8) Reorder histo list - Energy Time Waveform Pulsehight
//
// DONE 9) Fix Gains from config not midas
//
// odb-error
//
// Sr90[]
//
/////////////////////////////////////////////////////////////////////////////

/* new js fn similar to odbget - will have to decode args in server.c */
/* Test commands ...

wget 'http://grifstore1:9093/?cmd=addDatafile&filename=/tig/grifstore1b/grifalt/schedule145/Dec2023/run21834_000.mid'

wget 'http://grifstore1:9093/?cmd=getDatafileList&dir=/tig/grifstore1b/grifalt/schedule145/Dec2023'

wget 'http://localhost:9093/?cmd=getDatafileList&dir=.'
wget 'http://localhost:9093/?cmd=getHistofileList&dir=.'
wget 'http://localhost:9093/?cmd=getSortStatus'
wget 'http://localhost:9093/?cmd=addDatafile&filename=/tig/grifstore1b/grifalt/schedule145/Dec2023/run21834_000.mid'
wget 'http://localhost:9093/?cmd=openHistofile&filename=histos.tar'


wget 'http://localhost:9093/?cmd=addGate&filename=histos.cfg'
wget 'http://localhost:9093/?cmd=saveConfig&filename=histos.cfg'



wget 'http://panther:9093/?cmd=addDatafile&filename=/tig/grifstore0b/griffin/schedule140/Calibrations-July2021/run18909_020.mid'

wget 'http://panther:9093/?cmd=endCurrentFile'

wget 'http://panther:9093/?cmd=getSpectrumList&filename=run18909_020.tar'

wget 'http://panther:9093/?cmd=callspechandler&spectum0=Hitpattern_Energy&spectrum1=GRG01BN00A_Energy'

*/
//////////////////////////////////////////////////////////////////////////////
//"spectrum1d/index.html"
/*    check for javascript cmds or send file               */
/*    file should be under custom pages - check url for .. */
// Format of CALLSPECHANDLER call is:
// host:PORT/?cmd=callSpectrumHandler&spectum0=firstSpec&spectrum1=secondSpec
// CURRENTLY DEFINED COMMANDS ...
//     /?cmd=getDatafileList&dir=XXXX
//     /?cmd=getHistofileList
//     /?cmd=getSortStatus
//     /?cmd=getSpectrumList
//     /?cmd=callSpectrumHandler
//     /?cmd=addDatafile&filename=XXXX
//
// ALSO FIXED URLs ...
//     /filter-status.html
//     /report
//     /*.css
//     /*.js
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <ctype.h>
#include "config.h"
#include "midas2smol.h"

///////////////////////////////////////////////////////////////////////////
//////////////         Url  Command  interpreter          /////////////////
///////////////////////////////////////////////////////////////////////////

extern int coinc_events_cutoff;

///////////////////////////////////////////////////////////////////////////
/////////////////            Config  I/O            ///////////////////////
///////////////////////////////////////////////////////////////////////////
Config *configs[MAX_CONFIGS];

int write_config(Config *cfg, FILE *fp)
{
   Cal_coeff *calib;  Global *global;
   Sortvar *var;      Cond *cond;        Gate *gate;
   int i, j;
   char tmp[64];

   fprintf(fp,"{\n   \"Analyzer\" : [\n");
   fprintf(fp,"      {\"Variables\" : [\n");
   for(i=0; i<cfg->nsortvar; i++){ var = &cfg->varlist[i];
      sprintf(tmp, "\"%s\"", var->name );
      fprintf(fp,"%9s{ \"name\" : %-16s,", "", tmp );
      fprintf(fp,"   \"title\" : \"%s\"}", var->title  );
      fprintf(fp,"%s", ( i<cfg->nsortvar-1 ) ? ",\n" : "\n" );
   }
   fprintf(fp,"      ]},\n");
   fprintf(fp,"      {\"Globals\" : [\n");
   for(i=0; i<cfg->nglobal; i++){ global = cfg->globals[i];
      fprintf(fp,"%9s{\"name\" : \"%s\" , \"min\" : %d , \"max\" : %d ",
              "", global->name, global->min, global->max );
      fprintf(fp, "%s", ( i<cfg->nglobal-1 ) ? "},\n" : "}\n" );
   }
   fprintf(fp,"      ]},\n");
   fprintf(fp,"      {\"Calibrations\" : [\n");
   for(i=0; i<cfg->ncal; i++){ calib = cfg->calib[i];
      fprintf(fp,"%9s{\"name\" : \"%s\" , \"address\" : %d , \"datatype\" : %d , \"offset\" : %f , \"gain\" : %f , \"quad\" : %e ", "", calib->name, calib->address, calib->datatype, calib->offset, calib->gain, calib->quad );
      if(calib->pileupk1[0] != -1 && !isnan(calib->pileupk1[0])){
        fprintf(fp,", \"pileupk1\" : [ %f , %f , %e , %e , %e , %e , %e ]",calib->pileupk1[0],calib->pileupk1[1],calib->pileupk1[2],calib->pileupk1[3],calib->pileupk1[4],calib->pileupk1[5],calib->pileupk1[6]);
        fprintf(fp,", \"pileupk2\" : [ %f , %f , %e , %e , %e , %e , %e ]",calib->pileupk2[0],calib->pileupk2[1],calib->pileupk2[2],calib->pileupk2[3],calib->pileupk2[4],calib->pileupk2[5],calib->pileupk2[6]);
        fprintf(fp,", \"pileupE1\" : [ %f , %f , %e , %e , %e , %e , %e ]",calib->pileupE1[0],calib->pileupE1[1],calib->pileupE1[2],calib->pileupE1[3],calib->pileupE1[4],calib->pileupE1[5],calib->pileupE1[6]);
      }else{
        fprintf(fp,", \"pileupk1\" : [ %d , %d , %d , %d , %d , %d , %d ]",1,0,0,0,0,0,0);
        fprintf(fp,", \"pileupk2\" : [ %d , %d , %d , %d , %d , %d , %d ]",1,0,0,0,0,0,0);
        fprintf(fp,", \"pileupE1\" : [ %d , %d , %d , %d , %d , %d , %d ]",0,0,0,0,0,0,0);
      }
      fprintf(fp, "%s", ( i<cfg->ncal-1 ) ? "},\n" : "}\n" );
   }
   fprintf(fp,"      ]},\n");
   fprintf(fp,"      {\"Directories\" : [\n");
   {
      fprintf(fp,"%9s{\"name\" : \"Data\", ", "");
      fprintf(fp,"\"Path\" : \"%s\"},\n", cfg->data_dir);
      fprintf(fp,"%9s{\"name\" : \"Config\", ", "");
      fprintf(fp,"\"Path\" : \"%s\"}\n", cfg->config_dir); // NO COMMA
   }
   fprintf(fp,"      ]},\n");
   fprintf(fp,"      {\"Midas\" : [\n");
   {
      fprintf(fp,"%9s{\"name\" : \"Title\", ", "");
      fprintf(fp,"\"Value\" : \"%s\"},\n", cfg->midas_title);
      fprintf(fp,"%9s{\"name\" : \"StartTime\", ", "");
      fprintf(fp,"\"Value\" : \"%d\"},\n", cfg->midas_start_time);
      fprintf(fp,"%9s{\"name\" : \"Duration\", ", "");
      fprintf(fp,"\"Value\" : \"%d\"}\n", cfg->midas_runtime);  // NO COMMA
   }
   fprintf(fp,"      ]}\n"); // NO COMMA AFTER FINAL TABLE
   fprintf(fp,"   ]\n}\n"); // Analyser,File
   return(0);
}

static char config_data[1024*1024];
int load_config(Config *cfg, char *filename, char *buffer)
{
   int i,j, len, value, val2, val3, val4, val5, val6, address, type, instring;
   char *ptr, *name, *valstr, *title, *path, *var, *var2, op[8], tmp[80];
   float gain, offset, quad;
   float puk1[7], puk2[7], puE1[7];
   Config *tmp_cfg;
   Cond *cond;
   FILE *fp;

   if( filename != NULL ){ len = strlen(filename);
      if( (fp=fopen(filename,"r")) == NULL ){
         fprintf(stderr,"load_config: cant open %s to read\n", filename);
         return(-1);
      }
      fprintf(stderr,"load_config: reading %s\n", filename);
      instring = len = 0;// read line by line, copy to conf_data, skip space
      while( fgets(tmp, 80, fp) != NULL ){ // tmp always null-terminated
         for(i=0; i<strlen(tmp); i++){ // DO NOT SKIP SPACE WITHIN STRINGS
            if( tmp[i] == '"' ){ instring = 1-instring; }
            if( !isspace(tmp[i]) || instring ){ config_data[len++] = tmp[i]; }
         }
      }
      fclose(fp);
   } else if( buffer != NULL ){
      instring = len = 0; // read line by line, copy to conf_data, skip space
      for(i=0; i<strlen(buffer); i++){    // DO NOT SKIP SPACE WITHIN STRINGS
         if( buffer[i] == '"' ){ instring = 1-instring; }
         if( !isspace(buffer[i])||instring ){ config_data[len++]=buffer[i]; }
      }
   } else {
      fprintf(stderr,"load_config: no file or buffer specified\n");
   }
   clear_config(cfg);
   ptr=config_data;
   if( strncmp(ptr,"{\"Analyzer\":[", 13) != 0 ){
      fprintf(stderr,"load_config: err1 byte %ld\n", ptr-config_data);
       return(-1);
   } ptr += 13;
   if( strncmp(ptr,"{\"Variables\":[", 14) != 0 ){
      fprintf(stderr,"load_config: err2 byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 14;
   while( 1 ){ // variables
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: err3 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"title\":\"", 10) != 0 ){
         fprintf(stderr,"load_config: err4 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 10;
      title = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      // VARIABLES CURRENTLY HARDCODED
      //  - this section of config is only used for passing to viewer
      //cfg->lock=1; add_variable(cfg, name, title);  cfg->lock=0;
      ++ptr; // skip '}'
      if( *ptr++ == ',' ){ continue; }
      ptr+=2; break; // skip '},'
   }
   if( strncmp(ptr,"{\"Globals\":[", 12) != 0 ){
      fprintf(stderr,"load_config: errL byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 12;
   while( 1 ){ // Globals
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errM byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"min\":",7) != 0 ){
         fprintf(stderr,"load_config: errN byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 7;
      valstr = ptr; while( isdigit(*ptr) || *ptr=='-' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%d", &value) < 1 ){
         fprintf(stderr,"load_config:errO byte %ld\n", ptr-config_data);
         return(-1);
      }
      if( strncmp(ptr,"\"max\":",6) != 0 ){
         fprintf(stderr,"load_config: errP byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 6;
      valstr = ptr; while( isdigit(*ptr) || *ptr=='-' ){++ptr;} *ptr++ = 0;
      // last char written over with 0 was '}'
      if( sscanf( valstr, "%d", &val2) < 1 ){
         fprintf(stderr,"load_config:errQ byte %ld\n", ptr-config_data);
         return(-1);
      }
      cfg->lock=1; add_global(cfg, name, value, val2); cfg->lock=0;
      if( *ptr++ == ',' ){ continue; } // have skipped ']' if not
      ptr+=2; break; // skip '},'
   }
   if( strncmp(ptr,"{\"Calibrations\":[", 17) != 0 ){
      fprintf(stderr,"load_config: errR byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 17;
   while( 1 ){ // Calibrations
      if( strncmp(ptr,"]},", 3) == 0 ){ 
         ptr+=3; 
         //fprintf(stdout,"Calibrations section empty so breaking here.\n"); 
         break; }// empty section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errS byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"address\":",11) != 0 ){
         fprintf(stderr,"load_config: errV byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 11; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%d", &address) < 1 ){
         fprintf(stderr,"load_config:errX byte %ld\n", ptr-config_data);
         return(-1);
      }
      if( strncmp(ptr,"\"datatype\":",11) != 0 ){
         fprintf(stderr,"load_config: errV byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 11; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%d", &type) < 1 ){
        fprintf(stderr,"load_config:errX byte %ld\n", ptr-config_data);
        return(-1);
      }
      if( strncmp(ptr,"\"offset\":",9) != 0 ){
        fprintf(stderr,"load_config: errT byte %ld\n", ptr-config_data);
        return(-1);
      } ptr += 9; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%f", &offset) < 1 ){
        fprintf(stderr,"load_config:errU byte %ld\n", ptr-config_data);
        return(-1);
      }
      if( strncmp(ptr,"\"gain\":",7) != 0 ){
        fprintf(stderr,"load_config: errV byte %ld\n", ptr-config_data);
        return(-1);
      } ptr += 7; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%f", &gain) < 1 ){
        fprintf(stderr,"load_config:errX byte %ld\n", ptr-config_data);
        return(-1);
      }
      if( strncmp(ptr,"\"quad\":",7) != 0 ){
        fprintf(stderr,"load_config: errY byte %ld\n", ptr-config_data);
        return(-1);
      } ptr += 7; valstr = ptr;
      while(isdigit(*ptr)||*ptr=='.'||*ptr=='-'||*ptr=='+'||*ptr=='e'){++ptr;}
      *ptr++ = 0; // last char written over with 0 was '}'
      if( sscanf( valstr, "%f", &quad) < 1 ){
        fprintf(stderr,"load_config:errZ byte %ld\n", ptr-config_data);
        return(-1);
      }
      // The pileup correction parameters were introduced in Feb 2025.
      // Config files before this date will not have pileup correcitons, and after this they are optional
      if( strncmp(ptr,"\"pileupk1\":",11) == 0 ){
        ptr += 11; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%e,%e,%e,%e,%e], ", &puk1[0],&puk1[1],&puk1[2],&puk1[3],&puk1[4],&puk1[5],&puk1[6]) != 7 ){
          fprintf(stderr,"load_config:errPUA byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!='\"'){++ptr;}
        if( strncmp(ptr,"\"pileupk2\":",11) != 0 ){
          fprintf(stderr,"load_config:errPUB byte %ld\n", ptr-config_data);
          return(-1);
        } ptr += 11; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%e,%e,%e,%e,%e], ", &puk2[0],&puk2[1],&puk2[2],&puk2[3],&puk2[4],&puk2[5],&puk2[6]) != 7 ){
          fprintf(stderr,"load_config:errPUC byte %ld\n", ptr-config_data);
          return(-1);
        }
        //while(isdigit(*ptr)||*ptr=='.'||*ptr=='-'||*ptr=='+'||*ptr=='e'||*ptr==','||*ptr=='['||*ptr==']'){++ptr;}
        while(*ptr!='\"'){++ptr;}
        if( strncmp(ptr,"\"pileupE1\":",11) != 0 ){
          fprintf(stderr,"load_config:errPUD byte %ld\n", ptr-config_data);
          return(-1);
        } ptr += 11; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%e,%e,%e,%e,%e]", &puE1[0],&puE1[1],&puE1[2],&puE1[3],&puE1[4],&puE1[5],&puE1[6]) != 7 ){
          fprintf(stderr,"load_config:errPUE byte %ld\n", ptr-config_data);
          return(-1);
        }
        //while(isdigit(*ptr)||*ptr=='.'||*ptr=='-'||*ptr=='+'||*ptr=='e'||*ptr==','||*ptr=='['||*ptr==']'){++ptr;}
        while(*ptr!='}'){++ptr;}
        ++ptr;
      }else{
           // Initialize pileup parameters to default values
           for(i=0; i<7; i++){
             puk1[i] = puk2[i] = puE1[i] = 0;
           }
           puk1[0] = puk2[0] = 1; // set default factor as 1 not zero
      }
      cfg->lock=1; edit_calibration(cfg,name,offset,gain,quad,puk1,puk2,puE1,address,type,1); cfg->lock=0;
      if( *ptr++ == ',' ){ continue; } // have skipped ']' if not
      ptr+=2; break; // skip '},'
    }
   if( strncmp(ptr,"{\"Directories\":[", 16) != 0 ){
      fprintf(stderr,"load_config: errZA byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 16;
   while(1){
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errZB byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"Path\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errZC byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      path = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      cfg->lock=1; set_directory(cfg, name, path);  cfg->lock=0;
      ++ptr; // skip '}'
      if( *ptr++ == ',' ){ continue; }
      ptr+=2; break; // skip '},'
   }
   if( strncmp(ptr,"{\"Midas\":[", 10) != 0 ){
      fprintf(stderr,"load_config: errZD byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 10;
   while(1){
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errZE byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"Value\":\"", 10) != 0 ){
         fprintf(stderr,"load_config: errZF byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 10;
      valstr = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      cfg->lock=1; set_midas_param(cfg, name, valstr);  cfg->lock=0;
      ++ptr; // skip '}'
      if( *ptr++ == ',' ){ continue; }
      ++ptr; break; // skip '}'
   }
   if( strncmp(ptr,"]}", 2) != 0 || ptr+2-config_data != len ){
      fprintf(stderr,"load_config: errR near %ld\n", ptr-config_data);
      return(-1);
   }
   return(0);
}

///////////////////////////////////////////////////////////////////////////
/////////////////          Config Commands          ///////////////////////
///////////////////////////////////////////////////////////////////////////
// config changes during sorting ... (sort probably running most of time)
// - sort will need to grab config and copy so any changes do not affect it
//


// init_config() call on server startup ...
//   get most recent config and initialise all variables, gates, histos
// *In case no recent config exists, a default config is setup first*
//  (this is usually immediately overwritten by the recent config)
//  (could change this to only setup the default, if recent fails)
int init_config()
{
  int tmp;
   Config *cfg = configs[0];
   if( cfg == NULL ){ // not yet alloc'd live set
      if( (cfg=configs[0]=add_config("live")) == NULL ){ return(-1); }
      if( (configs[1]=add_config("sort")) == NULL ){ return(-1); }
      configs[0]->type = configs[1]->type = MEM_CONFIG;
   }
   init_default_config(cfg);  // populate default "test" config during testing
   load_config(cfg, DEFAULT_CONFIG, NULL); // attempt to load, ignore any error
   clear_calibrations(cfg); // Clear the calibrations to default values following server restart
   fprintf(stdout,"Initial setup complete :-)\n\n");
   return(0);
}

int clear_config(Config *cfg)
{
   time_t current_time = time(NULL);
   int i;
   memset(cfg, 0, sizeof(Config));      // delete any current vars etc.
   for(i=0; i<MAX_CALIB;     i++){cfg->calib[i]    = &cfg->calib_array[i]; }
   for(i=0; i<MAX_GLOBALS;   i++){cfg->globals[i]  = &cfg->global_array[i]; }
   for(i=0; i<MAX_CONDS;     i++){cfg->condlist[i] = &cfg->cond_array[i]; }
   for(i=0; i<MAX_GATES;     i++){cfg->gatelist[i] = &cfg->gate_array[i]; }
   // used_sortvars, and user histos start empty and are filled randomly
   cfg->mtime = current_time;
   return(0);
}

int clear_calibrations(Config *cfg)
{
  float offset=0, gain=1, quad=0;
  float puk1[7], puk2[7], puE1[7];
  int i, address=-1, datatype=-1;

  // Initialize values to defaults
  for(i=0; i<7; i++){
    puk1[i] = puk2[i] = puE1[i] = 0;
  }
  puk1[0] = puk2[0] = 1; // set default factor as 1 not zero

  // delete any calibration values
  for(i=0; i<cfg->ncal;     i++){
    edit_calibration(cfg, cfg->calib[i]->name, offset, gain, quad, puk1, puk2, puE1, address, datatype, 1);
  }
  return(0);
}

int copy_config(Config *src, Config *dst)
{
   Global *global;
   Sortvar *srcvar, *dstvar;
   char *tmp, *tmp2, *ptr;
   long *tmp3;
   int i, j;

   src->lock = 1;
   // below is wrong - src array lists can contain holes if things were deleted
   //    use same offsets as in src arrays:
   //          offset = src->calib[i] - &src->calib_array[0];
   //   dst->calib[i] = offset + &dst->calib_array[0]
   for(i=0; i<MAX_CALIB;     i++){dst->calib[i]    = &dst->calib_array[i]; }
   for(i=0; i<MAX_GLOBALS;   i++){dst-> globals[i] = &dst->global_array[i]; }

   src->lock = 0; return(0);
}

Config *add_config(char *name)
{
   Config *cfg;
   int i, len;

   for(i=0; i<MAX_CONFIGS; i++){
      if( configs[i] == NULL ){ break; }
   }
   if( i == MAX_CONFIGS ){
      fprintf(stderr,"Exceeded Max histogram sets\n");
   }
   if( (cfg = configs[i] = calloc(1, sizeof(Config))) == NULL ){
      fprintf(stderr,"can't alloc new config\n");
   }
   if( (len=strlen(name)+1) >  SYS_PATH_LENGTH ){
      fprintf(stderr,"truncating configname: %s\n", name);
      len = SYS_PATH_LENGTH;
   }
   memcpy(cfg->name, name, len);
   for(i=0; i<MAX_CALIB;   i++){cfg->calib[i]    = &cfg->calib_array[i]; }
   for(i=0; i<MAX_GLOBALS; i++){cfg->globals[i] = &cfg->global_array[i]; }
   for(i=0; i<MAX_CONDS;   i++){cfg->condlist[i] = &cfg->cond_array[i]; }
   for(i=0; i<MAX_GATES;  i++){cfg->gatelist[i] = &cfg->gate_array[i]; }
   return( cfg );
}

int remove_config(Config *cfg)
{
   int i;

   if( cfg == configs[0] || cfg == configs[0] ){ return(0); }// do not remove set#0 or #1
   for(i=1; i<MAX_CONFIGS; i++){
      if( cfg == configs[i] ){ break; }
   }
   if( i == MAX_CONFIGS ){
      fprintf(stderr,"can't find config to remove\n");
   }
   clear_config(cfg);
   free(cfg);
   configs[i] = NULL;
   return(0);
}

int save_config(Config *cfg, char *filename, int overwrite)
{
   FILE *fp;
   if( cfg->lock ){ return(-1); } // config file currently being read
   if( !overwrite ){
      if( (fp=fopen(filename,"r")) != NULL ){
         fprintf(stderr,"save_config: file %s already exists\n", filename);
         fclose(fp); return(-1);
      }
   }
   if( (fp=fopen(filename,"w")) == NULL ){
      fprintf(stderr,"save_config: cant open %s to write\n", filename);
      return(-1);
   }
   write_config(cfg, fp);
   fclose(fp);
   return(0);
}

/////////////////////////////   Variable   /////////////////////////////////

// THERE IS CURRENTLY NO WAY TO CALCULATE THE VALUE OF A VARIABLE
//   GIVEN ONLY ITS TEXT DESCRIPTION
//      SO THIS FUNCTION WOULD BE POINTLESS AT THE MOMENT

//int add_variable(Config *cfg, char *name, char *title)
//{
//   time_t current_time = time(NULL);
//   int i = cfg->nsortvar;
//   if( i == MAX_SORT_VARS ){
//      fprintf(stderr,"too many variables when adding %s\n", name); return(-1);
//   }
//   memcpy(cfg->varlist[i].name, name, strlen(name)+1);
//   memcpy(cfg->varlist[i].title, title, strlen(title)+1);
//   ++cfg->nsortvar;
//   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
//  return(0);
//}

Sortvar *find_sortvar(Config *cfg, char *name)
{
   int i;
   for(i=0; i<cfg->nsortvar; i++){
      if( strcmp(cfg->varlist[i].name, name) == 0 &&
          strlen(cfg->varlist[i].name) == strlen(name) ){
         return( &cfg->varlist[i] );
      }
   }
   return(NULL);
}

/////////////////////////  Global Condition   //////////////////////////////

// every array member has a pointer to itself
//    [the pointer list entries are swapped etc, but not cleared]
// same with gates, conditions, variables
// but *not* used-variables, user-histos
int add_global(Config *cfg, char *name, int value, int val2)
{
   time_t current_time = time(NULL);
   Global *global;
   int i, len;
   for(i=0; i<cfg->nglobal; i++){
      if( strcmp(cfg->globals[i]->name, name) == 0 &&
          strlen(cfg->globals[i]->name) == strlen(name) ){ break; }
   }
   global = cfg->globals[i];
   if( i == cfg->nglobal ){ // new global - i is first unused ptr
      if( cfg->nglobal >= MAX_GLOBALS ){
         fprintf(stderr,"too many globals for %s\n", name); return(-1);
      }
      if( (len=strlen(name)+1) > STRING_LEN ){
         fprintf(stderr,"truncating globalname: %s\n", name);
         len = STRING_LEN;
      }
      memcpy(global->name, name, len);
       ++cfg->nglobal;
   }
   global->min = value; global->max = val2;
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

int remove_global(Config *cfg, char *name)
{
   Global *global, *lastglobal = cfg->globals[cfg->nglobal];
   time_t current_time = time(NULL);
   int i;
   for(i=0; i<cfg->nglobal; i++){
      if( strcmp(cfg->globals[i]->name, name) == 0 &&
          strlen(cfg->globals[i]->name) == strlen(name) ){
         global = cfg->globals[i]; break;
      }
   }
   if( i == cfg->nglobal ){
      fprintf(stderr,"can't find global: %s to remove\n", name); return(-1);
   }
//   if( global->use_count != 0 ){
//      fprintf(stderr,"global[%s] still in use[%d]\n", name, global->use_count);
//      return(-1);
//   }
   // if removing final entry - no rearrangement needed
   if( i != cfg->nglobal-1 ){  // otherwise - swap pointers with last
      cfg->globals[i             ] = lastglobal;
      cfg->globals[cfg->nglobal-1] = global;
   }
   --cfg->nglobal;
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

/////////////////////////////   CALIBRATION   /////////////////////////////

int set_calibration(Config *cfg, int num, char url_args[][STRING_LEN], int fd)
{
   float offset, gain, quad;
   float puk1[7], puk2[7], puE1[7];
   int i, address=-1, datatype=-1;
   char tmp[128];

   // Initialize values to -1
   for(i=0; i<7; i++){
     puk1[i] = puk2[i] = puE1[i] = -1;
   }

   for(i=2; i<num; i+=8){
      if( strncmp(url_args[i], "channelName", 10) != 0 ){
         sprintf(tmp,"set_calibration: expected \"channelName\" at %s\n",url_args[i]);
         
         fprintf(stderr,"expected \"channelName\" at %s\n", url_args[i]);
         return(-1);
      }
      if( strncmp(url_args[i+2], "quad", 4) != 0 ){
         sprintf(tmp,"set_calibration: expected \"quad\" at %s\n",url_args[i+2]);
         
         fprintf(stderr,"expected \"quad\" at %s\n", url_args[i+2]);
         return(-1);
      }
      if( sscanf(url_args[i+3], "%f", &quad) < 1 ){
         sprintf(tmp,"set_calibration: can't read quad value, %s\n",url_args[i+3]);
         
         fprintf(stderr,"can't read quad: %s\n", url_args[i+3]);
         return(-1);
      }
      if( strncmp(url_args[i+4], "gain", 4) != 0 ){
         sprintf(tmp,"set_calibration: expected \"gain\" at %s\n",url_args[i+4]);
         
         fprintf(stderr,"expected \"gain\" at %s\n", url_args[i+4]);
         return(-1);
      }
      if( sscanf(url_args[i+5], "%f", &gain) < 1 ){
         sprintf(tmp,"set_calibration: can't read gain value, %s\n",url_args[i+5]);
         
         fprintf(stderr,"can't read gain: %s\n", url_args[i+5]);
         return(-1);
      }
      if( strncmp(url_args[i+6], "offset", 6) != 0 ){
         sprintf(tmp,"set_calibration: expected \"offset\" at %s\n",url_args[i+6]);
         
         fprintf(stderr,"expected \"offset\" at %s\n", url_args[i+6]);
         return(-1);
      }
      if( sscanf(url_args[i+7], "%f", &offset) < 1 ){
         sprintf(tmp,"set_calibration: can't read offset value, %s\n",url_args[i+7]);
         
         fprintf(stderr,"can't read offset: %s\n", url_args[i+7]);
         return(-1);
      }
//      if( strncmp(url_args[i+8], "address", 4) != 0 ){
//         fprintf(stderr,"expected \"address\" at %s\n", url_args[i+8]);
//         return(-1);
//      }
//      if( sscanf(url_args[i+9], "%d", &address) < 1 ){
//         fprintf(stderr,"can't read address: %s\n", url_args[i+9]);
//         return(-1);
//      }
//      if( strncmp(url_args[i+10], "datatype", 6) != 0 ){
//         fprintf(stderr,"expected \"datatype\" at %s\n", url_args[i+10]);
//         return(-1);
//      }
//      if( sscanf(url_args[i+11], "%d", &datatype) < 1 ){
//         fprintf(stderr,"can't read datatype: %s\n", url_args[i+11]);
//         return(-1);
//      }
      edit_calibration(cfg, url_args[i+1], offset, gain, quad, puk1, puk2, puE1,
                       address, datatype, 1);
   }
   return(0);
}

int set_pileup_correction(Config *cfg, int num, char url_args[][STRING_LEN], int fd)
{
   float offset=-1, gain=-1, quad=-1;
   float puk1[7], puk2[7], puE1[7];
   int i, address=-1, datatype=-1;
   char tmp[128];

// Initialize values to defaults
for(i=0; i<7; i++){
  puk1[i] = puk2[i] = puE1[i] = 0;
}
puk1[0] = puk2[0] = 1; // set default factor as 1 not zero

   for(i=2; i<num; i+=8){
      if( strncmp(url_args[i], "channelName", 10) != 0 ){
         sprintf(tmp,"set_pileup_correction: expected \"channelName\" at %s\n",url_args[i]);
         
         fprintf(stderr,"expected \"channelName\" at %s\n", url_args[i]);
         return(-1);
      }
      if( strncmp(url_args[i+2], "pileupk1", 8) != 0 ){
         sprintf(tmp,"set_pileup_correction: expected \"pileupk1\" at %s\n",url_args[i+2]);
         
         fprintf(stderr,"expected \"pileupk1\" at %s\n", url_args[i+2]);
         return(-1);
      }
      if( sscanf(url_args[i+3], "%f,%f,%f,%f,%f,%f,%f", puk1,puk1+1,puk1+2,puk1+3,puk1+4,puk1+5,puk1+6) != 7 ){
         sprintf(tmp,"set_pileup_correction: can't read pileup k1 value (expected 7 values), %s\n",url_args[i+3]);
         
         fprintf(stderr,"can't read pileup k1, expected 7 values: %s\n", url_args[i+3]);
         return(-1);
      }
      if( strncmp(url_args[i+4], "pileupk2", 8) != 0 ){
         sprintf(tmp,"set_pileup_correction: expected \"pileupk2\" at %s\n",url_args[i+4]);
         
         fprintf(stderr,"expected \"pileupk2\" at %s\n", url_args[i+4]);
         return(-1);
      }
      if( sscanf(url_args[i+5], "%f,%f,%f,%f,%f,%f,%f", puk2,puk2+1,puk2+2,puk2+3,puk2+4,puk2+5,puk2+6) != 7 ){
         sprintf(tmp,"set_pileup_correction: can't read pileup k2 value (expected 7 values), %s\n",url_args[i+5]);
         
         fprintf(stderr,"can't read pileup k2, expected 7 values: %s\n", url_args[i+5]);
         return(-1);
      }
      if( strncmp(url_args[i+6], "pileupE1", 8) != 0 ){
         sprintf(tmp,"set_pileup_correction: expected \"pileupE1\" at %s\n",url_args[i+6]);
         
         fprintf(stderr,"expected \"pileupE1\" at %s\n", url_args[i+6]);
         return(-1);
      }
      if( sscanf(url_args[i+7], "%f,%f,%f,%f,%f,%f,%f", puE1,puE1+1,puE1+2,puE1+3,puE1+4,puE1+5,puE1+6) != 7 ){
         sprintf(tmp,"set_pileup_correction: can't read pileup E1 value (expected 7 values), %s\n",url_args[i+7]);
         
         fprintf(stderr,"can't read pileup E1, expected 7 values: %s\n", url_args[i+7]);
         return(-1);
      }
//      if( strncmp(url_args[i+14], "address", 4) != 0 ){
//         fprintf(stderr,"expected \"address\" at %s\n", url_args[i+14]);
//         return(-1);
//      }
//      if( sscanf(url_args[i+15], "%d", &address) < 1 ){
//         fprintf(stderr,"can't read address: %s\n", url_args[i+15]);
//         return(-1);
//      }
//      if( strncmp(url_args[i+16], "datatype", 6) != 0 ){
//         fprintf(stderr,"expected \"datatype\" at %s\n", url_args[i+16]);
//         return(-1);
//      }
//      if( sscanf(url_args[i+17], "%d", &datatype) < 1 ){
//         fprintf(stderr,"can't read datatype: %s\n", url_args[i+17]);
//         return(-1);
//      }
    edit_calibration(cfg, url_args[i+1], offset, gain, quad, puk1, puk2, puE1, address, datatype, 1);
   }
   return(0);
}

int edit_calibration(Config *cfg, char *name, float offset, float gain, float quad, float puk1[7], float puk2[7], float puE1[7], int address, int type, int overwrite)
{
   time_t current_time = time(NULL);
   int i,j, len, arg;
   Cal_coeff *cal;

   for(i=0; i<cfg->ncal; i++){ cal = cfg->calib[i];
      if( strncmp(name, cfg->calib[i]->name, strlen(name)) == 0 &&
          strlen(name) == strlen(cfg->calib[i]->name)       ){ break; }
   }
   if( i < cfg->ncal ){ // calib already exists
      if( overwrite ){
      if( offset != -1 ){ cal->offset = offset; }
      if( gain   != -1 ){ cal->gain = gain; }
      if( quad   != -1 ){ cal->quad = quad; }
      if( puk1[0] != -1 ){ for(j=0; j<7; j++){cal->pileupk1[j] = puk1[j];} }
      if( puk2[0] != -1 ){ for(j=0; j<7; j++){cal->pileupk2[j] = puk2[j];} }
      if( puE1[0] != -1 ){ for(j=0; j<7; j++){cal->pileupE1[j] = puE1[j];} }
      }
      if( address != -1 ){
         cal->address = address;  cal->datatype = type;
      }
   } else { // not already there ... add new one
      if( cfg->ncal >= MAX_CALIB ){
         fprintf(stderr,"too many calibs: %d\n", MAX_CALIB); return(-1);
      }
      cal = cfg->calib[i]; ++cfg->ncal;
      if( (len=strlen(name)+1) > 64 ){
         fprintf(stderr,"truncating calibname: %s\n", name);
         len = 64;
      }
      memcpy(cal->name, name, len);
      if( offset != -1 ){ cal->offset = offset; }else{ cal->offset = 0; }
      if( gain   != -1 ){ cal->gain = gain; }else{ cal->gain = 1; }
      if( quad   != -1 ){ cal->quad = quad; }else{ cal->quad = 0; }
      if( puk1[0] != -1 ){
        for(j=0; j<7; j++){ cal->pileupk1[j] = puk1[j]; }
      }else{
        cal->pileupk1[0]=1; cal->pileupk1[1]=0; cal->pileupk1[2]=0; cal->pileupk1[3]=0; cal->pileupk1[4]=0; cal->pileupk1[5]=0; cal->pileupk1[6]=0;
      }
      if( puk2[0] != -1 ){
        for(j=0; j<7; j++){cal->pileupk2[j] = puk2[j]; }
      }else{
        cal->pileupk2[0]=1; cal->pileupk2[1]=0; cal->pileupk2[2]=0; cal->pileupk2[3]=0; cal->pileupk2[4]=0; cal->pileupk2[5]=0; cal->pileupk2[6]=0;
      }
      if( puE1[0] != -1 ){
        for(j=0; j<7; j++){cal->pileupE1[j] = puE1[j]; }
      }else{
        cal->pileupE1[0]=0; cal->pileupE1[1]=0; cal->pileupE1[2]=0; cal->pileupE1[3]=0; cal->pileupE1[4]=0; cal->pileupE1[5]=0; cal->pileupE1[6]=0;
      }
      cal->address = address;  cal->datatype = type;
    }
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

int set_directory(Config *cfg, char *name, char *path)
{
   time_t current_time = time(NULL);
   int len;


   if( (len=strlen(path)) >= SYS_PATH_LENGTH ){
      fprintf(stderr,"set_directory: path too long[%s]\n", path);
      return(-1);
   }
   if( strncmp(name, "Data",   4) == 0 ){
      memcpy(cfg->data_dir, path, len+1);
   } else if( strncmp(name, "Config", 6) == 0 ){
      memcpy(cfg->config_dir, path, len+1);
   } else {
      fprintf(stderr,"set_directory: Unknown directory:%s\n", name);
      return(-1);
   }
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

int set_midas_param(Config *cfg, char *name, char *value)
{
   time_t current_time = time(NULL);
   int len;

   if( (len=strlen(value)) >= SYS_PATH_LENGTH ){
      fprintf(stderr,"set_midas_param: value too long[%s]\n", value);
      return(-1);
   }
   if(        strncmp(name, "Title",   5) == 0 ){
      memcpy(cfg->midas_title, value, len+1);
   } else if( strncmp(name, "StartTime",  9) == 0 ){
      if( sscanf(value, "%d", &cfg->midas_start_time) < 1 ){
         fprintf(stderr,"set_midas_param: can't read starttime: %s\n", value);
      }
   } else if( strncmp(name, "Duration", 8) == 0 ){
      if( sscanf(value, "%d", &cfg->midas_runtime) < 1 ){
         fprintf(stderr,"set_midas_param: can't read runtime: %s\n", value);
      }
   } else {
      fprintf(stderr,"set_midas_param: Unknown param:%s\n", name);
      return(-1);
   }
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}


#define READ_LIMIT (5*1024*1024) // don't read more than this looking for odb
int read_datafile_info(Sort_status *sort, char *path)
{
   int i, expt=0, ppg=0, flt=0, done=0;
   char tmp[256], *ptr;
   FILE *fp;
   if( (fp=fopen(path,"r")) == NULL ){
      fprintf(stderr,"read_datafile_info: can't read file: %s\n", path);
      return(-1);
   }
   while( (ptr=fgets(tmp, 256, fp)) != NULL ){
      if( ftell(fp) > READ_LIMIT ){ break; }
      while( isspace(*ptr) ){ ++ptr; }
      if( strncmp(ptr,"</dir>",6) == 0 ){ expt = 0; }
      if( strncmp(ptr,"<dir name=\"Run parameters\">",27) == 0 ){ expt = 2; }
      // The Title/comment are in Run parameters
      //    Experiment/Status-items and Edit-on-start contain LINKS
      //if( strncmp(ptr,"<dir name=\"Experiment\">",23)     == 0 ){ expt = 1; }
      //if( strncmp(ptr,"<dir name=\"Edit on start\">", 26) == 0 ){ expt = 3; }
      if( expt && strncmp(ptr,"<key name=\"Run Title\"",21) == 0 ){
         ptr+=21; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[0][i++] = *ptr++; }
         sort->file_info[0][i]=0; done |= 1;
      }
      if( expt && strncmp(ptr,"<key name=\"Comment\"",19) == 0 ){
         ptr+=19; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[1][i++] = *ptr++; }
         sort->file_info[1][i]=0; done |= 2;
      }
      if( strncmp(ptr,"<dir name=\"PPG\">",16)            == 0 ){ ppg = 1; }
      if( ppg && strncmp(ptr,"<key name=\"Current\"", 19) == 0 ){
         ptr+=19; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[2][i++] = *ptr++; }
         sort->file_info[2][i]=0; ppg = 0; done |= 4;
      }
      if( strncmp(ptr,"<dir name=\"Filter\">",19)         == 0 ){ flt = 1; }
      if( flt && strncmp(ptr,"<key name=\"Current\"", 19) == 0 ){
         ptr+=19; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[3][i++] = *ptr++; }
         sort->file_info[3][i]=0; flt = 0; done |= 8;
      }
      if( done == 15 ){ break; } // 1st 4 bits of done all set
   }
   fclose(fp);
   return(0);
}

int send_file_details(char *path, int fd)
{
   char tmp2[256];
   Sort_status *tmp;
   if( (tmp = calloc(sizeof(Sort_status), 1)) == NULL){
      
      fprintf(stderr,"send_file_details: failed alloc\n");
   }
   
   read_datafile_info(tmp, path);
   sprintf(tmp2,"%s\n%s\n%s\n%s\n", tmp->file_info[0], tmp->file_info[1],
                                   tmp->file_info[2], tmp->file_info[3] );
   free(tmp);
   return(0);
}

static struct stat statbuf;

// split path into name,dir  then get histo and config names
// also do stat to get size
int run_number(Sort_status *arg, char *name);
char *subrun_filename(Sort_status *sort, int subrun);
int add_sortfile(char *path)
{
   int i, plen, dlen, ext_len, hlen, clen;
   char ptr, tmp[256], *fname;
   Sort_status *sort;
   Config *cfg = configs[0];

   sort = get_sort_status();
   plen=strlen(path);
   ext_len = ( strncmp(path+plen-4, ".mid", 4) == 0 ) ? 4 : 0;
   for(i=plen; i>=0; i--){ if( path[i] == '/' ){ ++i; break; } }
   if( (dlen = i) == -1 ){ dlen = 0; } // no directory separator in path
   if( (sort->data_dir = malloc((size_t)(dlen + 2))) == NULL ){
      fprintf(stderr,"can't alloc string for data_dir");
      free_sortfile(sort); return(-1);
   }
   if( dlen == 0 ){ sprintf(sort->data_dir, ".");
   } else {
      memcpy((char *)sort->data_dir, path, dlen-1);
      *(sort->data_dir+dlen-1) = 0; // overwrite trailing '/'
      set_directory(cfg, "Data", sort->data_dir);
   }
   if( (sort->data_name = malloc(plen-dlen+1)) == NULL ){
      fprintf(stderr,"can't alloc string for :%s", path);
      free_sortfile(sort); return(-1);
   }
   memcpy((char *)sort->data_name, path+dlen, plen-dlen);
   *(sort->data_name+plen-dlen) = 0;
   if( run_number(sort, sort->data_name) ){ return(-1); }
   memset(&statbuf, 0, sizeof(struct stat)); sort->data_size = 0;

   //get number of subruns
   sort->subrun = -1;
   for(i=0; ; i++){
      fname = subrun_filename(sort, i);
         if( stat(fname, &statbuf) != 0 ){
            //fprintf(stderr,"can't stat multi-subrun: %s\n", path);
            //fprintf(stderr,"fname: %s\n", fname);
            break;
         }
         sort->data_size += (long)statbuf.st_size;
   }
   sort->num_subruns = i; sort->subrun = 0;

   return(0);
}

//////////////////////  sorting loop /////////////////////
// if( currentfile == finalfile ){ sleep(1); continue; }
// open_sortfilelist();
// sortnextfile()
//     create midas thread
//     read gains + init histos
//     sort data
//     close thread
//     write histos
// close_sortfilelist();
//////////////////////////////////////////////////////////
int open_next_sortfiles(Sort_status *sort)
{
   char ptr, tmp[256];
   if( sort->num_subruns == 0 ){
      sprintf(tmp, "%s", sort->data_name);
   } else {
      sprintf(tmp, "%s", subrun_filename(sort, 0) );
   }
   if( (sort->data_fp=fopen(tmp,"r")) == NULL ){
      fprintf(stderr,"can't open %s to read\n", tmp);  return(-1);
   }
   fprintf(stdout,"sorting file: %s\n", tmp);
   sort->midas_bytes = 0;
   sort->cal_overwrite = 1;
   return(0);
}

// Only first subrun contains odb event
// *could* just sort single subrun if subrun of specified data file is nonzero
//         otherwise sort all subruns
int open_next_subrun(Sort_status *sort)
{
   char *filename;
   fclose(sort->data_fp); sort->data_fp = NULL;

   if( sort->subrun >= sort->num_subruns-1 ){ // last or no subruns
      fprintf(stderr,"Final Subrun[%d] completed\n", sort->subrun); return(-1);
   }
   filename = subrun_filename(sort, ++sort->subrun);
   if( (sort->data_fp=fopen(filename,"r")) == NULL ){
      fprintf(stderr,"can't open %s to read\n", filename);  return(-1);
   }
   fprintf(stdout,"\n================================================\n");
   fprintf(stdout,"=============  SORTING SUBRUN %3d  =============\n", sort->subrun);
   fprintf(stdout,"================================================\n\n");
   return(0);
}

char *subrun_filename(Sort_status *sort, int subrun)
{
   static char name[256];
   int len, digits;
   char tmp[64];

   //printf("run: %i, run digits: %i, subrun digits: %i\n",sort->run, sort->run_digits, sort->subrun_digits);

   sprintf(name, "%s/run", sort->data_dir);

   sprintf(tmp,"%d", sort->run); digits = strlen(tmp);
   len = strlen(name);
   while( digits++ < sort->run_digits ){ name[len] = '0'; name[1+len++] = 0; }
   sprintf(name+strlen(name),"%d_", sort->run);

   sprintf(tmp,"%d", subrun);  digits = strlen(tmp);
   len = strlen(name);
   while( digits++ < sort->subrun_digits ){
      name[len++]='0'; name[len] = 0;
   }
   sprintf(name+strlen(name),"%d.mid", subrun);

   return(name);
}

int close_sortfiles(Sort_status *sort)
{  // data_fp usually closed in nxtsubrun
   if( sort->online_mode ){ return(0); } // no files in this mode
   if( sort->data_fp != NULL ){ fclose(sort->data_fp); }\
   free_sortfile(sort);
   return(0);
}

int free_sortfile(Sort_status *sort)
{
   if( sort->data_dir   != NULL ){ free(sort->data_dir);   }
   if( sort->data_name  != NULL ){ free(sort->data_name);  }
   memset(sort->data_dir, 0, sizeof(sort->data_dir));
   memset(sort->data_name, 0, sizeof(sort->data_name));
   return(0);
}

// read sun/subrun from filename, then count #digits in each
int run_number(Sort_status *arg, char *name)
{
   char *ptr = name, fmt[16], tmp[256];
   FILE *fp;
   int i;

   if( strncmp(ptr, "run", 3) != 0 ){
      fprintf(stderr,"datafilename:%s does not being with \"run\"\n", name);
      return(-1);
   } ptr += 3;
   while( 1 ){
      if( !isdigit(*ptr) ){ arg->run_digits = ptr-name-3; break; }
      ++ptr;
   }
   if( *ptr != 0 ){               // name contains stuff after run number ...
      if( *ptr == '_' ){
         if( sscanf(name, "run%d_%d.mid", &arg->run, &arg->subrun) != 2 ){
            fprintf(stderr,"can't read run and subrun number in %s\n", name);
            return(-1);
         }
         arg->subrun_digits = -1; ++ptr; while( 1 ){
            if( !isdigit(*ptr) ){
               arg->subrun_digits = ptr-name-4-arg->run_digits;  break;
            }
            ++ptr;
         }
         if( arg->subrun_digits == -1 ){
            fprintf(stderr,"can't read subrun number in %s\n", name);
            return(-1);
         }
         if( strncmp(ptr, ".mid", 4) == 0 ){ return(0); }
         else {
            fprintf(stderr,"no .mid extension in datafile: %s\n", name);
            return(-1);
         }
      } else if( strncmp(ptr, ".mid", 4) == 0 ){
         fprintf(stderr,"subrun number missing in datafile: %s\n", name);
         return(-1);
      } else {
         fprintf(stderr,"bad data filename format in %s\n", name);
         return(-1);
      }
   } else {               // name only contains run number - sort all subruns
      arg->subrun = -1;
      if( sscanf(name, "run%d.mid", &arg->run) != 1 ){
         fprintf(stderr,"can't read run number in %s\n", name);
         return(-1);
      }
      for(i=1; i<5; i++){ // look for up to 5 subrun digits (usually 3)
         sprintf(fmt, "%%s/%%s_%%0%dd.mid", i);
         sprintf(tmp, fmt, arg->data_dir, name, 0 );
         if( (fp = fopen(tmp,"r")) == NULL ){ continue; }
         arg->subrun_digits = i; fclose(fp); return(0);
      }
      fprintf(stderr,"can't open subrun0 for datafile:%s\n", name);
      return(-1);
   }
   return(0);
}

int end_current_sortfile(int fd)
{
   Sort_status *arg;

   arg = get_sort_status();
   arg->end_of_data = 1; //arg-> shutdown_midas = 1;
   return(0);
}


///////////////////////////////////////////////////////////////////////////
/////////////////          Directory reading          /////////////////////
///////////////////////////////////////////////////////////////////////////
#include <dirent.h>

int send_datafile_list(char *path, int fd, int type)
{
   char tmp[256];  Sort_status *tmp_srt;
   int nlen, run, subrun, entry=0;
   struct dirent *d_ent;
   DIR *d;

   if( (d=opendir(path)) == NULL ){
      sprintf(tmp,"can't open directory %s\n",path);
      fprintf(stderr,"can't open directory %s\n", path);
      return(-1);
   }
   set_directory(configs[0], "Data", path);
   if( type == 1 ){
      if( (tmp_srt  = calloc(sizeof(Sort_status),    1)) == NULL ){
         fprintf(stderr,"send_datafile_list: failed alloc\n");
      }
   } else { tmp_srt = NULL; }
   
   while( (d_ent = readdir(d)) != NULL ){
      //fprintf(stdout,"File[%s] ...\n", d_ent->d_name);
      if( strncmp(d_ent->d_name, ".", 1) == 0 ){
         continue; // Ignore
      }
      if( sscanf(d_ent->d_name, "run%d_%d.mid", &run, &subrun) != 2 ){
         if( sscanf(d_ent->d_name, "run%d.mid", &run) != 1 ){
            continue; // Not Midas Data File
         }
         subrun=0;
      }
      nlen = strlen(d_ent->d_name);
      if( strncmp(d_ent->d_name+nlen-4, ".mid", 4)     != 0 ){ // &&
          //strncmp(d_ent->d_name+nlen-7, ".mid.gz", 7)  != 0 &&
          //strncmp(d_ent->d_name+nlen-8, ".mid.bz2", 8) != 0 ){
         continue; // Not Midas DataFilename Extension
      }
      sprintf(tmp,"%s/%s", path, d_ent->d_name);
      if( stat(tmp, &statbuf) != 0 ){
         fprintf(stderr,"can't stat %s\n", tmp); statbuf.st_size = 1;
      }
      //put_line(fd, d_ent->d_name, strlen(d_ent->d_name) );
      sprintf(tmp," , %ld ", (long)statbuf.st_size);
      //put_line(fd, tmp, strlen(tmp) );
      if( (entry % 1000) == 0 ){ printf("Entry: %d\n", entry); }
      if( type == 0 ){ continue; }

      if( subrun == 0 ){
         sprintf(tmp,"%s/%s", path, d_ent->d_name);
         read_datafile_info(tmp_srt, tmp);
      } else { tmp_srt->file_info[0][0] = tmp_srt->file_info[1][0] = 0; }
      if( strlen(tmp_srt->file_info[0]) > 0 ){
         sprintf(tmp," , %s ", tmp_srt->file_info[0] );
      } else {
         sprintf(tmp," , %s ", tmp_srt->file_info[1] );
      }
      //put_line(fd, tmp, strlen(tmp) );
   }
   //put_line(fd, " ]\n", 3 );
   if( tmp_srt != NULL ){ free(tmp_srt); }
   return(0);
}

