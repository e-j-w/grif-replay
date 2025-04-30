#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "config.h"
#include "grif-format.h"
#include "smol-format.h"
#include "grif-angles.h"
#include "default_sort.h"

int          odb_daqsize;// number of daq channels currently defined in the odb
int         subsys_table[MAX_DAQSIZE];
int        crystal_table[MAX_DAQSIZE]; // Ge/BGO have 4 "crystals" per clover
int        element_table[MAX_DAQSIZE]; // BGO have 5 elements per crystal
int       polarity_table[MAX_DAQSIZE]; // 1 is negative, 0 is positive, -1 is unset
short       address_chan[MAX_ADDRESS];
static short  addr_table[MAX_DAQSIZE]; short   *addrs = addr_table;
       char    chan_name[MAX_DAQSIZE][CHAN_NAMELEN];
static int   dtype_table[MAX_DAQSIZE]; int    *dtypes = dtype_table;
static float  gain_table[MAX_DAQSIZE]; float   *gains = gain_table;
static float  offs_table[MAX_DAQSIZE]; float *offsets = offs_table;
static float  quad_table[MAX_DAQSIZE]; float   *quads = quad_table;
float  pileupk1[MAX_DAQSIZE][7];
float  pileupk2[MAX_DAQSIZE][7];
float  pileupE1[MAX_DAQSIZE][7];
static short *chan_address = addr_table;
static int subsys_initialized[MAX_SUBSYS];
extern Grif_event grif_event[MAX_COINC_EVENTS];

// Default sort function declarations
extern uint8_t fill_smol_entry(FILE *out, const int win_idx, const int frag_idx);

//generates a random double value on the interval [0,1]
//for smoothing purposes
double randomDbl(){
  return ((double)(rand()) / (double)(RAND_MAX));
}

//get the CFD corrected time
double getGrifTime(Grif_event *ptr){
  long cfdCorrTs = (ptr->ts & 0xFFFFFFFFFFFC0000); //timestamp, excluding lower 18 bits
  return ((double)cfdCorrTs + ((double)ptr->cfd/16.0) + randomDbl())*10.0;
  //return ((double)ptr->ts + randomDbl())*10.0;
}

double getGrifEnergy(Grif_event *ptr){
  double correctedE = (ptr->energy*1.0 + randomDbl());
  return ((double)offsets[ptr->chan]) + correctedE*(((double)gains[ptr->chan]) + correctedE*quads[ptr->chan]);
}

// odb tables need to be transferred into config, which is saved with histos
int init_default_sort(Config *cfg, Sort_status *arg)
{
  Cal_coeff *cal;
  int i, j;

  // Initialize all pileup parameters to unset values
  for(i=0; i<odb_daqsize; i++){
    for(j=0; j<7; j++){
      pileupk1[i][j] = pileupk2[i][j] = pileupE1[i][j] = -1;
    }
  }

  cfg->odb_daqsize = odb_daqsize;
  for(i=0; i<odb_daqsize; i++){ // update config with odb info
    edit_calibration(cfg, chan_name[i], offsets[i], gains[i], quads[i], pileupk1[i], pileupk2[i], pileupE1[i],
      chan_address[i],  dtype_table[i], arg->cal_overwrite );
  }
  // ALSO need to transfer config info to the arrays that are used in sort
  for(i=0; i<odb_daqsize; i++){

    cal = cfg->calib[i];
    if( strcmp(chan_name[i], cal->name) != 0 ){ // conf not in odb order
      for(j=0; j<cfg->ncal; j++){ cal = cfg->calib[j];
        if( strcmp(chan_name[i], cal->name) == 0 ){ break; }
      }
      if( j == cfg->ncal ){ continue; } // not found in config
    }

    // overwrite = 0 => USE CONFIG NOT ODB for offset, gain, quads
    if( arg->cal_overwrite == 0 ){
      offsets[i]=cal->offset; gains[i]=cal->gain;  quads[i]=cal->quad;
    }

    // Pileup parameters do not exist in the MIDAS ODB so must always be copied from the config
    for(j=0; j<7; j++){
      pileupk1[i][j] = (isnan(cal->pileupk1[j])) ? 0.0 : cal->pileupk1[j];
      pileupk2[i][j] = (isnan(cal->pileupk2[j])) ? 0.0 : cal->pileupk2[j];
      pileupE1[i][j] = (isnan(cal->pileupE1[j])) ? 0.0 : cal->pileupE1[j];
    }
  }

  return(0);
}

//#######################################################################
//#####        BASIC DEFAULT SORT (common to most experiments)      #####
//#######################################################################

float spread(int val){ return( val + rand()/(1.0*RAND_MAX) ); }
int GetIDfromAddress(unsigned short addr){ // address must be an unsigned short
  return(address_chan[addr]);
}

// this is the first function to be called on processing an event -
// before any presort/singles/coinc-sorting
int apply_gains(Grif_event *ptr)
{
  float energy;
  int chan = ptr->chan;

  // Protect against invalid channel numbers
  if( chan < 0 || chan >= odb_daqsize ){
    fprintf(stderr,"unpack_event: ignored event in chan:%d [0x%04x]\n", chan, ptr->address );
    return(-1);
  }

  // Calculate the energy and calibrated energies
  ptr->energy = energy = ( ptr->integ == 0 ) ? ptr->q : spread(ptr->q)/ptr->integ;
  ptr->ecal=ptr->esum = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);

  // NOBODY CURRENTLY USES e2,e3,e4 ...
  if( ptr->integ2 != 0 ){
    energy = ptr->energy2 = spread(ptr->q2)/ptr->integ2;
    ptr->e2cal = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
  }
  if( ptr->integ3 != 0 ){
    energy = ptr->energy3 = spread(ptr->q3)/ptr->integ3;
    ptr->e3cal = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
  }
  if( ptr->integ4 != 0 ){
    energy = ptr->energy4 = spread(ptr->q4)/ptr->integ4;
    ptr->e4cal = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
  }

  // Assign the subsys type
  if( (ptr->subsys = subsys_table[chan]) == -1 ){ return(-1); }

  // HPGe pileup
  if( ptr->subsys == SUBSYS_HPGE_A){
    ptr->psd = 14; // Pileup class - default value of 12 for all HPGe events
    if(ptr->pileup==1 && ptr->nhit ==1){
      // Single hit events
      // no pileup, this is the most common type of HPGe event
      ptr->psd = 1; // Pileup class, default for single hit events
    }
  }

  return(0);
}

// Presort - do Suppression and Addback here
//  - frag_idx is about to leave coinc window (which ends at end_idx)
//    check other frags in window for possible suppression and/or summing
//  also calculate multiplicities[store in frag_idx only]
int pre_sort(int frag_idx, int end_idx)
{
  Grif_event *alt2, *alt, *ptr = &grif_event[frag_idx];
  int i, j, dt, dt13;
  float q1,q2,q12,k1,k2,k12,e1,e2,e12,m,c;
  int chan,found,pos;
  float energy,correction;
  float correction12, correction23;

  // Assign chan local variable and check it is a valid channel number
  if( (chan=ptr->chan)<0 || ptr->chan >= odb_daqsize ){
    fprintf(stderr,"presort error: ignored event in chan:%d\n",ptr->chan );
    return(-1);
  }
  i = frag_idx; ptr->fold = 1;
  while( i != end_idx ){ // need at least two events in window
    if( ++i >=  MAX_COINC_EVENTS ){ i=0; } alt = &grif_event[i]; // WRAP
    if( alt->chan<0 || alt->chan >= odb_daqsize ){
      fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
      return(-1);
    }

    // Determine fold
    if( alt->subsys == ptr->subsys ){ ++ptr->fold; }

    // Determine absolute time difference between timestamps
    dt = ptr->ts - alt->ts; if( dt < 0 ){ dt = -1*dt; }

    // SubSystem-specific pre-processing
    switch(ptr->subsys){
      case SUBSYS_HPGE_A:

      // HPGe pile-up corrections
      // THE PRE_SORT WINDOW SHOULD BE EXTENDED TO COVER THE FULL POSSIBLE TIME DIFFERENCE BETWEEN PILE-UP events
      // THIS IS EQUAL TO THE DIFF PERIOD OF HPGE TYPE
      // First assign the pileup class type, then correct the energies
      if(alt->subsys == SUBSYS_HPGE_A && alt->chan == ptr->chan){
        if(ptr->pileup==1 && ptr->nhit ==1){
          // no pileup, this is the most common type of HPGe event
          ptr->psd = 1; // Pileup class
        }else if(ptr->pileup==0){
          // pileup error
          ptr->psd = 0; // Pileup class, error
        }else if((ptr->pileup==1 && ptr->nhit==2) && (alt->pileup==2 && alt->nhit==1)){
          ptr->psd = alt->psd = 9; // Pileup class, error for 2Hits
          ptr->ts_int = alt->ts_int = dt;
          if(ptr->q>0 && ptr->integ>0 && ptr->q2>0 && ptr->integ2>0 && alt->q>0 && alt->integ>0){
            // 2 Hit, Type A
            // The (ptr) fragement is the first Hit of a two Hit pile-up event.
            // Assign the pileup class numbers to the two hits and calculate the time difference between them
            ptr->psd = 3; // Pileup class, 1st of 2Hits
            alt->psd = 4; // Pileup class, 2nd of 2Hits
            ptr->ts_int = alt->ts_int = dt;

          }else{
            // 2 Hit, Type B
            // 2Hit pileup where 2nd Hit integration region starts after 1st Hit integration ends but before 2nd Hit CFD has completed
            // Assign the pileup class numbers to the two hits and calculate the time difference between them
            ptr->psd = 7; // Pileup class, 1st of 2Hits
            alt->psd = 8; // Pileup class, 2nd of 2Hits
            ptr->ts_int = alt->ts_int = dt;

          }
        }else if((ptr->pileup==1 && ptr->nhit==2) && (alt->pileup==1 && alt->nhit==1)){
          // 2 Hit, Type C
          // 2Hit pileup where 2nd Hit integration region starts after 1st Hit integration and 2nd Hit CFD have ended
          // Assign the pileup class numbers to the two hits and calculate the time difference between them
          ptr->psd = 5; // Pileup class, 1st of 2Hits
          alt->psd = 6; // Pileup class, 2nd of 2Hits
          ptr->ts_int = alt->ts_int = dt; // Save the time difference between pileup hits into both hits

        }else if((ptr->pileup==1 && ptr->nhit==3) && (alt->pileup==2 && alt->nhit==2)){ // 3Hit pileup
          ptr->psd = alt->psd = 13; // Pileup class, error for 3Hits
          if(ptr->q>0 && ptr->integ>0 && ptr->q2>0 && ptr->integ2>0 && alt->q>1 && alt->integ>0 && alt->q2>0 && alt->integ2>0){
            /*
            found=0;
            //  if(ptr->ecal > 1160 && ptr->ecal < 1180 && alt->ecal > 1325 && alt->ecal < 1350){
            fprintf(stdout,"Found a pileup 3 hit group for chan %d with dt %d\n",ptr->chan,dt);
            fprintf(stdout,"ptr: %ld %d, PU=%d, nhits=%d, q: %d %d %d %d, k: %d %d %d %d, ecal: %d %d %d %d, %lf %lf %lf\n",ptr->ts,ptr->cfd,ptr->pileup,ptr->nhit,ptr->q,ptr->q2,ptr->q3,ptr->q4,ptr->integ,ptr->integ2,ptr->integ3,ptr->integ4,ptr->ecal,ptr->e2cal,ptr->e3cal,ptr->e4cal,offsets[ptr->chan],gains[ptr->chan],quads[ptr->chan]);
            fprintf(stdout,"alt: %ld %d, PU=%d, nhits=%d, q: %d %d %d %d, k: %d %d %d %d, ecal: %d %d %d %d, %lf %lf %lf\n",alt->ts,alt->cfd,alt->pileup,alt->nhit,alt->q,alt->q2,alt->q3,alt->q4,alt->integ,alt->integ2,alt->integ3,alt->integ4,alt->ecal,alt->e2cal,alt->e3cal,alt->e4cal,offsets[alt->chan],gains[alt->chan],quads[alt->chan]);
            fprintf(stdout,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",ptr->q,ptr->q2,alt->q,ptr->integ,ptr->integ2,alt->integ,ptr->energy,ptr->energy2,alt->energy,ptr->ecal,ptr->e2cal,alt->ecal);
            found=1;
            //    }
            */
            j=i+1;
            while( j != end_idx ){ // need to find the third events in window associated with this channel
              if( ++j >=  MAX_COINC_EVENTS ){ break; } alt2 = &grif_event[j]; // WRAP

              if(alt2->chan == ptr->chan){ // It must also be a HPGe if the channel number is the same
                /*
                fprintf(stdout,"alt2: %ld %d, PU=%d, nhits=%d, q: %d %d %d %d, k: %d %d %d %d, ecal: %d %d %d %d, %lf %lf %lf\n",alt2->ts,alt2->cfd,alt2->pileup,alt2->nhit,alt2->q,alt2->q2,alt2->q3,alt2->q4,alt2->integ,alt2->integ2,alt2->integ3,alt2->integ4,alt2->ecal,alt2->e2cal,alt2->e3cal,alt2->e4cal,offsets[alt2->chan],gains[alt2->chan],quads[alt2->chan]);
                */
                if(alt2->pileup==3 && alt2->nhit==1){
                  alt2->psd = 12; // Pileup class
                  if(alt2->q>1 && alt2->integ>0){

                    // Determine absolute time difference between timestamps for Hit 1 and 3
                    dt13 = ptr->ts - alt2->ts; if( dt13 < 0 ){ dt13 = -1*dt13; }


                    // Three-Hit pile-up case...
                    // there are two types depending on the relative timing of the third Hit...
                    // Triple pileup case A ... in which the 3rd pulse occurs more than L samples after the first
                    //                          there are 5 regions, 3 of which are not piled up (1 per pulse)
                    //                      ____                ____
                    //    |   |        |  /|    |\  |      |  /|    |\  |      |          :
                    //    |   |        | / |    | \ |      | / |    | \ |      |          :
                    //    |   *________|/__|____|  \|______|/__|____|  \|______|          :
                    //    |  /                   \                  \          \          :
                    //    | /   K1           K12  \    K2        K23 \     K3   \         :
                    //  __*/                       \_____             \______    \________:
                    //    0            S                   X
                    //
                    // ------------------------------------------------------------------------------------------
                    // Triple pileup case B ... in which the 3rd pulse occurs less than L samples after the first
                    //                          again 5 regions, only 2 of which are not piled up (first and last pulse)
                    //                          There is no region to obtain the height of pulse 2
                    //                          so the event contains K12, the sum of pulse 1+2, in place of pulseheight2
                    //                                    ________
                    //    |   |        |   |         |  /|        |\  |      |   |      |   :
                    //    |   |        |   |         | / |        | \ |      |   |      |   :
                    //    |   |        |   |_________|/__|________|  \|______|   |      |   :
                    //    |   |        |  /|          :           |\  |      |\  |      |   :
                    //    |   |        | / |          :           | \ |      | \ |      |   :
                    //    |   *________|/__|__________:___________|  \|______|__\|______|   :
                    //    |  /                        :           \                     \   :
                    //    | /     K1            K12   :    K123    \     K23        K3   \  :
                    //  __*/                          :             \_____                \_:
                    //    0            S              X           L
                    //

                    // The Differencitation period of the HPGe Pulse Height evaluation is L = 5000ns.
                    if(dt13>500){
                      // Triple pileup case A ... in which the 3rd pulse occurs more than L samples after the first
                      //                          there are 5 regions, 3 of which are not piled up (1 per pulse)
                      correction23 = (alt->q/alt->integ)-((alt->q2/alt->integ2)-(alt2->q/alt2->integ));
                      correction12 = (ptr->q/ptr->integ)-((ptr->q2/ptr->integ2)-(alt->q/alt->integ)-correction23);
                      // Hit 1
                      ptr->psd = 10; // Pileup class
                      ptr->energy = energy = (spread(ptr->q)/ptr->integ) + correction12;
                      ptr->ecal=ptr->esum = offsets[ptr->chan]+energy*(gains[ptr->chan]+energy*quads[ptr->chan]);
                      // Hit 2
                      alt->ts_int = dt;
                      alt->psd = 11; // Pileup class
                      alt->energy = energy = (spread(alt->q)/alt->integ) - correction12 + correction23;
                      alt->ecal=alt->esum = offsets[alt->chan]+energy*(gains[alt->chan]+energy*quads[alt->chan]);
                      // Hit 3
                      alt2->ts_int = dt13;
                      alt2->psd = 12; // Pileup class
                      alt2->energy = energy = (spread(alt2->q)/alt2->integ) - correction23;
                      alt2->ecal=alt2->esum = offsets[alt2->chan]+energy*(gains[alt2->chan]+energy*quads[alt2->chan]);
                    }else{
                      // Triple pileup case B ... in which the 3rd pulse occurs less than L samples after the first
                      //                          again 5 regions, only 2 of which are not piled up (first and last pulse)
                      //                          There is no region to obtain the height of pulse 2
                      //                          so the event contains K12, the sum of pulse 1+2, in place of pulseheight2
                      correction23 = (alt->q/alt->integ)-((alt->q2/alt->integ2)-(alt2->q/alt2->integ));
                      correction12 = (ptr->q/ptr->integ)-((ptr->q2/ptr->integ2)-(alt->q/alt->integ)-correction23);
                      // Hit 1
                      ptr->psd = 10; // Pileup class
                      ptr->energy = energy = (spread(ptr->q)/ptr->integ) + correction12;
                      ptr->ecal=ptr->esum = offsets[ptr->chan]+ptr->energy*(gains[ptr->chan]+ptr->energy*quads[ptr->chan]);
                      // Hit 2
                      alt->ts_int = dt;
                      alt->psd = 11; // Pileup class
                      alt->energy = energy = (spread(alt->q)/alt->integ) - correction12 + correction23;
                      alt->ecal=alt->esum = offsets[alt->chan]+energy*(gains[alt->chan]+energy*quads[alt->chan]);
                      // Hit 3
                      alt2->ts_int = dt13;
                      alt2->psd = 12; // Pileup class
                      alt2->energy = energy = (spread(alt2->q)/alt2->integ) - correction23;
                      alt2->ecal=alt2->esum = offsets[alt2->chan]+energy*(gains[alt2->chan]+energy*quads[alt2->chan]);
                    }
                    /*
                    fprintf(stdout,"Complete 3Hit PU event, dt13=%d: %d,%d,%d,%d,%d, %d,%d,%d,%d,%d, %d,%d,%d,%d,%d, %d,%d,%d,%d,%d\n\n",dt13,ptr->q,ptr->q2,alt->q,alt->q2,alt2->q,ptr->integ,ptr->integ2,alt->integ,alt->integ2,alt2->integ,ptr->energy,ptr->energy2,alt->energy,alt->energy2,alt2->energy,ptr->ecal,ptr->e2cal,alt->ecal,alt->e2cal,alt2->ecal);
                    */
                    break; // Break the while if we found the third Hit
                  }
                }
              }
            } // end of while for triple coincidence
          }
        } // end of 3Hit pileup type assignments

        // Now apply hit-specific energy corrections
        if(pileupk1[chan][0] != 1){
          if(ptr->psd>=3 && ptr->psd<=8){ // 2-Hit pileup events
            // Apply the k1 dependant correction to the energy of the first hit
            // It was already checked that chan for ptr and alt are the same for pileup events
            pos  = crystal_table[ptr->chan];
            k1 = ptr->integ;
            ptr->ecal=ptr->esum = ptr->ecal*( pileupk1[chan][0]+(k1*pileupk1[chan][1])+(k1*k1*pileupk1[chan][2])+(k1*k1*k1*pileupk1[chan][3])
            +(k1*k1*k1*k1*pileupk1[chan][4])+(k1*k1*k1*k1*k1*pileupk1[chan][5])+(k1*k1*k1*k1*k1*k1*pileupk1[chan][6]));
            alt->e4cal=ptr->ecal; // Remember the ecal of the first Hit in this second Hit

            // Apply the E1 and k2 dependant offset correction to the energy of the second hit
            // Apply the k2 dependant correction to the energy of the second hit
            k2 = alt->integ;
            correction = ptr->ecal*( pileupE1[chan][0]+(k2*pileupE1[chan][1])+(k2*k2*pileupE1[chan][2])+(k2*k2*k2*pileupE1[chan][3])
            +(k2*k2*k2*k2*pileupE1[chan][4])+(k2*k2*k2*k2*k2*pileupE1[chan][5])+(k2*k2*k2*k2*k2*k2*pileupE1[chan][6]));
            alt->ecal=alt->esum = (alt->ecal*( pileupk2[chan][0]+(k2*pileupk2[chan][1])+(k2*k2*pileupk2[chan][2])+(k2*k2*k2*pileupk2[chan][3])
            +(k2*k2*k2*k2*pileupk2[chan][4])+(k2*k2*k2*k2*k2*pileupk2[chan][5])+(k2*k2*k2*k2*k2*k2*pileupk2[chan][6])))+correction;
          }
        }
      } // end of if(alt->subsys == SUBSYS_HPGE_A && alt->chan == ptr->chan)

      // BGO suppression of HPGe
      if( (dt >= bgo_window_min && dt <= bgo_window_max) && alt->subsys == SUBSYS_BGO && !ptr->suppress ){
        // could alternatively use crystal numbers rather than clover#
        //    (don't currently have this for BGO)
        if( crystal_table[ptr->chan]/16 == crystal_table[alt->chan]/16 ){ ptr->suppress = 1; }
      }
      break;
      default: break; // Unrecognized or unprocessed subsys type
    }// end of switch
  }// end of while
  return(0);
}

//writes data for a single sorted_evt in a SMOL tree
int lastWinIdx = -1;
uint8_t fill_smol_entry(FILE *out, const int win_idx, const int frag_idx)
{
  //fprintf(stdout,"Called fill entry for win: %i, frag: %i, last win: %i\n",win_idx,frag_idx,lastWinIdx);
  if(win_idx == 0){
    lastWinIdx = -1; //wrap
  }
  if(frag_idx < win_idx){
    lastWinIdx = -1; //wrap
    return 0;
  }
  if(win_idx <= lastWinIdx){
    return 0; //don't double fill
  }
  Grif_event *ptr;
  int i;

  // initialize SMOL tree event
  sorted_evt *sortedEvt = calloc(1,sizeof(sorted_evt));
  sortedEvt->header.numHPGeHits = 0;

  for(i=win_idx; ; i++){
    
    ptr = &grif_event[i];

    if( i >= MAX_COINC_EVENTS ){ i=0; } //wrap 
    //if( i != win_idx && flag == SORT_ONE ){ break; }
    if( ptr->dtype == 15 ){ if( i==frag_idx ){ break; } continue; } // scalar
    if( ptr->chan == -1 ){
        printf("DefSort: UNKNOWN_CHAN type=%d\n", ptr->dtype);
        if( i==frag_idx ){ break; } continue;
    } //  ????
    
    /*if( ptr->chan<0 || ptr->chan >= odb_daqsize ){
      //fprintf(stderr,"SmolSort: UNKNOWN_CHAN=%i type=%d\n",ptr->chan,ptr->dtype);
      if( i==frag_idx ){ break; } continue;
    }
    if( i != win_idx && i==frag_idx ){ break; }
    if( ptr->dtype == 15 ){ if( i==frag_idx ){ break; } continue; } // scalar*/
    //fprintf(stdout,"  checking index %i (ts=%li)\n",i,ptr->ts);
    
    switch(ptr->subsys){
      case SUBSYS_HPGE_A: // Ge
        // Only use GRGa
        if(ptr->suppress == 0){
          //passes Compton suppression
          int c1 = crystal_table[ptr->chan];
          if( c1 >= 0 && c1 < 64){
            if(sortedEvt->header.numHPGeHits >= MAX_EVT_HIT){
              fprintf(stderr,"WARNING: too many hits in win_idx %i, frag_idx %i\n",win_idx,frag_idx);
              break;
            }
            double grifT = getGrifTime(ptr);
            //fprintf(stdout,"ts: %li, time: %f ns\n",ptr->ts,grifT);
            if(sortedEvt->header.evtTimeNs == 0){
              sortedEvt->header.evtTimeNs = grifT;
            }
            sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].energy = (float)getGrifEnergy(ptr);
            sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].timeOffsetNs = (float)(grifT - sortedEvt->header.evtTimeNs);
            sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core = (uint8_t)(c1);
            if(sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core >= 64){
              fprintf(stderr,"WARNING: invalid GRIFFIN core: %u",sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core);
              break;
            }
            //filter out duplicate data
            uint8_t dupFound = 0;
            for(int i = 0; i<sortedEvt->header.numHPGeHits;i++){
              if(sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core == sortedEvt->hpgeHit[i].core){
                if(fabs(sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].timeOffsetNs - sortedEvt->hpgeHit[i].timeOffsetNs) < 20.0){
                  dupFound = 1; //duplicate hit found
                }
              }
            }
            if(dupFound == 0){
              //hit is not a duplicate, so store it
              sortedEvt->header.numHPGeHits++;
            }
          }
        }
        break; // outer-switch-case-GE
      case SUBSYS_BGO:
        //at least one suppressor fired
        sortedEvt->header.metadata |= (uint8_t)(1U << 1);
        break;
      default:
        break; // Unrecognized or unprocessed subsys type
    }// end of switch(ptr)

    lastWinIdx = i;
    if( i==frag_idx ){ break; }
  }
  
  if((sortedEvt->header.numHPGeHits > 0)&&(sortedEvt->header.numHPGeHits <= MAX_EVT_HIT)){

    //finalize sorted event data
    sortedEvt->header.metadata |= (uint8_t)(1U << 7); //set data validation bit

    //write sorted event to SMOL tree
    fwrite(&sortedEvt->header,sizeof(evt_header),1,out);
    //write hits
    for(uint8_t j = 0; j<sortedEvt->header.numHPGeHits;j++){
      fwrite(&sortedEvt->hpgeHit[j].timeOffsetNs,sizeof(float),1,out);
      fwrite(&sortedEvt->hpgeHit[j].energy,sizeof(float),1,out);
      fwrite(&sortedEvt->hpgeHit[j].core,sizeof(uint8_t),1,out);
      //fprintf(stdout,"Hit %u - core: %u, energy: %0.2f, time offset: %0.2f, win: %i, frag: %i\n",j,sortedEvt->hpgeHit[j].core,(double)sortedEvt->hpgeHit[j].energy,(double)sortedEvt->hpgeHit[j].timeOffsetNs,win_idx,frag_idx);
    }
    free(sortedEvt);
    return 1; //sorted an event
  }

  free(sortedEvt);
  return 0;
  
}

//#######################################################################
//###########   READ XML ODB DUMP FROM START OF DATA FILE   #############
//#######################################################################

// Note - the odb names do not distinguish between subtypes of detectors
// e.g Ge A and B channels
// the subsystem names will be extended to include this information
// (and the odb-specific names are only used below)

#define MAX_ODB_SUBSYS 24
#define ODBHANDLE_GRG   0
#define ODBHANDLE_GRS   1
#define ODBHANDLE_SEP   2
#define ODBHANDLE_PAC   3
#define ODBHANDLE_LBS   4
#define ODBHANDLE_LBT   5
#define ODBHANDLE_LBL   6
#define ODBHANDLE_DSC   7
#define ODBHANDLE_ART   8
#define ODBHANDLE_ZDS   9
#define ODBHANDLE_RCS  10
#define ODBHANDLE_XXX  11
#define ODBHANDLE_DSW  12
#define ODBHANDLE_DSG  13
#define ODBHANDLE_DAL  14
#define ODBHANDLE_DAT  15
#define ODBHANDLE_UNK  23
static char odb_handle[MAX_ODB_SUBSYS][8] = {
   "GRG", "GRS", "SEP", "PAC",  //  0- 3
   "LBS", "LBT", "LBL", "DSC",  //  4- 7
   "ART", "ZDS", "RCS", "XXX",  //  8-11
   "DSW", "DSG", "DAL", "DAT",  //  12-15
   "",    "",    "",    "",
   "",    "",    "",    "UNK"
};

static char   path[256];
static char dirname[64],value[32],type[32];
extern char midas_runtitle[SYS_PATH_LENGTH];

static void *arrayptr;
int read_odb_items(int len, int *bank_data)
{
   char *path_ptr, *ptr, *str, *odb_data = (char *)bank_data, posn[2];
   int i, c = '<', d = '>', dtype=0, active=0, index=0;
   ptr = odb_data;  path_ptr = path;
   while(1){
      if( (str = strchr(ptr,c)) == NULL ){ break; }
      ptr = str;
      if( (str = strchr(ptr,d)) == NULL ){ break; }

      if( strncmp(ptr,"<!--",4) == 0 || strncmp(ptr,"<odb", 4) == 0 ||
                                        strncmp(ptr,"</odb",5) == 0 ){ // comment - skip
      } else if( strncmp(ptr,"<dir ",5) == 0 ){
         if( strncmp(ptr,"<dir name=\"",11) == 0 ){
            i=11; while( *(ptr+i) != '"' && *(ptr+i) != d ){ ++i; }
         }
         memcpy(dirname, ptr+11, i-11); dirname[i-11] = '\0';
         if( *(ptr+1+i) == '/' ){ ptr=str+1; continue; }
         //if( sscanf(ptr,"<dir name=\"%s\">", dirname) < 1 ){
         //   fprintf(stderr,"can't read dirname\n"); ptr=str+1; continue;
         //}
         //if( strncmp(dirname+strlen(dirname)-3,"\"/>",3) == 0 ){
         //   ptr=str+1; continue;
         //}
         //if( dirname[strlen(dirname)-1]=='>'  ){
         //   dirname[strlen(dirname)-1]='\0';
         //}
         //if( dirname[strlen(dirname)-1]=='\"' ){
         //  dirname[strlen(dirname)-1]='\0';
         //}
         *path_ptr = '/'; strcpy(path_ptr+1, dirname);
         path_ptr += strlen(dirname)+1;
         *path_ptr = '\0';
      } else if( strncmp(ptr,"</dir>",6) == 0 ){
         while(1){
            if( --path_ptr < path ){ path_ptr = path;  *path_ptr = '\0';  break; }
            if( *path_ptr == '/' ){ *path_ptr = '\0';  break; }
         }
         index=0; // for debugger to stop here
      } else if( strncasecmp(ptr,"<key name=\"Run Title\" type=\"STRING\"", 35) == 0 ){
         ptr = str+1;
         if( (str = strchr(ptr,c)) == NULL ){ break; }
         i = (str-ptr) > SYS_PATH_LENGTH-1 ? SYS_PATH_LENGTH-1 : (str-ptr);
         memcpy( midas_runtitle, ptr, i ); midas_runtitle[i] = 0;
         ptr += i+1;
         if( (str = strchr(ptr,d)) == NULL ){ break; }
      } else if( strncmp(ptr,"</keyarray>",10) == 0 ){ active = 0; arrayptr = (void *)('\0');
      } else if( strncmp(ptr,"<keyarray ",10) == 0 ){
         if( strcmp(path,"/DAQ/params/MSC") != 0 &&
             strcmp(path,"/DAQ/MSC")        != 0 &&
             strcmp(path,"/DAQ/PSC")        != 0 ){  ptr=str+1; continue; }
         if( sscanf(ptr,"<keyarray name=\"%s", value) < 1 ){
            fprintf(stderr,"can't read keyarray entry\n"); ptr=str+1; continue;
         }
         if( value[strlen(value)-1]=='\"' ){ value[strlen(value)-1]='\0'; }
         if( strcmp(value,"PSC") == 0 || strcmp(value,"MSC") == 0 ){
            active = 1; arrayptr = (void *)addr_table; dtype=1;
         }
         if( strcmp(value,"chan") == 0 ){
            active = 1; arrayptr = (void *)chan_name; dtype=3;
         }
         if( strcmp(value,"datatype") == 0 ){
          //active = 1; arrayptr = (void *)dtype_table; dtype=1;
            active = 1; arrayptr = (void *)dtype_table; dtype=0;
         }
         if( strcmp(value,"gain") == 0 ){
            active = 1; arrayptr = (void *)gain_table; dtype=2;
         }
         if( strcmp(value,"offset") == 0 ){
            active = 1; arrayptr = (void *)offs_table; dtype=2;
         }
         if( strcmp(value,"quadratic") == 0 ){
            active = 1; arrayptr = (void *)quad_table; dtype=2;
         }
      } else if( strncmp(ptr,"<value index=",13) == 0 ){
         if( !active ){ ptr=str+1; continue; }
         // remove the >< surrounding the value, and move str to the end of the line
         *str = ' '; if( (str = strchr(str,c)) == NULL ){ break; }
         *str = ' '; if( (str = strchr(str,d)) == NULL ){ break; }
         if( sscanf(ptr,"<value index=\"%d\" %s /value>", &index, value) < 2 ){
            fprintf(stderr,"can't read value entry\n");
         }
         if( index < 0 || index >= MAX_DAQSIZE ){
            fprintf(stderr,"index %d out of range\n", index);
         }
         // index starts at zero, odb_daqsize is count
         if( index >= odb_daqsize ){ odb_daqsize = index+1; }
         if(        dtype == 0 ){  // int
            if( sscanf(value,"%d", (((int *)arrayptr)+index)) < 1 ){
               fprintf(stderr,"can't read value %s\n", value);
            }
         } else if( dtype == 1 ){  // short int
            if( sscanf(value,"%hd", (((short *)arrayptr)+index)) < 1 ){
               fprintf(stderr,"can't read value %s\n", value);
            }
         } else if( dtype == 2 ){  // float
            if( sscanf(value,"%f", (((float *)arrayptr)+index)) < 1 ){
               fprintf(stderr,"can't read value %s\n", value);
            }
         } else {                 // string
            strncpy(arrayptr+index*CHAN_NAMELEN, value, CHAN_NAMELEN);
            *((char *)arrayptr+(index+1)*CHAN_NAMELEN - 1) = '\0';
         }
      }
      ptr=str+1;
   }
   fprintf(stdout,"odb record: %d bytes\n", len);

   // arrays typically around 500 entries [one per "chan"] each entry with ...
   //   daq-address, name, type, gains etc.
   //
   gen_derived_odb_tables();

   return(0);
}

// original odb arrays were read into {addr_table,chan_name,dtype_table(+gains)}
// extract extra details stored in channel names (and record for later)
// (these details include crystal/element numbers and polarities)
// [use above for subsystem (no longer use datatype to determine subsystems)]
extern int read_caen_odb_addresses(int odb_daqsize, unsigned short *addr_table);
int gen_derived_odb_tables()
{
   int i, j, tmp, subsys, pos, element, output_type;
   char sys_name[64], crystal, polarity, type;

   read_caen_odb_addresses(odb_daqsize, (unsigned short *)addrs);

   // generate reverse mapping of address to channel number
   //  (most of this array is undefined and stays at -1)
   memset(address_chan, 0xFF, sizeof(address_chan)); // set to -1
   for(i=0; i<MAX_ADDRESS && i<odb_daqsize; i++){
      address_chan[ (unsigned short)chan_address[i] ] = i;
   }

   memset(crystal_table,  0xff, MAX_DAQSIZE*sizeof(int)); // initialise all to -1
   memset(element_table,  0xff, MAX_DAQSIZE*sizeof(int));
   memset(polarity_table, 0xff, MAX_DAQSIZE*sizeof(int));
   memset(subsys_table,   0xff, MAX_DAQSIZE*sizeof(int));
   for(i=0; i<MAX_DAQSIZE && i<odb_daqsize; i++){
      if( (tmp=sscanf(chan_name[i], "%3c%d%c%c%d%c", sys_name, &pos, &crystal, &polarity, &element, &type)) != 6 ){
         fprintf(stderr,"can't decode name[%s] decoded %d of 6 items\n", chan_name[i], tmp );
         continue;
      }
      for(j=0; j<MAX_ODB_SUBSYS; j++){
         if( strncmp(sys_name, odb_handle[j], 3) == 0 ){ subsys = j; break; }
      }
      if( j == MAX_ODB_SUBSYS ){ subsys = j-1; // use final entry: "unknown"
         fprintf(stderr,"Unknown subsystem[%s] in %s\n", sys_name, chan_name[i]);
      }

      // Mention bad detector types (no longer relied on for subsystem id)
      if( dtype_table[i] < 0 || dtype_table[i] >= 16 ){
         fprintf(stderr,"bad datatype[%d] at table position %d\n", dtype_table[i], i);
      }

      // Some detector elements have more than one output (HPGe A and B)
      // 1 is A, 0 is B, -1 is X or unknown
      output_type = type=='A' ? 1 : (type=='B' ? 0 : -1);

      // Polarity: 1 is N, 0 is P or T or S, -1 is anything else
      if(        polarity == 'N' ){ polarity_table[i] = 1;
      } else if( polarity == 'P' ){ polarity_table[i] = 0;
      } else if( polarity == 'T' ){ polarity_table[i] = 0; // TAC signal
      } else if( polarity == 'S' ){ polarity_table[i] = 1; // ARIES Standard Ouput signal
      } else if( polarity == 'F' ){ polarity_table[i] = 0; // ARIES Fast Output signal
      } else if( polarity == 'X' ){ polarity_table[i] = 0; // XXX type
      } else { fprintf(stderr,"unknown polarity[=%c] in %s\n", polarity, chan_name[i]); }

      // Record crystal and element numbers [** Naming schemes are subsystem-dependant **]
      switch(subsys){
      case ODBHANDLE_LBL: case ODBHANDLE_LBS: // LaBr,Paces, Ares and Zds
      case ODBHANDLE_LBT: case ODBHANDLE_ART:
      case ODBHANDLE_DAL: case ODBHANDLE_DAT:
      case ODBHANDLE_PAC: case ODBHANDLE_ZDS: case ODBHANDLE_DSW:
         crystal_table[i] = pos;
         if(        crystal == 'A' ){ element_table[i] = 1;
         } else if( crystal == 'B' ){ element_table[i] = 2;
         } else if( crystal == 'C' ){ element_table[i] = 3;
         } else if( crystal == 'X' ){ element_table[i] = -1; // just one crystal for LaBr3, ZDS, ART
         } else {
            fprintf(stderr,"unknown crystal for ancillary[=%c] in %s\n", crystal, chan_name[i]);
         } break;
      case ODBHANDLE_RCS:
         crystal_table[i] = pos;
         element_table[i] = element;
         break;
      case ODBHANDLE_GRG: case ODBHANDLE_GRS:
         element_table[i] = element;
         pos -= 1; pos *=4;
         if(        crystal == 'B' ){ crystal_table[i] = pos;
         } else if( crystal == 'G' ){ crystal_table[i] = pos+1;
         } else if( crystal == 'R' ){ crystal_table[i] = pos+2;
         } else if( crystal == 'W' ){ crystal_table[i] = pos+3;
         } else if( crystal == 'X' ){ crystal_table[i] = -1; // crystal undefined
         } else {
            fprintf(stderr,"unknown crystal[=%c] in %s\n", crystal, chan_name[i]);
         } break;
      default: break;
      }

      // set full subsystem id (including polarity/output-type etc)
      switch(subsys){
      case ODBHANDLE_GRS: subsys_table[i] = SUBSYS_BGO;       break;
      case ODBHANDLE_SEP: subsys_table[i] = SUBSYS_SCEPTAR;   break;
      case ODBHANDLE_PAC: subsys_table[i] = SUBSYS_PACES;     break;
      case ODBHANDLE_LBS: subsys_table[i] = SUBSYS_LABR_BGO;  break;
      case ODBHANDLE_LBL: subsys_table[i] = SUBSYS_LABR_L;    break;
      case ODBHANDLE_DAL: subsys_table[i] = SUBSYS_LABR_L;    break;
      case ODBHANDLE_DSC: subsys_table[i] = SUBSYS_DESCANT;   break;
      case ODBHANDLE_RCS: subsys_table[i] = SUBSYS_RCMP;      break;
      case ODBHANDLE_DSW: subsys_table[i] = SUBSYS_DESWALL;  break;
      case ODBHANDLE_DSG: subsys_table[i] = SUBSYS_DSG;  break;
      case ODBHANDLE_GRG: subsys_table[i] = (output_type == 1) ? SUBSYS_HPGE_A :SUBSYS_HPGE_B; break;
      case ODBHANDLE_ZDS: subsys_table[i] = (output_type == 1) ? SUBSYS_ZDS_A  :SUBSYS_ZDS_B;  break;
      case ODBHANDLE_ART: subsys_table[i] = (polarity_table[i] == 1) ? SUBSYS_ARIES_A:SUBSYS_ARIES_B;break;
      case ODBHANDLE_XXX: subsys_table[i] = SUBSYS_IGNORE;    break;
      case ODBHANDLE_UNK: subsys_table[i] = SUBSYS_UNKNOWN;   break;
      case ODBHANDLE_DAT: if(crystal_table[i]<8){ subsys_table[i] = SUBSYS_TAC_LABR;
                          }else{ subsys_table[i] = SUBSYS_TAC_ZDS; }
                          break;
      case ODBHANDLE_LBT: if(crystal_table[i]<8){ subsys_table[i] = SUBSYS_TAC_LABR;
                          }else if(crystal_table[i]>8){ subsys_table[i] = SUBSYS_TAC_ART;
                          }else{ subsys_table[i] = SUBSYS_TAC_ZDS; }
                          break;
      }
   }
   memset(subsys_initialized, 0, sizeof(int)*MAX_SUBSYS );

   return(0);
}