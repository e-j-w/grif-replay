// TODO ..
//   commented "inuse check for remove global
// remove_gate_from_group
//  gate/grouplist ordering may need to be reset on clearing
// config mtime to ALL new functions
//   commented free(sort->global[i]->name) due to crash

#ifndef CONFIG_H
#define CONFIG_H

#include "midas2smol.h"

int add_sortfile(char *path);
int open_next_sortfiles(Sort_status *arg);
int free_sortfile(Sort_status *sort);
int close_sortfiles(Sort_status *arg);
int end_current_sortfile(int fd);
int pre_sort(int frag_idx, int end_idx);
uint8_t fill_smol_entry(FILE *out, const int win_idx, const int frag_idx);
uint64_t sort_built_event(int window_start, int win_end, FILE *out);

#define DEFAULT_CONFIG "last.json"

typedef struct global_struct {
   char name[STRING_LEN]; int min; int max;
} Global;

typedef struct cal_coeff_struct {
   char name[CHAN_NAMELEN]; float offset; float gain; float quad;
   float pileupk1[7], pileupk2[7], pileupE1[7];
   short address; short datatype;
} Cal_coeff;

typedef struct sortvar_struct {     // sortvars used by histos AND conditions
   int value;      int offset;    int dtype;     //  but "use_count" only used
   int use_count_x;  int valid;   int local;     //  to refer to histo_use
   char name[STRING_LEN]; char title[STRING_LEN];
} Sortvar;

typedef struct cond_struct {      // use count inc/dec when un/applying gates
   char name[STRING_LEN]; Sortvar *var; int op; int value;      // (to histos)
   int use_count;  int veto;
   int passed; int pass_count;  // passed only used during sort
} Cond;

typedef struct gate_struct {
   char name[STRING_LEN]; int nconds; Cond *conds[MAX_GATE_CONDS];
   int use_count;  int passed;
} Gate;

typedef struct histo_folder_struct {
   char name[HISTO_FOLDER_LENGTH];
   struct histo_folder_struct *next_subfolder;  // 1 level lower
   struct histo_folder_struct *next_folder;     // same level
   struct th1i_struct *first_histo;             // same level
} Folder;

#define DISK_CONFIG 0
#define MEM_CONFIG  1

// group config element stuff together, with add/delete/copy fns
//   to allow quick switch between storing data or pointers (malloc/free)
typedef struct config_set_struct { int  type; // memory(live,sort) or disk
   char name[SYS_PATH_LENGTH];     // "live", "sort", or pathname of tar file
   char data_dir[SYS_PATH_LENGTH];   // most recent datafile directory
   char config_dir[SYS_PATH_LENGTH]; // most recent configfile directory
   char midas_title[SYS_PATH_LENGTH];// title of midas run(histo_config)
   char out_file[SYS_PATH_LENGTH]; // output file path
   int  midas_start_time;  int midas_runtime;  int mtime;  int lock;
   int folders_valid;    Folder first_folder;
   int current_depth;    Folder *treepath[HISTO_FOLDER_LENGTH]; // length 2296
   int ncal;             Cal_coeff *calib[MAX_CALIB];
   int nglobal;          Global *globals[MAX_GLOBALS];
   int nconds;           Cond  *condlist[MAX_CONDS];       //   sorted list
   int ngates;           Gate  *gatelist[MAX_GATES];//(inactive at end)           47384
   int nusedvar;         Sortvar *usedvars[MAX_SORT_VARS];
   int nuser;            
   int nsortvar;         Sortvar varlist[MAX_SORT_VARS];                     // 33921336
   Cond cond_array[MAX_GATES];  Gate gate_array[MAX_GATES]; // unsorted
   Global global_array[MAX_GLOBALS];                                         // 34115896
   Cal_coeff calib_array[MAX_CALIB];  int odb_daqsize;
} Config;

extern Config *configs[MAX_CONFIGS]; // start unallocated

extern int init_default_sort(Config *cfg, Sort_status *arg);

extern Config *add_config(char *name);
extern int remove_config(Config *cfg);
extern int next_condname(Config *cfg);
extern int set_calibration(Config *cfg, int num, char url_args[][STRING_LEN], int fd);
extern int set_pileup_correction(Config *cfg, int num, char url_args[][STRING_LEN], int fd);

/////////////////////////////////////////////////////////////////////////
/////////////////////          Gains         ////////////////////////////
/////////////////////////////////////////////////////////////////////////
extern int edit_calibration(Config *cfg, char *name, float offset, float gain, float quad, float pileupk1[7], float pileupk2[7], float pileupE1[7], int address, int type, int overwrite);

/////////////////////////////////////////////////////////////////////////
/////////////////////       Variables        ////////////////////////////
/////////////////////////////////////////////////////////////////////////

/////// sort variables are currently predefined and fixed at compilation time
/////// maybe later will add ability to add/remove them during runtime
//extern int add_variable(char *name);
//extern int remove_variable(char *name);
extern Sortvar *find_sortvar(Config *cfg, char *name);

extern int add_global(Config *cfg, char *name, int value, int val2);
extern int remove_global(Config *cfg, char *name);

/////////////////////////////////////////////////////////////////////////
/////////////////          Config Files         /////////////////////////
/////////////////////////////////////////////////////////////////////////

extern int init_config();
extern int init_user_config(Config *cfg);
extern int init_default_config(Config *cfg);
extern int write_config(Config *cfg, FILE *fp);
extern int copy_config(Config *src, Config *dst);
extern int clear_config(Config *cfg);
extern int clear_calibrations(Config *cfg);
extern int delete_config(Config *cfg);
extern int load_config(Config *cfg, char *filename, char *buffer);
extern int save_config(Config *cfg, char *filename, int overwrite);

extern int set_directory(Config *cfg, char *name, char *path);
extern int set_midas_param(Config *cfg, char *name, char *value);

#define HPGeE        0 // HPGE
#define HPGeA        1 //
#define HPGeT        2 //
#define HPGeTS       3 //
#define HPGePH       4 //
#define HPGeC        5 //
#define HPGeCL       6 //
#define HPGePU       7 //
#define HPGeIP       8 //
#define HPGeDT       9 //
#define HPGeEU      10 // #
#define HPGeAU      11 //
#define GRGTHETA    12 //
#define GRGPHI      13 //
#define CLVTHETA    14 //
#define CLVPHI      15 //
#define SEPE        16 // SCEPTAR
#define SEPT        17 //
#define SEPTS       18 //
#define EPPH        19 //
#define SEPPU       20 //
#define SEPTHETA    21 //
#define SEPPH       22 //
#define SEPNUM      23 //
#define PACE        24 // PACES
#define PACT        25 //
#define PACTS       26 //
#define PACPH       27 //
#define PACPU       28 //
#define PACTHETA    29 //
#define PACPHI      30 //
#define PACNUM      31 //
#define LBLE        32 // LaBr3
#define LBLT        33 //
#define LBLTS       34 //
#define LBLPH       35 //
#define LBLPU       36 //
#define LBLTHETA    37 //
#define LBLPHI      38 //
#define LBLNUM      39 //
#define LBT         40 // TACs
#define LBTT        41 //
#define LBTTS       42 //
#define LBTPH       43 //
#define LBTPU       44 //
#define LBTNUM      45 //
#define GRSE        46 // Clover BGO
#define GRST        47 //
#define GRSTS       48 //
#define GRSPH       49 //
#define GRSPU       50 //
#define GRSNUM      51 //
#define GRSPOS      52 //
#define GRSTYPE     53 //
#define LBSE        54 // Ancillary BGO
#define LBST        55 //
#define LBSTS       56 //
#define LBSPH       57 //
#define LBSPU       58 //
#define LBSNUM      59 //
#define LBSPOS      60 //
#define MIDAS_Time  61 // Time Differences
#define TD_GRG_GRG  62 //
#define TD_SEP_SEP  63 //
#define TD_PAC_PAC  64 //
#define TD_LBL_LBL  65 //
#define TSD_GRG_GRG 66 //
#define TSD_SEP_SEP 67 //
#define TSD_PAC_PAC 68 //
#define TSD_LBL_LBL 69 //
#define TD_GRG_SEP  70 //
#define TSD_GRG_SEP 71 //
#define TD_GRG_PAC  72 //
#define TSD_GRG_PAC 73 //
#define TD_SEP_PAC  74 //
#define TSD_SEP_PAC 75 //
#define TD_GRG_LBL  76 //
#define TSD_GRG_LBL 77 //
#define TD_SEP_LBL  78 //
#define TSD_SEP_LBL 79 //
#define ANG_GRG_GRG 80 // Angular Differences [HpGeDistDependant]
#define ANG_CLV_CLV 81 //
#define ANG_SEP_SEP 82 //
#define ANG_PAC_PAC 83 //
#define ANG_LBL_LBL 84 //
#define ANG_GRG_SEP 85 //
#define ANG_GRG_PAC 86 //
#define ANG_GRG_LBL 87 //
#define ANG_PAC_LBL 88 //
#define ANG_SEP_LBL 89 //
#define PPG_NUM     90 // Cycle Timing (PPG events)
#define PPG_TIME    91 //
#define PPG_PAT     92 //

#endif
