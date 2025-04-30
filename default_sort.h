// default_sort.h
// Header file for default_sort.c


//#######################################################################
//########         Subsystem and Detector definitions          ##########
//#######################################################################

// do not alter order without also changing subsys_e_vs_e, subsys_dt
#define MAX_SUBSYS      24
#define SUBSYS_HPGE_A    0
#define SUBSYS_PACES     1
#define SUBSYS_LABR_L    2
#define SUBSYS_RCMP      3
#define SUBSYS_ARIES_A   4
#define SUBSYS_ZDS_A     5 // GRIF16
#define SUBSYS_TAC_LABR  6
#define SUBSYS_LABR_BGO  7
#define SUBSYS_BGO       8
#define SUBSYS_SCEPTAR   9
#define SUBSYS_DESCANT  10
#define SUBSYS_DESWALL  11
#define SUBSYS_DSG      12
#define SUBSYS_IGNORE   13
#define SUBSYS_HPGE_B   16
#define SUBSYS_ARIES_B  17 // CAEN
#define SUBSYS_ZDS_B    18 // CAEN
#define SUBSYS_TAC_ZDS  19
#define SUBSYS_TAC_ART  20
#define SUBSYS_UNKNOWN  23
static char subsys_handle[MAX_SUBSYS][8] = {
  "GRGA", "PAC",  "LBL",  "RCS",
  "ARTA", "ZDSA", "LBT",  "LBS",
  "BGO",  "SEP",  "DSC",  "DSW",
  "DSG", "XXX1", "XXX2", "XXX3",
  "GRGB", "ARTB", "ZDSB", "", // secondary names start after #16
  "",     "",     "",     "UNK"
};
static char subsys_name[MAX_SUBSYS][STRING_LEN] = {
  "Griffin", "PACES",   "LaBrX",   "RCMP",     //  0- 3
  "ARIES",   "ZDSA",    "TAC_LBL",   "LaBrS",    //  4- 7
  "BGO",     "Sceptar", "Descant", "DES_WALL", //  8-11
  "Des_Ancil", "Ignore1", "Ignore2", "Ignore3",  // 12-15
  "Grif_B",  "ARS_B",   "ZDS_B",   "TAC_ZDS",      // 16-19
  "TAC_ART",      "",        "",        "Unknown"   // 20-23
}; // final entry will be used if not found - make sure it is not empty
// #####################################################################

#define N_CLOVER 16
#define N_HPGE 64
#define N_ARIES 76
#define N_LABR 8
#define N_TACS 12
#define N_RCMP_POS 6
#define N_RCMP_STRIPS 32
#define N_DES_WALL 60

//#######################################################################
//########                PRESORT Time Gates                   ##########
//#######################################################################

// The definition of the time difference gate in 10 nanosecond units.
// The value is the maximum time difference in 10 nanosecond units.
// The default values set here are replaced by the Global value at start of sorting.
static int bgo_window_min = 0;
static int bgo_window_max = 20;
