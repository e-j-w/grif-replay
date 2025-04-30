#ifndef SMOLFMT_H
#define SMOLFMT_H

//structs defining the SMOL tree format used by Jonathan for sorting data,
//which can then be processed by external sort codes

//this format is optimized for the smallest possible filesize, as disk read is the
//typical bottleneck when sorting

//SMOL tree files consist of evt_header structs written to disk, followed
//by N hpge_hit structs, where N=evt_header->numHPGeHits.
//Each even is written sequentially to disk.
//At the start of the SMOL tree file, there is a single uint64_t header value,
//formatted as follows:
// - bits 0 - 48:  totalnumber of events in the file
// - bits 49 - 63: SMOL tree type (default 0, can be used to specify different tree formats)
//Ideally, the SMOL tree type information can be used to specify different sorted_evt formats,
//so that different experiment types can have optimized tree formats only containing the info
//relevant for that experiment (so that eg. ancillary detector metadata isn't saved for 
//experiments that only use GRIFFIN)

#include <stdint.h> //allows uint8_t and similiar types

#define MAX_EVT_HIT 64

typedef struct
{
	uint8_t numHPGeHits; //number of HPGe array hits
	uint8_t metadata; //bit 0: TIGRESS (on) or GRIFFIN (off), bit 1: any suppressor fired, bit 7: always set, for data validation
	double evtTimeNs; //time of the first hit, in ns - float doesn't have enough precision
}evt_header;

typedef struct
{
	float timeOffsetNs; //relative to evtTimeNs
	float energy; //calibrated energy in keV
	uint8_t core; //0-indexed
}hpge_hit;

//struct to hold an event resident in memory,
//not used in the actual file written to disk
typedef struct
{
	evt_header header;
	hpge_hit hpgeHit[MAX_EVT_HIT]; //non-addback
}sorted_evt;

#endif
