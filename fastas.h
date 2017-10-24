#pragma once
#include <htslib/faidx.h>

typedef struct{
  char *fastaname;
  faidx_t *fai;//contains the faidx structure
  char *seqs;//contains the reference for the current chr;
  int curChr;//the exact chromosome name for the seqs above
  int chrLen;//length of chromosome
}perFasta;


//this will initialize our data
perFasta *perFasta_init(const char *fname);
char *loadChr(perFasta *f, char*chrName,int chrId);
