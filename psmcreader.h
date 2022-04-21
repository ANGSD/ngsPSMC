#pragma once
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <map>
#include <htslib/bgzf.h>
#include "fastas.h"
#include "header.h"

//#define GL_AS_CHAR

#ifdef GL_AS_CHAR
typedef char mygltype;
#else
typedef double mygltype;
#endif 

typedef struct{
  mygltype *gls;//2*len homo or hetero read
  int *pos;//positions of start for considering chromosome
  size_t len;
  size_t firstp;//if we have specified a region, then this is the first index to use
  size_t lastp;//if we have specified a region, then this is the last index to use
}rawdata;

// information about input file(FOR saf.idx)
typedef struct{
  size_t nSites;
  myMap mm;// Dictionary for (chromosome : [datm = [nSites, Pos, sAF (Sample allel frequency)]]
  char *bgzf_pos;//name of input.psmc.gz file
  char *bgzf_gls;//name of input.psmc.pos.gz file
  int version;//1 is PSMC, 2 is VCF, otherwise assuming fasta
  char *fname;//input.saf.idx? NOW VCF MAY BE IMPLEMENTET
  perFasta *pf;
}infstruct;

infstruct* infstruct_init(char *fname,int nChr);
int psmcversion(const char *fname);
int vcfversion(const char *fname);
void writepsmc_header(FILE *fp,infstruct *pp,int onlysubset);
void infstruct_destroy(infstruct *pp);
rawdata readstuff(infstruct *pp,char *chr,int blockSize,int start,int stop);
