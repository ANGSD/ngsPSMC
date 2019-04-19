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
  mygltype *gls;//2*len
  int *pos;//len
  size_t len;
  size_t firstp;//if we have specified a region, then this is the first index to use
  size_t lastp;//if we have specified a region, then this is the last index to use
}rawdata;

typedef struct{
  size_t nSites;
  myMap mm;
  char *bgzf_pos;
  char *bgzf_gls;
  int version;//1 is gl, otherwise assuming fasta
  char *fname;//input.saf.idx?
  perFasta *pf;
}perpsmc;

perpsmc* perpsmc_init(char *fname,int nChr);
void writepsmc_header(FILE *fp,perpsmc *pp,int onlysubset);
void perpsmc_destroy(perpsmc *pp);
rawdata readstuff(perpsmc *pp,char *chr,int blockSize,int start,int stop);
