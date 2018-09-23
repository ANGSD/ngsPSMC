#pragma once
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <map>
#include <htslib/bgzf.h>
#include "fastas.h"
#include "header.h"

typedef double mygltype;
typedef struct{
  size_t nSites;
  myMap mm;
  BGZF *bgzf_pos;
  BGZF *bgzf_gls;
  int version;//1 is gl, otherwise assuming fasta
  int *pos;//contains the positions
  size_t pos_l;
  char *fname;//input.saf.idx?
  size_t first;//if we have specified a region, then this is the first index to use
  size_t last;//if we have specified a region, then this is the last index to use
  //  double *gls;//2*nSites long homogl_1 hetgl_1 homogl_2 hetgl_2
  mygltype *gls;//one gl. Postive means its the gl for the hom, negative means its th gl for the het. The other gl is zero. This value is therefore very beautiful. This pointer should be copied to HMM. And the HMM should free it after use.
  perFasta *pf;
  double *tmpgls;//temparray to holddata
  size_t tmpgls_l;//length of temparray to hold data
}perpsmc;

perpsmc* perpsmc_init(char *fname,int nChr);
void writepsmc_header(FILE *fp,perpsmc *pp,int onlysubset);
void perpsmc_destroy(perpsmc *pp);
myMap::iterator iter_init(perpsmc *,char *,int,int,int);
