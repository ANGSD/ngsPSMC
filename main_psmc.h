#pragma once
#include "psmcreader.h"
#include "fastas.h"

//from psmc
typedef struct{
  int n; // $n$ in psmc.tex. number of intervals equals to $n+1$
  int n_free; // number of free lambdas
  int *par_map; // parameter groups
  char *pattern;
  double *times;//<-splittimes
  double *params;//<- effective population sizes
  double TR[2];
}psmc_par;

typedef struct{
  std::vector<double> eN1;
  std::vector<double> eN2;
  double rho[2];
  double theta;
}msarg;



typedef struct {
  int nChr;
  int smartsize;
  char *chooseChr;
  int start;
  int stop;
  size_t nSites;
  int maxIter;
  double tole;
  perpsmc * perc;
  char *fname;
  int onlyOnce;
  long seed;//<-seed=-1 old version;seed=0 means time; othervise it will be used as seed
  int blocksize; //in bp
  psmc_par *par;
  perFasta *pf;
  int RD;
  int nThreads;
  int nIter;
  int doLinear;
  char *outname;
  double init;
  char *msstr;
  msarg msstr_arg;
  char *psmc_infile;
  double init_theta;
  double init_rho;
  double init_max_t;
}args;
args * getArgs(int argc,char **argv,int dontprint);
void destroy_args(args *p);
int main_psmc(int argc,char **argv);
