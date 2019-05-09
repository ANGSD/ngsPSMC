/*
  GNU license or whetever its called

  vlshchur@gmail.com
  thorfinn@binf.ku.dk
*/

#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <signal.h>
#include <cassert>
#include <pthread.h>
#include <unistd.h>
#include <htslib/hts.h>
#include <htslib/kstring.h> //<- included directly in this source due to some string errors in latest htslib 17sep2019
#include <libgen.h>
#include "psmcreader.h"
#include "header.h"
#include "main_psmc.h"
#include "version.h"

int SIG_COND =1;
double ttol = 1e-16; 
int nThreads =1;
int really_kill =3;
int VERBOSE = 1;
void handler(int s) {
  if(s==13)//this is sigpipe
    exit(0);
  if(VERBOSE)
    fprintf(stderr,"\n\t-> Caught SIGNAL: Will try to exit nicely (no more threads are created.\n\t\t\t  We will wait for the current threads to finish)\n");
  
  if(--really_kill!=3)
  fprintf(stderr,"\n\t-> If you really want \'ngsPSMC\' to exit uncleanly ctrl+c: %d more times\n",really_kill+1);
  fflush(stderr);
  if(!really_kill)
    exit(0);
  VERBOSE=0;
  SIG_COND=0;

}

int print_header(int argc,char **argv){

  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx \n");
    return 0; 
  }
  
  args *pars = getArgs(argc,argv,1);
  if(!pars)
    return 1;
  
  writepsmc_header(stdout,pars->perc,0);
  destroy_args(pars);
  return 0;
}

void toGl(double ingl,double outgl[2]){
  if(std::isinf(ingl)){
    outgl[0]=outgl[1]=log(0);
  }else{
    if(ingl<0){
      outgl[0]=0;
      outgl[1]=ingl;
    }else{
      outgl[0]=-ingl;
      outgl[1]=0;
    }
  }
}


int print_main(int argc,char **argv){
  
  if(argc<1){
    fprintf(stderr,"\t-> Must supply afile.saf.idx files \n");
    fprintf(stderr,"\t-> Examples \n");
    fprintf(stderr,"\t-> ./ngsPSMC print pop1.saf.idx \n");
    fprintf(stderr,"\t-> ./ngsPSMC print pop1.saf.idx -r chr1:10000000-12000000\n");
    return 0; 
  }
  
  args *pars = getArgs(argc,argv,1);
  if(!pars)
    return 0;
  //  writepsmc_header(stderr,pars->perc,0);
  
  for(myMap::iterator it=pars->perc->mm.begin();it!=pars->perc->mm.end();++it){
    

    rawdata rd = readstuff(pars->perc,pars->chooseChr!=NULL?pars->chooseChr:it->first,pars->blocksize,-1,-1);
    fprintf(stderr,"\t-> gls %f %f\n",rd.gls[2],rd.gls[3]);
    double tmp[2];
    
    for(size_t s=rd.firstp;s<rd.lastp;s++){
      toGl(rd.gls[s],tmp);
      fprintf(stdout,"%s\t%d\t%e\t%e\n",it->first,rd.pos[s]+1,tmp[0],tmp[1]);
    }
    delete [] rd.pos; delete [] rd.gls;
    if(pars->chooseChr!=NULL)
      break;
  }
  
  destroy_args(pars);
  return 0;
}



double em(double &x,mygltype *gler,size_t nSites,double tol,int nIter){
  fprintf(stderr,"[%s] x:%f gls:%p nSites:%lu\n",__FUNCTION__,x,gler,nSites);
  fflush(stderr);

  double est = 0;
  double start = x;
  double lastllh =0;

  for(int iter=0;iter<nIter;iter++){
    size_t efsize=0;
    fprintf(stderr,"\r %d/%d     ",iter,nIter);fflush(stderr);
    double llh = 0;
    for(size_t i=0;i<nSites;i++) {
      double tmp[2];
      if(std::isinf(gler[i]))
	continue;
      double igl[2];
      toGl(gler[i],igl);
      //fprintf(stderr,"igl[%lu]=(%f,%f)\n",i,igl[0],igl[1]);
      tmp[0] = exp(igl[0])/4.0*(1-start);
      tmp[1] = exp(igl[1])/6.0*(start);
      est += tmp[1]/(tmp[0]+tmp[1]);
      llh -= log(tmp[0]+tmp[1]);
      efsize++;
    }
    est = est/((double) efsize);
    fprintf(stderr,"iter:%d est: %f llh: %f diffInLlh:%e diffInPars:%e efsize:%lu nsites:%lu\n",iter,est,llh,llh-lastllh,est-start,efsize,nSites);
    if(iter>0&&llh>lastllh){
      fprintf(stderr,"\t-> Problem with EM newllh is larger than lastllh, will break\n");
      break;
    }
    if(iter>0&&fabs(est-start)<tol){
      fprintf(stderr,"\t-> Difference in estimated pars is smaller than tol, convergence achieved\n");
      start=est;
      lastllh=llh;
      break;
    }
    start=est;
    lastllh=llh;
  }
  x=start;
  return lastllh;
}

void calcpost(double &x,mygltype *gler,int nSites,double *pp){
  double start =x;
  for(int i=0;i<nSites;i++) {
    //  fprintf(stderr,"gls=(%f,%f)\n",gls[2*i],gls[2*i+1]);
    double tmp[2];
    double gls[2];
    toGl(gler[i],gls);
    tmp[0] = exp(gls[0])/4.0*(1-start);
    tmp[1] = exp(gls[1])/6.0*(start);
    double est = tmp[1]/(tmp[0]+tmp[1]);
    pp[i] = est;
    //    fprintf(stdout,"%f\n",est);
  }

}

void writefa(kstring_t *kstr,double *positInt,int regLen,int block, int NBASE_PER_LINE,double cutoff){
  int at=0;
    for(int i=0;i<regLen;i+=block) {
      int isHet =0;
      for(int b=0;b<block;b++)
	if(positInt[i+b]>cutoff)
	  isHet++;
      if(at==NBASE_PER_LINE){
	ksprintf(kstr,"\n");
	at=0;
      }
      if(isHet)
	ksprintf(kstr,"K");
      else
	ksprintf(kstr,"T");
      at++;
    }
}


int makevcf2fq(int argc,char **argv){
  if(argc<1){
    fprintf(stderr,"\t-> output is a vcf2fq style file \n");
    return 0; 
  }
  args *pars = getArgs(argc,argv,1);
  if(!pars)
    return 0;
  writepsmc_header(stderr,pars->perc,1);
  
  fprintf(stderr,"nSize: %lu\n",pars->perc->nSites);
  mygltype *gls = new mygltype[pars->perc->nSites];
  size_t at=0;
  //first pull all the data
  for(myMap::iterator it=pars->perc->mm.begin();it!=pars->perc->mm.end();++it){//loop over chrs
    rawdata rd = readstuff(pars->perc,pars->chooseChr!=NULL?pars->chooseChr:it->first,pars->blocksize,-1,-1);
    //    fprintf(stderr,"it->first:%s\tlast:%lu\n",it->first,pars->perc->last);
    memcpy(gls+at,rd.gls,sizeof(mygltype)*rd.lastp);
    at += rd.lastp;
    delete [] rd.pos; delete [] rd.gls;
    if(pars->chooseChr!=NULL)
      break;
  }
  fprintf(stderr,"\t-> at: %lu\n",at);
  double opt = 0.01;
  double llh = em(opt,gls,at,1e-12,500);
  fprintf(stderr,"estimated het:%f with llh:%f\n",opt,llh);
  delete [] gls;

  kstring_t kstr;kstr.l=kstr.m=0;kstr.s=NULL;
  for(myMap::iterator it=pars->perc->mm.begin();it!=pars->perc->mm.end();++it){
    rawdata rd = readstuff(pars->perc,pars->chooseChr!=NULL?pars->chooseChr:it->first,pars->blocksize,-1,-1);

    ksprintf(&kstr,">%s\n",it->first);
    double *pp = new double[rd.lastp];
    calcpost(opt,rd.gls,rd.lastp,pp);
    writefa(&kstr,pp,rd.lastp,100,50,0.9);
    if(kstr.l>0&&kstr.s[kstr.l-1]!='\n')
      ksprintf(&kstr,"\n");
    fwrite(kstr.s,sizeof(char),kstr.l,stdout);
    kstr.l=0;
    delete [] pp;
    delete [] rd.pos; delete [] rd.gls;
    if(pars->chooseChr!=NULL)
      break;
  }
  free(kstr.s);

  destroy_args(pars);
  return 0;
}

int main(int argc,char **argv){
  fprintf(stderr,"\t-> ngsPSMC version: %s (htslib: %s) build(%s %s)\n",ngsPSMC_VERSION,hts_version(),__DATE__,__TIME__);
  fprintf(stderr,"\t-> sizeof(double): %lu sizeof(float): %lu sizeof(mygltype): %lu\n",sizeof(double),sizeof(float),sizeof(mygltype));
  int argc_orig=argc;
  char **argv_orig=argv;
  //start of signal handling
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags =  0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  
  clock_t t=clock();
  time_t t2=time(NULL);

  if(argc==1){
    fprintf(stderr, "\t-> ---./ngsPSMC\n");
    fprintf(stderr,"\t-> ./ngsPSMC [print print_header vcf2fq] afile.psmc.idx \n");
    fprintf(stderr,"\t-> ./ngsPSMC -tole -maxIter -winSize -RD -nThreads -nIter -p -tkfile -nSites -seed -infile -doLinear -nChr\n");
    return 0;
  }

  ++argv;
  --argc;
  if(!strcasecmp(*argv,"print"))
    print_main(--argc,++argv);
  else if(!strcasecmp(*argv,"print_header"))
    print_header(--argc,++argv);
  else if(!strcasecmp(*argv,"vcf2fq"))
    makevcf2fq(--argc,++argv);
  else {
    fprintf(stdout,"MM\t");
    for(int i=0;i<argc_orig;i++)
      fprintf(stdout," %s",argv_orig[i]);
    fprintf(stdout,"\n");
    fprintf(stdout,"MM\t");
    fprintf(stdout,"\t-> ngsPSMC version: %s (htslib: %s) build(%s %s)\n",ngsPSMC_VERSION,hts_version(),__DATE__,__TIME__); 
    fprintf(stdout,"MM\t TR is theta rho, theta in this context is given by the persite theta, and not the per window theta, this is different from the original PSMC\n");
    if(isatty(fileno(stdout))){
      fprintf(stderr,"\t-> You are printing results to the terminal consider dumping into a file\n");
      fprintf(stderr,"\t-> E.g.: \'./ngsPSMC ");
      for(int i=0;i<argc;i++)
	fprintf(stderr," %s",argv[i]);
      fprintf(stderr," >psmc.ml.txt\'\n");   
    }

    main_psmc(argc,argv);
  }

  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  return 0;
}
