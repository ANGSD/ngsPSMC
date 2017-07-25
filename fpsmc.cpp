#include <vector>
#include <cassert>
#include <cmath>
#include <pthread.h>
#include "psmcreader.h"
#include "main_psmc.h"
#include "hmm_psmc.h"
#include "bfgs.h"
#include <errno.h>

int nThreads =1;

int nChr =0;


typedef struct{
  double *tk;
  int tk_l;
  double *epsize;
  double rho;
}shared_forhmm;

shared_forhmm shmm;


fastPSMC **objs = NULL;

typedef struct{
  double **nP;
  double **PP;
  double *tk;
  int tk_l;
  double pix;
  int numWind;

  double rho;
  double *epsize;
}oPars;

/*
  objective function. Function to be optimized
*/

double qFunction(const double *params ,const void *d){
  oPars *data = (oPars*) d;
  void ComputeGlobalProbabilities(double *tk,int tk_l,double **P,double *epsize,double rho);
  ComputeGlobalProbabilities(data->tk,data->tk_l,data->nP,data->epsize,data->rho);
  double Q = 0;
  for (unsigned i = 0; i < data->tk_l; i++)
    Q += qkFunction(i, data->pix,data->numWind,data->nP,data->PP);
  return Q;
}


void runoptim(double *tk,int tk_l,double *epsize,double rho,double **PP,double pix,int numWin){
  double **nP = new double *[8];
  for(int i=0;i<8;i++)
    nP[i] = new double[tk_l];

  //get start
  double pars[tk_l];
  for(int i=0;i<tk_l;i++)
    pars[i] = drand48();
  //set bounds
  int nbd[tk_l];
  double lbd[tk_l];
  double ubd[tk_l];
  for(int i=0;i<tk_l;i++){
    nbd[i]=1;
    lbd[i]=0.000001;
    ubd[i]=PSMC_T_INF;
  }

  oPars data;
  data.nP = nP;
  data.PP = PP;
  data.tk = tk;
  data.tk_l = tk_l;
  data.pix = pix;
  data.numWind=numWin;
  data.rho= rho;
  data.epsize=epsize;
  
  double max_llh = findmax_bfgs(tk_l,pars,(void *)&data,qFunction,NULL,lbd,ubd,nbd,1);
}

void printarray(FILE *fp,double *ary,int l);
void printmatrix(FILE *fp,double **mat,int x,int y);



void printmatrixf(char *fname,double **m,int x,int y){
  FILE *fp = NULL;
  if(!(fp=fopen(fname,"wb"))){
    fprintf(stderr,"\t-> Problem writing file: \'%s\'\n",fname);
    exit(0);
  }
  printmatrix(fp,m,x,y);
  fclose(fp);
}

void printarrayf(char *fname,double *m,int x){

  FILE *fp = NULL;
  if(!(fp=fopen(fname,"wb"))){
    fprintf(stderr,"\t-> Problem writing file: \'%s\'\n",fname);
    exit(0);
  }
  printarray(fp,m,x);
  fclose(fp);
}

void *run_a_hmm(void *ptr){
  size_t at =(size_t) ptr;
  fprintf(stderr,"at:%lu\n",at);
  sleep(drand48()*10);
  objs[at]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.rho);
  pthread_exit(NULL);
}


void main_analysis(double *tk,int tk_l,double *epsize,double rho){
  shmm.tk=tk;
  shmm.tk_l=tk_l;
  shmm.rho=rho;
  shmm.epsize=epsize;

  pthread_t thread[nThreads];
  if(nThreads==1)
    for(int i=0;i<nChr;i++)
      objs[i]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.rho);
  else {
    int at=0;
    while(at<nChr){
      int thisround = std::min(nChr-at,nThreads);
      for(int t=0;t<thisround;t++){
	size_t index = at+t;
	if(pthread_create( &thread[t], NULL, run_a_hmm, (void*) index)){
	  fprintf(stderr,"[%s] Problem spawning thread\n%s\n",__FUNCTION__,strerror(errno));
	  exit(0);
	}
      }
      for(int t=0;t<thisround;t++){
	if(pthread_join( thread[t], NULL)){
	  fprintf(stderr,"[%s] Problem joining thread\n%s\n",__FUNCTION__,strerror(errno));
	  exit(0);
	}
      }
      at+=thisround;
    }
  }

}


int psmc_wrapper(args *pars,int block) {
#if 1 //print pars
  psmc_par *p=pars->par;
  fprintf(stderr,"par->n:%d\tpar->n_free:%d\tpar_map:%p\tpar->pattern:%s\tpar->times:%p\tpar->params:%p\n",p->n,p->n_free,p->par_map,p->pattern,p->times,p->params);
  for(int i=0;i<pars->par->n+1;i++)
    fprintf(stderr,"%i)\t%f\t%f\n",i,pars->par->times[i],pars->par->params[i]);
  
#endif
  int tk_l = pars->par->n+1;
  double *tk = new double [tk_l];
  double *epsize = new double [tk_l];
  setEPSize(epsize,tk_l,p->params);
  //(nelems,array,max_t,alpha,array with values from file, can be NULL)
  setTk(tk_l,tk,15,0.01,p->times);//<- last position will be infinity
#if 1
  for(int i=0;i<tk_l;i++)
    fprintf(stderr,"[%d]:(%f,%f)\n",i,tk[i],epsize[i]);
#endif
  
  //initialize all hmm (one for each chr), for now just a single
  int nobs = pars->chooseChr?1:pars->perc->mm.size();
  fprintf(stderr,"\t-> nobs/nchr: %d\n",nobs);
  objs = new fastPSMC*[nobs];
  for (myMap::const_iterator it = pars->perc->mm.begin() ;it!=pars->perc->mm.end();it++) {
    if(pars->chooseChr!=NULL)
      iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop);
    else
      iter_init(pars->perc,it->first,pars->start,pars->stop);
    fastPSMC *obj=objs[nChr++]=new fastPSMC;
    obj->setWindows(pars->perc->gls,pars->perc->pos,pars->perc->last,pars->block);
    //    obj->printWindows(stdout);exit(0);
    obj->allocate(tk_l);
    if(pars->chooseChr!=NULL)
      break;
  }
  fprintf(stderr,"\t-> We have now allocated hmm's for: %d chromosomes\n",nChr);
  objs[0]->make_hmm(tk,tk_l,epsize,0.1);
  /*
    printarrayf("stationary",obj.stationary,tk_l);
    printarrayf("tk",tk,tk_l);
    printarrayf("epsize",epsize,tk_l);
    printmatrixf("P",obj.P,7,tk_l);
    printmatrixf("emis",obj.emis,tk_l,obj.windows.size()+1);
    printmatrixf("fw",obj.fw,tk_l,obj.windows.size()+1);
    printarrayf("r1",obj.R1,tk_l);
    printarrayf("r2",obj.R2,tk_l);
  */
  
  //calculate window end points
  
  //  obj.printWindows(stdout);
  //allocate internal structures needed
 
  //make an hmm

  // objs[0].make_hmm(tk,tk_l,epsize);
  /*
    printarrayf("stationary",obj.stationary,tk_l);
    printarrayf("tk",tk,tk_l);
    printarrayf("epsize",epsize,tk_l);
    printmatrixf("P",obj.P,7,tk_l);
    printmatrixf("emis",obj.emis,tk_l,obj.windows.size()+1);
    printmatrixf("fw",obj.fw,tk_l,obj.windows.size()+1);
    printarrayf("r1",obj.R1,tk_l);
    printarrayf("r2",obj.R2,tk_l);
  */
  return 1;
}
