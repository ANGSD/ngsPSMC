#include <vector>
#include <cassert>
#include <cmath>
#include <pthread.h>
#include "psmcreader.h"
#include "main_psmc.h"
#include "hmm_psmc.h"
#include "bfgs.h"
#include "compute.h"
#include <errno.h>

extern int nThreads;

int nChr = 0;


int doQuadratic = 0; //<-only used in qFunction_wrapper


typedef struct{
  double *tk;
  int tk_l;
  double *epsize;
  double theta;
  double rho;
  double **trans;
}shared_forhmm;


typedef struct{
  double **nP;
  double **PP;
  double **baumwelch;
  double *tk;
  int tk_l;
  double pix;
  int numWind;
  double theta;
  double rho;
  double **trans;
  //  double *epsize;
  double llh;
  double *parsIn;
}oPars;

int remap_l;
int  *remap;
shared_forhmm shmm;


fastPSMC **objs = NULL;
oPars *ops = NULL;

/*
  objective function. Function to be optimized, for each chromo
*/

double qFunction_inner(double *tk,int tk_l,const double *epsize,double rho,double pix,int numWind,double **nP,double **PP);


double qFunction(const double *params ,const void *d){
  oPars *data = (oPars*) d;
  //  fprintf(stderr,"pix:%f\n",data->pix);
  return qFunction_inner(data->tk,data->tk_l,params,data->rho,data->pix,data->numWind,data->nP,data->PP);

}


double qFunction_inner2(double *tk,int tk_l,const double *epsize,double rho,double pix,int numWind,double **nP,double **PP,double **trans);


double qFunction2(const double *params ,const void *d){
  oPars *data = (oPars*) d;
  //  fprintf(stderr,"pix:%f\n",data->pix);
  return qFunction_inner2(data->tk,data->tk_l,params,data->rho,data->pix,data->numWind,data->nP,data->baumwelch,data->trans);

}

extern int SIG_COND;//<- used for killing program
static int ncals=0;
void *qFunction2_thd(void *ptr);
void *qFunction_thd(void *ptr);
double qFunction_wrapper(const double *pars,const void *){
  ncals++;
  //fprintf(stderr,"\t-> calling objective function: remap_l:%d [%d]\n",remap_l,ncals);
  double pars2[ops[0].tk_l];
  int at=0;
  for(int i=0;i<remap_l;i++)
    for(int j=0;j<remap[i];j++){
      pars2[at++] = pars[i]; 
      //      fprintf(stderr,"\t-> pars2: %e\n",pars2[at-1]);
    }
  double ret =0;
  if(nThreads==1){
    for(int i=0;SIG_COND&&i<nChr;i++)
      if(doQuadratic)
	ret += qFunction2(pars2,&ops[i]);
      else
	ret += qFunction(pars2,&ops[i]);
  }else {
    pthread_t thread[nThreads];
    for(int i=0;i<nChr;i++)
      ops[i].parsIn = pars2;
    int at=0;
    while(SIG_COND&&at<nChr){
      int thisround = std::min(nChr-at,nThreads);
      for(int t=0;t<thisround;t++){
	size_t index = at+t;
	if(doQuadratic){
	  if(pthread_create( &thread[t], NULL, qFunction2_thd, (void*) index)){
	    fprintf(stderr,"[%s] Problem spawning thread\n%s\n",__FUNCTION__,strerror(errno));
	    exit(0);
	  }
	}else{
	   if(pthread_create( &thread[t], NULL, qFunction_thd, (void*) index)){
	    fprintf(stderr,"[%s] Problem spawning thread\n%s\n",__FUNCTION__,strerror(errno));
	    exit(0);
	  }
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
    for(int i=0;i<nChr;i++)
      ret += ops[i].llh;

  }    
  if(std::isinf(ret)||SIG_COND==0)
    ret= -1000000000;
  
  //  fprintf(stderr,"qfun:%e\n",ret);
  return -ret;
}


void make_remapper(psmc_par *pp){
  //  fprintf(stderr,"makeadsfadsfasfasfadf\n");
  int at=0;
  remap = new int[pp->n_free];
  remap_l=pp->n_free;
  for(int i=0;i<pp->n_free;i++){
    int howMany=0;
    while(pp->par_map[at+howMany]==pp->par_map[at])
      howMany++;
    at+=howMany;
    remap[i] = howMany;
  }
  
#if 0
  for(int i=0;i<pp->n_free;i++)
    fprintf(stderr,"[%d]: %d\n",i,remap[i]);
#endif
  
}

static int  mysupercounter =0;

void runoptim3(double *tk,int tk_l,double *epsize,double theta,double rho,psmc_par *pp,FILE *FLOG,double &ret_qval){
  clock_t t=clock();
  time_t t2=time(NULL);
  fprintf(stderr,"\t-> Starting Optimization\n");
  fflush(stderr);
#if 0
  fprintf(stderr,"pp->n:%d\n",pp->n);
  fprintf(stderr,"pp->n_free:%d\n",pp->n_free);
  for(int i=0;i<pp->n+1;i++)
    fprintf(stderr,"pars_map[%d]\t%d\n",i,pp->par_map[i]);
#endif

  int ndim = pp->n_free;
  fprintf(stderr,"\t-> ndim:%d\n",ndim);

  double pars[ndim];
 

  int at=0;
  for(int i=0;i<remap_l;i++){
    pars[i] = epsize[at+remap[i]-1];
    at+=remap[i];
  }

  //set bounds
  int nbd[ndim];
  double lbd[ndim];
  double ubd[ndim];
  for(int i=0;i<ndim;i++){
    nbd[i]=2;
    lbd[i]=0.000001;
    ubd[i]=1000;//PSMC_T_INF;
  }

  for(int i=0;i<nChr;i++){
    ops[i].nP = objs[i]->nP;
    ops[i].PP = objs[i]->PP;
    ops[i].baumwelch  = objs[i]->baumwelch;
    ops[i].tk = tk;
    ops[i].tk_l = tk_l;
    ops[i].pix = objs[i]->pix;
    ops[i].numWind=objs[i]->windows.size();
    ops[i].rho= rho;
    ops[i].theta = theta;
    ops[i].trans = objs[i]->trans;
    //    fprintf(stderr,"trans[0][0]\n",ops[i].trans[0][0]);exit(0);
  }
  ncals=0;
  //  mysupercounter=0;
  fprintf(FLOG,"\tpreopt[%d]",mysupercounter);
  for(int i=0;i<ndim;i++)
     fprintf(FLOG,"\t%f",pars[i]);
  fprintf(FLOG,"\n");
  //we are not optimizing llh but qfunction
  double max_qval = findmax_bfgs(ndim,pars,NULL,qFunction_wrapper,NULL,lbd,ubd,nbd,-1);
  ret_qval=max_qval;
  fprintf(FLOG,"\tpostopt[%d]",mysupercounter);
  for(int i=0;i<ndim;i++)     
    fprintf(FLOG,"\t%f",pars[i]);
  fprintf(FLOG,"\n");
  fprintf(FLOG,"\t-> optim done: after ncalls:%d best total qval:%f\n",ncals,max_qval);
  for(int i=0;0&i<ndim;i++)
    fprintf(stderr,"newpars after optim[%d]: %f\n",mysupercounter, pars[i]);
  fflush(stderr);fflush(FLOG);
  mysupercounter++;
  at=0;
  for(int i=0;i<remap_l;i++){
    for(int j=0;j<remap[i];j++){
      //      fprintf(stderr,"epsize[%d]: %f\n",at,epsize[at]);
      epsize[at++] = pars[i];
    }
  }
  for(int i=0;0&&i<at;i++)
    fprintf(stderr,"[%d]: %f\n",i,epsize[i]);

  fflush(FLOG);;
  fprintf(stderr, "\t-> [RUNOPTIM3 TIME]:%s cpu-time used =  %.2f sec \n",__func__, (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t-> [RUNOPTIM3 Time]:%s walltime used =  %.2f sec \n",__func__, (float)(time(NULL) - t2));
 
  fflush(stdout);
  fflush(stderr);
}

void printarray(FILE *fp,double *ary,int l);
void printmatrix(FILE *fp,double **mat,int x,int y);



void printmatrixf(char *fname,double **m,int x,int y){
  //  return ;
  FILE *fp = NULL;
  if(!(fp=fopen(fname,"wb"))){
    fprintf(stderr,"\t-> Problem writing file: \'%s\'\n",fname);
    exit(0);
  }
  printmatrix(fp,m,x,y);
  fclose(fp);
}

void printmatrixf2(char *fname,int index,double **m,int x,int y){
  //  return ;
  char tmpnam[1024];
  snprintf(tmpnam,1024,"%s_%d",fname,index);
  fprintf(stderr,"\t-> tmpnam:%s\n",tmpnam);
  printmatrixf(tmpnam,m,x,y);
}


void printmatrixf3(char *fname,char *fname2,int index,double **m,int x,int y){
  //  return ;
  char tmpnam[1024];
  snprintf(tmpnam,1024,"%s_%s_%d",fname,fname2,index);
  fprintf(stderr,"\t-> tmpnam:%s\n",tmpnam);
  printmatrixf(tmpnam,m,x,y);
}

void printarrayf(char *fname,double *m,int x){
  //  return;
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
  //  fprintf(stderr,"at:%lu\n",at);
  //  sleep(drand48()*10);
  objs[at]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.theta,shmm.rho);
  pthread_exit(NULL);
}

void *qFunction2_thd(void *ptr){
  size_t at =(size_t) ptr;
  //  fprintf(stderr,"at:%lu\n",at);
  //  sleep(drand48()*10);
  ops[at].llh = qFunction2(ops[at].parsIn,&ops[at]);
  //  objs[at]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.theta,shmm.rho);
  pthread_exit(NULL);
}


void *qFunction_thd(void *ptr){
  size_t at =(size_t) ptr;
  //  fprintf(stderr,"at:%lu\n",at);
  //  sleep(drand48()*10);
  ops[at].llh = qFunction(ops[at].parsIn,&ops[at]);
  //  objs[at]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.theta,shmm.rho);
  pthread_exit(NULL);
}



void main_analysis_make_hmm(double *tk,int tk_l,double *epsize,double theta,double rho,double &ret_llh,double &ret_qval){

  fprintf(stderr,"\t-> [%s:%s:%d] nthreads:%d tk_l:%d theta:%f rho:%f\n",__FILE__,__FUNCTION__,__LINE__,nThreads,tk_l,theta,rho);
  shmm.tk=tk;
  shmm.tk_l=tk_l;
  shmm.theta=theta;
  shmm.rho=rho;
  shmm.epsize=epsize;

  pthread_t thread[nThreads];
  //  double qval =0;
  if(nThreads==1)
    for(int i=0;i<nChr;i++){
      objs[i]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.theta,shmm.rho);
  }else {
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

#if 1
  double fwllh,bwllh,qval;
  fwllh=bwllh=qval=0;
  for(int i=0;i<nChr;i++){
    //    fprintf(stderr,"\t-> hmm.fwllh for chr:%d\n",i);
    fwllh += objs[i]->fwllh;
    bwllh += objs[i]->bwllh;
    qval += objs[i]->qval;
  }
  fprintf(stderr,"\t[total llh]  fwllh:%f\n\t[total llh]  bwllh:%f\n\t[total qval] qval:%f\n",fwllh,bwllh,qval);
  ret_llh=fwllh;
  ret_qval=qval;
#endif
}


void main_analysis_optim(double *tk,int tk_l,double *epsize,double theta,double rho){

  shmm.tk=tk;
  shmm.tk_l=tk_l;
  shmm.theta=theta;
  shmm.rho=rho;
  shmm.epsize=epsize;

  pthread_t thread[nThreads];
  if(nThreads==1)
    for(int i=0;i<nChr;i++)
      objs[i]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.theta,shmm.rho);
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
  double fwllh,bwllh,qval;
  fwllh=bwllh=qval=0;
  for(int i=0;i<nChr;i++){
    //    fprintf(stderr,"\t-> hmm.fwllh for chr:%d\n",i);
    fwllh += objs[i]->fwllh;
    bwllh += objs[i]->bwllh;
    qval += objs[i]->qval;
  }
  fprintf(stderr,"\t[total llh]  fwllh:%f\n\t[total llh]  bwllh:%f\n\t[total qval] qval:%f\n",fwllh,bwllh,qval);

}

void smartsize(fastPSMC **myobjs,double *tk,int tk_l,double rho){
  fprintf(stderr,"[smartsize] start;\n");
  double **newepsize = new double*[nChr];
  for(int i=0;i<nChr;i++){
    newepsize[i] = new double[tk_l];
    smartsize1(tk_l,myobjs[i]->windows.size(),myobjs[i]->fw,myobjs[i]->bw,tk,myobjs[i]->P[1],myobjs[i]->emis,myobjs[i]->pix,rho,newepsize[i]);
  }

  fprintf(stderr,"[smartsize] done;\n");
}


void calculate_emissions(double *tk,int tk_l,double *gls,std::vector<wins> &windows,double theta,double **emis,double *epsize);
void main_analysis(double *tk,int tk_l,double *epsize,double theta,double rho,psmc_par *pp,int nIter,int doSmartsize,FILE *FRES,FILE *FLOG){
  assert(FLOG!=NULL);
  make_remapper(pp);
  //test fix:
  double ret_llh,ret_qval,ret_qval2;
#if 1
  fprintf(FLOG,"FLOGS:%p\n",FLOG);
  fprintf(FLOG,"ninter:%d\n",nIter);
  fprintf(FLOG,"domartsiez:%d\n",doSmartsize);
  fprintf(FLOG,"theta:%f\n",theta);
  fprintf(FLOG,"rho:%f\n",rho);
  fprintf(FLOG,"\t-> [%s]\t-> nIter:%d dosmartsize:%d theta:%f rho:%f\n",__FUNCTION__,nIter,doSmartsize,theta,rho);
  for(int i=0;i<tk_l;i++)
    fprintf(FLOG,"\t-> [%s]\t->\t%f\t%f\n",__FUNCTION__,tk[i],epsize[i]);
#endif
  //first make_hmm for all chrs;
  fprintf(FLOG,"[%s]\t-> nIter:%d dosmartsize:%d theta:%f rho:%f\n",__FUNCTION__,nIter,doSmartsize,theta,rho);
  //  rho=0.1;
  double dummyepsize[tk_l];
  for(int i=0;i<tk_l;i++)
    dummyepsize[i] = 1.0;
  
  int at_it=0;
  extern int SIG_COND;
  while(SIG_COND) {
    fprintf(stderr,"----------------------------------------\n");
    fprintf(FLOG,"\t-> Running analysis, RD:%d rho:%f theta:%f\n",at_it,rho,theta);
#if 0
    for(int ii=0;ii<tk_l;ii++) 
      fprintfFLOG,"\t[%d]\tmaking hmm with epsize:%d) %f %f\n",i,ii,tk[ii],epsize[ii]);
#endif
  fflush(FLOG);
    main_analysis_make_hmm(tk,tk_l,epsize,theta,rho,ret_llh,ret_qval);
    if(ncals>0)
      fprintf(stdout,"IT\t%d\n",ncals);
    fprintf(stdout,"RD\t%d\n",at_it);
    fprintf(stdout,"LK\t%f\n",ret_llh);
    fprintf(stdout,"QD\t%f -> %f\n",ret_qval,ret_qval2);
    fprintf(stdout,"RI\t?\n");
    fprintf(stdout,"TR\t%f\t%f\n",theta,rho);
    fprintf(stdout,"MT\t1000000.0\n");
    for(int i=0;i<tk_l;i++)
      fprintf(stdout,"RS\t%d\t%f\t%f\t1000000.0\t1000000.0\t1000000.0\n",i,tk[i],epsize[i]);
    fprintf(stdout,"PA\t%s %.9f %.9f 666.666666666",pp->pattern,theta,rho);
    int at=0;
    for(int i=0;i<remap_l;i++){
      fprintf(stdout," %.9f",epsize[at+remap[i]-1]);
      at+=remap[i];
    }
    fprintf(stdout,"//\n");
    fflush(stdout);

    if(at_it++>=nIter){
      fprintf(stderr,"\t-> Breaking since i>nIter\n");
      break;
    }
    if(doSmartsize==0)
      runoptim3(tk,tk_l,epsize,theta,rho,pp,FLOG,ret_qval2);
    else
      smartsize(objs,tk,tk_l,rho);

    
}

  for(int i=0;i<tk_l;i++) 
    fprintf(stderr,"epsize_after_all_rounds:\t%d\t%f\t%f\n",i,tk[i],epsize[i]);
  
}

int psmc_wrapper(args *pars,int blocksize) {
  fprintf(stderr,"\t-> we are in file: %s function: %s line:%d blocksize:%d\n",__FILE__,__FUNCTION__,__LINE__,blocksize);
  psmc_par *p=pars->par;
#if 1 //print pars
  fprintf(stderr,"\t-> par->n:%d\tpar->n_free:%d\tpar_map:%p\tpar->pattern:%s\tpar->times:%p\tpar->params:%p\n",p->n,p->n_free,p->par_map,p->pattern,p->times,p->params);
  for(int i=0;0&&i<pars->par->n+1;i++)
    fprintf(stderr,"[psmc_wrapper]:%i)\t%f\t%f\n",i,pars->par->times[i],pars->par->params[i]);
  //  exit(0);
#endif
  int tk_l = pars->par->n+1;
  fprintf(stderr,"tk_l in psmc_wrapper pars->par->n+1 tk_l:%d\n",tk_l);
  double *tk = new double [tk_l];
  double *epsize = new double [tk_l];
  setEPSize(epsize,tk_l,p->params);
  //(nelems,array,max_t,alpha,array with values from file, can be NULL)
  setTk(tk_l,tk,15,0.01,p->times);//<- last position will be infinity
  //  fprintf(stderr,"[%s] tk=(%f,%f)\n",__FUNCTION__,tk[0],tk[1]);//exit(0);
#if 0
  for(int i=0;i<tk_l;i++)
    fprintf(stderr,"psmc_wrapper: (tk,epsize)[%d]:(%f,%f)\n",i,tk[i],epsize[i]);
#endif
  
  //initialize all hmm (one for each chr), for now just a single
  int nobs = pars->chooseChr?1:pars->perc->mm.size();
  fprintf(stderr,"\t-> nobs/nchr: %d\n",nobs);
  objs = new fastPSMC*[nobs];
  ops = new oPars[nobs];
  for (myMap::const_iterator it = pars->perc->mm.begin() ;it!=pars->perc->mm.end();it++) {
    myMap::const_iterator it2;
    if(pars->chooseChr!=NULL)
      it2 = iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop);
    else
      it2 = iter_init(pars->perc,it->first,pars->start,pars->stop);
    //    fprintf(stderr,"\t-> Parsing chr:%s \n",it2->first);
    fastPSMC *obj=objs[nChr++]=new fastPSMC;
    //    fprintf(stderr,"gls1:%f %f %f %f\n",pars->perc->gls[0],pars->perc->gls[1],pars->perc->gls[2],pars->perc->gls[3]);
    obj->setWindows(pars->perc->gls,pars->perc->pos,pars->perc->last,pars->blocksize);
    //fprintf(stderr,"gls2:%f %f %f %f\n",pars->perc->gls[0],pars->perc->gls[1],pars->perc->gls[2],pars->perc->gls[3]);
    //  obj->printWindows(stdout);exit(0);
    obj->allocate(tk_l);
    if(pars->chooseChr!=NULL)
      break;
  }
  objs[0]->outnames = strdup(pars->outname);
  assert(pars->flog!=NULL);
  main_analysis(tk,tk_l,epsize,pars->par->TR[0],pars->par->TR[1],pars->par,pars->nIter,pars->smartsize,pars->fres,pars->flog);
  fprintf(pars->fres,"%f\t%f\n",pars->par->TR[0],pars->par->TR[1]);
  for(int i=0;i<tk_l;i++)
    fprintf(pars->fres,"%f\t%f\n",tk[i],epsize[i]);
  return 1;
}
