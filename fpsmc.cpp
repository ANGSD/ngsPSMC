#include <vector>
#include <cassert>
#include <cmath>
#include <pthread.h>
#include "psmcreader.h"
#include "main_psmc.h"
#include "hmm_psmc.h"
#include "bfgs.h"
#include <errno.h>

extern int nThreads;

int nChr = 0;

int doQuadratic = 0;


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

static int ncals=0;
double qFunction_wrapper(const double *pars,const void *){
  ncals++;
  //  fprintf(stderr,"\t-> calling objective function: remap_l:%d\n\n",remap_l);
  double pars2[ops[0].tk_l];
  int at=0;
  for(int i=0;i<remap_l;i++)
    for(int j=0;j<remap[i];j++){
      pars2[at++] = pars[i]; 
      //      fprintf(stderr,"\t-> pars2: %e\n",pars2[at-1]);
    }
  double ret =0;
  for(int i=0;i<nChr;i++)
    if(doQuadratic)
      ret += qFunction2(pars2,&ops[i]);
    else
      ret += qFunction(pars2,&ops[i]);

  if(std::isinf(ret))
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

}


void runoptim3(double *tk,int tk_l,double *epsize,double theta,double rho,psmc_par *pp){
  fprintf(stderr,"STARNG OPTIM\n");
#if 0
  fprintf(stderr,"pp->n:%d\n",pp->n);
  fprintf(stderr,"pp->n_free:%d\n",pp->n_free);
  for(int i=0;i<pp->n+1;i++)
    fprintf(stderr,"pars_map[%d]\t%d\n",i,pp->par_map[i]);
#endif
  make_remapper(pp);
  int ndim = pp->n_free;

  double pars[ndim];
  for(int i=0;0&&i<ndim;i++)
    pars[i] =5.0*drand48()+1e-3;

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
  for(int i=0;i<ndim;i++)
     fprintf(stderr,"pre opt: %f\n",pars[i]);
  double max_llh = findmax_bfgs(ndim,pars,NULL,qFunction_wrapper,NULL,lbd,ubd,nbd,-1);
  for(int i=0;i<ndim;i++)     
    fprintf(stderr,"post opt: %f\n",pars[i]);
  fprintf(stderr,"\t-> optim done: after ncalls:%d best total qval:%f\n",ncals,max_llh);
  for(int i=0;1&i<ndim;i++)
    fprintf(stderr,"newpars after optim: %f\n", pars[i]);
  at=0;
  for(int i=0;i<remap_l;i++)
    for(int j=0;j<remap[i];j++)
      epsize[at++] = pars[i]; 
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


void main_analysis_make_hmm(double *tk,int tk_l,double *epsize,double theta,double rho){

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
void calculate_emissions(double *tk,int tk_l,double *gls,std::vector<wins> &windows,double theta,double **emis,double *epsize);
void main_analysis(double *tk,int tk_l,double *epsize,double theta,double rho,psmc_par *pp,int nIter){
  fprintf(stderr,"\t-> nIter:%d\n",nIter);
  //first make_hmm for all chrs;
#if 1
  theta=0.046797/2.0;
  rho = 0.006664;
#endif
  //  rho=0.1;
  double dummyepsize[tk_l];
  for(int i=0;i<tk_l;i++)
    dummyepsize[i] = 1.0;

  for(int i=0;0&&i<nChr;i++)
    calculate_emissions(tk,tk_l,objs[i]->gls,objs[i]->windows,theta,objs[i]->emis,dummyepsize);

  int i=0;
  while(1){
    fprintf(stderr,"\t-> Running analysis, RD:%d rho:%f theta:%f\n",i,rho,theta);
#if 1
    for(int ii=0;ii<tk_l;ii++) 
      fprintf(stderr,"making hmm with epsize:%d) %f %f\n",ii,tk[ii],epsize[ii]);
#endif
    main_analysis_make_hmm(tk,tk_l,epsize,theta,rho);
    
    if(i++>=nIter)
      break;
    runoptim3(tk,tk_l,epsize,theta,rho,pp);
    //    break;
  }
  
  for(int i=0;i<tk_l;i++) 
    fprintf(stderr,"epsize_after_all_rounds:%d) %f\n",i,epsize[i]);

}

int psmc_wrapper(args *pars,int block) {
  fprintf(stderr,"\t-> we are in file: %s function: %s line:%d\n",__FILE__,__FUNCTION__,__LINE__);
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
#if 1
  for(int i=0;i<tk_l;i++)
    fprintf(stderr,"(tk,epsize)[%d]:(%e,%e)\n",i,tk[i],epsize[i]);
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
    obj->setWindows(pars->perc->gls,pars->perc->pos,pars->perc->last,pars->block);
    //fprintf(stderr,"gls2:%f %f %f %f\n",pars->perc->gls[0],pars->perc->gls[1],pars->perc->gls[2],pars->perc->gls[3]);
    //  obj->printWindows(stdout);exit(0);
    obj->allocate(tk_l);
    if(pars->chooseChr!=NULL)
      break;
  }
  main_analysis(tk,tk_l,epsize,pars->par->TR[0],pars->par->TR[1],pars->par,pars->nIter);
  return 1;
}
