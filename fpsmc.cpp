#include <vector>
#include <cassert>
#include <cmath>
#include <pthread.h>
#include <errno.h>
#include "main_psmc.h"
#include "psmcreader.h"
#include "msArg_toPars.h"
#include "hmm_psmc.h"
#include "bfgs.h"
#include "compute.h"
#include "fpsmc.h"
#include "splineEPsize.h"
//#include 
extern int nThreads;

int nChr = 0;

int doQuadratic = 1; //<-only used in qFunction_wrapper

int DOSPLINE=0;

splineEPSize *spl=NULL;

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
  double llh;
  double *parsIn;
}oPars;

int remap_l;
int  *remap;
shared_forhmm shmm;


fastPSMC **objs = NULL;
oPars *ops = NULL;



/*
  This functions either set the tk, NOT the intervals.
  n, is the true length of tk. First entry zero, last entry INF
 */
//remember to subract by one
void setTk(int n, double *t, double max_t, double alpha, double *inp_ti){
  //  assert(inp_ti!=NULL);
  assert(max_t!=-1);
  fprintf(stderr,"\t[%s] (n,tk,max_t,alpha,inp_ti)=(%d,%p,%f,%f,%p)\n",__FUNCTION__,n,t,max_t,alpha,inp_ti);
  int k;
  if (inp_ti == 0) {
    double beta;
    beta = log(1.0 + max_t / alpha) / n; // beta controls the sizes of intervals
    for (k = 0; k < n; ++k)
      t[k] = alpha * (exp(beta * k) - 1);
    t[n] = max_t;
    //    t[n] = PSMC_T_INF; // the infinity: exp(PSMC_T_INF) > 1e310 = inf
  } else {
    memcpy(t, inp_ti, (n+1) * sizeof(double));
  }
}
void setEPSize(double *ary,int tk_l,double *from_infile){
  //  fprintf(stderr,"tk_l:%d\n",tk_l);
  if(!from_infile)
    for (int i = 0; i <tk_l; i++)
      ary[i]=1;
  else{
    memcpy(ary,from_infile,tk_l*sizeof(double));
  }
  for(int i=0;0&&i<tk_l;i++)
    fprintf(stderr,"%d) %f\n",i,from_infile[i]);
  //  exit(0);
}



/*
  objective function. Function to be optimized, for each chromo
*/

double qFunction_inner(double *tk,int tk_l,const double *epsize,double rho,double pix,int numWind,double **nP,double **PP);


double qFunction(const double *params ,const void *d){
  oPars *data = (oPars*) d;
  return qFunction_inner(data->tk,data->tk_l,params,data->rho,data->pix,data->numWind,data->nP,data->PP);
}


double qFunction_inner2(double *tk,int tk_l,const double *epsize,double rho,double pix,int numWind,double **nP,double **PP,double **trans);


double qFunction2(const double *params ,const void *d){
  oPars *data = (oPars*) d;
  //  fprintf(stderr,"pix:%f\n",data->pix);
  return qFunction_inner2(data->tk,data->tk_l,params,data->rho,data->pix,data->numWind,data->nP,data->baumwelch,data->trans);

}

void convert_pattern(const double *pars,double *pars2,int tofull){
  if(tofull==0){
    int at=0;
    for(int i=0;i<remap_l;i++)
      for(int j=0;j<remap[i];j++)
	pars2[at++] = pars[i]; 
  }else{
    int at=0;
    for(int i=0;i<remap_l;i++){
      pars2[i] = pars[at+remap[i]-1];
      at+=remap[i];
    }
  }
}




static int ncals=0;
double qFunction_wrapper(const double *pars,const void *d){
  ncals++;
  double pars2[ops[0].tk_l];
  if(DOSPLINE==0)
    convert_pattern(pars,pars2,0);
  else{
    spl->convert(pars,pars2,0);
    for(int i=0;i<ops[0].tk_l;i++){
      if(pars2[i]<0)
	 return 1000000000;

    }
  }

  for(int i=0;0&&i<ops[0].tk_l;i++)
    fprintf(stderr,"after scaling:%d %f\n",i,pars2[i]);
  //  exit(0);

  if(doQuadratic){
    oPars *data = (oPars*) d;
    ComputeGlobalProbabilities(ops[0].tk,ops[0].tk_l,ops[0].nP,pars2,ops[0].rho);
    double calc_trans(int,int,double**);
    for(int i=0;i<ops[0].tk_l;i++)
      for(int j=0;j<ops[0].tk_l;j++)
	objs[0]->trans[i][j] = calc_trans(i,j,ops[0].nP);
  }
  //fprintf(stderr,"\t-> calling objective function: remap_l:%d [%d]\n",remap_l,ncals);
  double ret =0;
  
  for(int i=0;i<nChr;i++)
    if(doQuadratic)
      ret += qFunction2(pars2,&ops[i]);
    else
      ret += qFunction(pars2,&ops[i]);
  
  if(std::isinf(ret))
    ret= -1000000000;
  
  return -ret;
}


void make_remapper(psmc_par *pp){
  fprintf(stderr,"\t-> [%s] pp->n_free:%d\n",__FUNCTION__,pp->n_free);
  int at=0;
  remap = new int[pp->n_free];
  remap_l=pp->n_free;
  for(int i=0;i<pp->n_free;i++){
    int howMany=0;
    while(pp->par_map[at+howMany]==pp->par_map[at]){
      howMany++;
      if(at+howMany==pp->n+1)//important for breaking properly
	break;
    }
    at+=howMany;
    remap[i] = howMany;
  }
  
#if 0
  for(int i=0;i<pp->n_free;i++)
    fprintf(stderr,"[%d]: %d\n",i,remap[i]);
#endif
  
}

timer starttimer(){
  timer t;
  t.t=clock();
  t.t2=time(NULL);
  return t;
}

void stoptimer(timer &t){
  t.tids[0]= ((float)(clock() - t.t) / CLOCKS_PER_SEC)/60.0;
  t.tids[1]= ((float)(time(NULL) - t.t2))/60.0;
}

//tk is full
void runoptim3(double *tk,int tk_l,double *epsize,double theta,double rho,int ndim,double &ret_qval){
  clock_t t=clock();
  time_t t2=time(NULL);
  fprintf(stderr,"\t-> Starting Optimization ndim:%d\n",ndim);

  double pars[ndim];
  if(DOSPLINE==0)
    convert_pattern(epsize,pars,1);
  else{
    spl->convert(epsize,pars,1);
  }

  //set bounds
  int nbd[ndim];
  double lbd[ndim];
  double ubd[ndim];
  if(DOSPLINE==0){
    for(int i=0;i<ndim;i++){
      nbd[i]=2;
      lbd[i]=0.000001;
      ubd[i]=1000;//PSMC_T_INF;
    }
  }else{
    int at=0;
    for(int i=0;i<ndim/2;i++){

      nbd[i]=2;
      lbd[i]=0.000001;
      ubd[i]=1000;//PSMC_T_INF;
      //      fprintf(stderr,"fv[%d][%d/%d] bd[%d]:%d:(%f,%f)\n",at++,i,ndim/2,i,nbd[i],lbd[i],ubd[i]);
    }
    for(int i=ndim/2;i<ndim;i++){
      nbd[i]=0;
      lbd[i]=-1000;
      ubd[i]=1000;//PSMC_T_INF;
      //fprintf(stderr,"gv[%d][%d/%d] bd[%d]:%d:(%f,%fd)\n",at++,i,ndim/2,i,nbd[i],lbd[i],ubd[i]);
    }
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
  timer opt_timer = starttimer();
  //we are not optimizing llh but qfunction
  double max_qval = findmax_bfgs(ndim,pars,NULL,qFunction_wrapper,NULL,lbd,ubd,nbd,-1);
  stoptimer(opt_timer);
  fprintf(stdout,"MM\toptimization: (wall(min),cpu(min)):(%f,%f)\n",opt_timer.tids[1],opt_timer.tids[0]);
  ret_qval=max_qval;

  if(DOSPLINE==0)
    convert_pattern(pars,epsize,0);
  else{
    spl->convert(pars,epsize,0);
  }
  
  fprintf(stderr, "\t-> [RUNOPTIM3 TIME]:%s cpu-time used =  %.2f sec \n",__func__, (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t-> [RUNOPTIM3 Time]:%s walltime used =  %.2f sec \n",__func__, (float)(time(NULL) - t2));
 
  fflush(stdout);
  fflush(stderr);
}

void *run_a_hmm(void *ptr){
  size_t at =(size_t) ptr;
  objs[at]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.theta,shmm.rho);
  pthread_exit(NULL);
}

void *qFunction2_thd(void *ptr){
  size_t at =(size_t) ptr;
  pthread_exit(NULL);
}


void *qFunction_thd(void *ptr){
  size_t at =(size_t) ptr;
  ops[at].llh = qFunction(ops[at].parsIn,&ops[at]);
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
  objs[0]->make_hmm_pre(shmm.tk,shmm.tk_l,shmm.epsize,shmm.theta,shmm.rho);
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


//tk_l is dimension of transistionsspace ndim is size of dimension
//tk is tk_l long, epsize is tk_l long
void main_analysis(double *tk,int tk_l,double *epsize,double theta,double rho,char *pattern,int ndim,int nIter,double maxt){
  fprintf(stderr,"%s\n",pattern);//exit(0);
  //test fix:
  double ret_llh,ret_qval,ret_qval2;
  ret_qval=ret_qval2=0;
  
  int at_it=0;
  extern int SIG_COND;
  while(SIG_COND) {
    fprintf(stderr,"----------------------------------------\n");

    timer hmm_t = starttimer();
    main_analysis_make_hmm(tk,tk_l,epsize,theta,rho,ret_llh,ret_qval);
    stoptimer(hmm_t);
    
    if(ncals>0)
      fprintf(stdout,"IT\t%d\n",ncals);
    fprintf(stdout,"RD\t%d\n",at_it);
    fprintf(stdout,"LK\t%f\n",ret_llh);
    fprintf(stdout,"QD\t%f -> %f\n",ret_qval,ret_qval2);
    //    exit(0);
    fprintf(stdout,"RI\t?\n");
    fprintf(stdout,"TR\t%f\t%f\n",theta,rho);
    fprintf(stdout,"MT\t%f\n",maxt);
    fprintf(stdout,"MM\tbuildhmm(wall(min),cpu(min)):(%f,%f) tk_l:%d\n",hmm_t.tids[1],hmm_t.tids[0],tk_l);
    for(int i=0;i<tk_l;i++)//this prints out all
      fprintf(stdout,"RS\t%d\t%f\t%f\t1000000.0\t1000000.0\t1000000.0\n",i,tk[i],epsize[i]);
    fprintf(stdout,"PA\t%s %.9f %.9f 666.666666666",pattern,theta,rho);
    int at=0;
    //
    if(strcmp(pattern,"spline")!=0){//this only prints each from each grouping
      //only when not using splines
      for(int i=0;i<remap_l;i++){
	fprintf(stdout," %.9f",epsize[at+remap[i]-1]);
	at+=remap[i];
      }
    }
    //no need to print more in PA line with using spline
    fprintf(stdout,"\n//\n");
    fflush(stdout);
    
    if(at_it++>=nIter){
      fprintf(stderr,"\t-> Breaking since i>nIter\n");
      break;
    }
    runoptim3(tk,tk_l,epsize,theta,rho,ndim,ret_qval2);
    
  }
}

#define DEFAULT_PATTERN "4+5*3+4"
void setpars( char *fname,psmc_par *pp,int which) ;
int *psmc_parse_pattern(const char *pattern, int *n_free, int *n_pars);
int psmc_wrapper(args *pars,int blocksize) {
  DOSPLINE=pars->dospline;
  if(pars->psmc_infile)
    setpars(pars->psmc_infile,pars->par,pars->RD);

  if(pars->par->pattern==NULL)
    pars->par->pattern = strdup(DEFAULT_PATTERN);

  if(pars->par->pattern!=NULL){
    if(pars->par->par_map)
      free(pars->par->par_map);
    pars->par->par_map = psmc_parse_pattern(pars->par->pattern, &pars->par->n_free, &pars->par->n);
  }

  fprintf(stderr,"\t-> we are in file: %s function: %s line:%d blocksize:%d\n",__FILE__,__FUNCTION__,__LINE__,blocksize);

  fprintf(stderr,"\t-> par->n:%d\tpar->n_free:%d\tpar_map:%p\tpar->pattern:%s\tpar->times:%p\tpar->params:%p\n",pars->par->n,pars->par->n_free,pars->par->par_map,pars->par->pattern,pars->par->times,pars->par->params);

  if(pars->msstr){
    pars->msstr_arg = parse_msStr(pars->msstr); 
    msarg_toPars(pars->msstr_arg,pars->par,pars->perc->version==0?pars->blocksize:1);
  }else
    make_remapper(pars->par);

  //adjust theta:
  pars->par->TR[0] = pars->par->TR[0]/2.0;
  fprintf(stderr,"\t-> p->perc->version:%d (one is gls, otherwise fasta)\n",pars->perc->version);
  if(pars->perc->version==1){//if it is gls
    fprintf(stderr,"\t-> Adjusing theta with blocksize: %d\n",pars->blocksize);
    pars->par->TR[0] = pars->par->TR[0]/(1.0*pars->blocksize);
  }

  
  int tk_l = pars->par->n+1;
  double *tk,*epsize;
  tk=epsize=NULL;
  int ndim=-1;
  char *pattern = NULL;
  double max_t =pars->psmc_infile?pars->par->MT: pars->init_max_t;
  if(DOSPLINE!=0){
    spl = new splineEPSize(14,3,7,pars->init_max_t);
    fprintf(stderr,"\t-> spl.tk_l:%d spl->ndim:%d\n",spl->tk_l,spl->ndim);
    tk_l = spl->tk_l;
    tk=spl->tk;
    pattern=strdup("spline");
    ndim=spl->ndim;
    spl->fillit();
    spl->computeSpline();
    epsize = new double [tk_l];
    assert(pars->init==-1);
    spl->computeEPSize(epsize);
    for(int i=0;0&&i<tk_l;i++)
      fprintf(stderr,"%d %f\n",i,epsize[i]);
  }else{
    fprintf(stderr,"\t-> allocating array with length: tk_l:%d\n",tk_l);
    tk = new double [tk_l];
    setTk(tk_l-1,tk,max_t,0.1,pars->par->times);//<- last position will be infinity
    pattern=pars->par->pattern;
    ndim=pars->par->n_free;
    epsize = new double [tk_l];
    setEPSize(epsize,tk_l,pars->par->params);
  }
  fprintf(stderr,"tk_l:%d\n",tk_l);

  double theta=pars->par->TR[0];
  double rho=pars->par->TR[1];

 

  if(pars->init!=-1)
    for(int i=0;i<tk_l;i++)
      epsize[i] = pars->init;
  if(pars->init_theta!=-1)
    theta=pars->init_theta;
  if(pars->init_rho!=-1)
    rho=pars->init_rho;
  
  assert(theta!=-1&&rho!=-1&&tk_l>0);
  
#if 0
  for(int i=0;i<tk_l;i++)
    fprintf(stderr,"psmc_wrapper: (tk,epsize)[%d]:(%f,%f)\n",i,tk[i],epsize[i]);
  exit(0);
#endif
  fprintf(stderr,"\t-> tk_l in psmc_wrapper pars->par->n+1 tk_l:%d p->times:%p\n",tk_l,pars->par->times);  
  int nobs = pars->chooseChr?1:pars->perc->mm.size();
  fprintf(stderr,"\t-> nobs/nchr: %d\n",nobs);
  objs = new fastPSMC*[nobs];
  ops = new oPars[nobs];
  for (myMap::const_iterator it = pars->perc->mm.begin() ;it!=pars->perc->mm.end();it++) {
    myMap::const_iterator it2;
    if(pars->chooseChr!=NULL)
      it2 = iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop,pars->blocksize);
    else
      it2 = iter_init(pars->perc,it->first,pars->start,pars->stop,pars->blocksize);
    //    fprintf(stderr,"\t-> Parsing chr:%s \n",it2->first);
    fastPSMC *obj=objs[nChr++]=new fastPSMC;
    //    fprintf(stderr,"gls1:%f %f %f %f\n",pars->perc->gls[0],pars->perc->gls[1],pars->perc->gls[2],pars->perc->gls[3]);
    obj->setWindows(pars->perc->gls,pars->perc->pos,pars->perc->last,pars->blocksize);
    //fprintf(stderr,"gls2:%f %f %f %f\n",pars->perc->gls[0],pars->perc->gls[1],pars->perc->gls[2],pars->perc->gls[3]);

    obj->allocate(tk_l);
    //    fprintf(stderr,"transer:%p\n",obj[0].trans);
    if(pars->chooseChr!=NULL)
      break;
  }
  objs[0]->outnames = strdup(pars->outname);
  main_analysis(tk,tk_l,epsize,theta,rho,pattern,ndim,pars->nIter,max_t);

  free(objs[0]->outnames);
  for (int i=0;i<nChr;i++)
    delete objs[i];
  delete [] objs;
  delete [] remap;
  delete [] tk;
  delete [] epsize;
  delete [] ops;
  return 1;
}
