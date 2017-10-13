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

const double rho = 0.1;

typedef struct{
  double *tk;
  int tk_l;
  double *epsize;
  double rho;
}shared_forhmm;


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



shared_forhmm shmm;


fastPSMC **objs = NULL;
oPars *ops = NULL;

/*
  objective function. Function to be optimized, for each chromo
*/

double qFunction_inner(double *tk,int tk_l,double *epsize,double rho,double pix,int numWind,double **nP,double **PP);


double qFunction(const double *params ,const void *d){
  fprintf(stderr,"asdfasdfadsfasdfdsfdsfasdfdsafadsf\n");
  oPars *data = (oPars*) d;

  return qFunction_inner(data->tk,data->tk_l,data->epsize,data->rho,data->pix,data->numWind,data->nP,data->PP);

}

double qFunction_wrapper(const double *pars,const void *){
  fprintf(stderr,"calling objective function\n");
  for(int i=0;i<nChr;i++){
    /*
      fprintf(stderr,"tk_l:%d\n",ops[i].tk_l);
      fprintf(stderr,".epsize[i]:%f\n",ops[i].epsize[0]);
      fprintf(stderr,".parsze[i]:%f\n",pars[0]);
    */
    for(int j=0;j<ops[i].tk_l;j++){
      ops[i].epsize[j] = pars[j];
    }
  }
  for(int j=0;j<ops[0].tk_l;j++)
    fprintf(stderr," %f",ops[0].epsize[j]);
  fprintf(stderr,"\n");
  for(int i=0;i<nChr;i++)
    qFunction(pars,&ops[i]);
  
  double ret =0;
  for(int i=0;i<nChr;i++)
    ret += objs[i]->qval;
  fprintf(stderr,"qfun:%f\n",ret);
  exit(0);
  return -ret;
}


void runoptim2(double *tk,int tk_l,double *epsize,double rho){
  double **nP = new double *[8];
  for(int i=0;i<8;i++)
    nP[i] = new double[tk_l];

  //get start
  double pars[tk_l];
  for(int i=0;i<tk_l;i++)
    pars[i] = drand48()*5;
  //set bounds
  int nbd[tk_l];
  double lbd[tk_l];
  double ubd[tk_l];
  for(int i=0;i<tk_l;i++){
    nbd[i]=1;
    lbd[i]=0.000001;
    ubd[i]=PSMC_T_INF;
  }

  for(int i=0;i<nChr;i++){
    ops[i].nP = objs[i]->nP;
    ops[i].PP = objs[i]->PP;
    ops[i].tk = tk;
    ops[i].tk_l = tk_l;
    ops[i].pix = objs[i]->pix;
    ops[i].numWind=objs[i]->windows.size();
    ops[i].rho= rho;		// 
    ops[i].epsize=new double[tk_l];
  }
  double max_llh = findmax_bfgs(tk_l,pars,NULL,qFunction_wrapper,NULL,lbd,ubd,nbd,1);
}

void printarray(FILE *fp,double *ary,int l);
void printmatrix(FILE *fp,double **mat,int x,int y);



void printmatrixf(char *fname,double **m,int x,int y){
  return ;
  FILE *fp = NULL;
  if(!(fp=fopen(fname,"wb"))){
    fprintf(stderr,"\t-> Problem writing file: \'%s\'\n",fname);
    exit(0);
  }
  printmatrix(fp,m,x,y);
  fclose(fp);
}

void printarrayf(char *fname,double *m,int x){
  return;
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
  objs[at]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.rho);
  pthread_exit(NULL);
}


void main_analysis_make_hmm(double *tk,int tk_l,double *epsize,double rho){

  fprintf(stderr,"\t-> [%s:%s:%d] nthreads:%d\n",__FILE__,__FUNCTION__,__LINE__,nThreads);
  shmm.tk=tk;
  shmm.tk_l=tk_l;
  shmm.rho=rho;
  shmm.epsize=epsize;

  pthread_t thread[nThreads];
  double qval =0;
  if(nThreads==1) {
    for(int i=0;i<nChr;i++){
      //      fprintf(stderr,"i:%d\n",i);
      qval += objs[i]->make_hmm(shmm.tk,shmm.tk_l,shmm.epsize,shmm.rho);
    }
    
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
  printarrayf((char*)"tk",tk,tk_l);
  printmatrixf((char*)"fw",objs[0]->fw,tk_l,objs[0]->windows.size()+1);
  printmatrixf((char*)"bw",objs[0]->bw,tk_l,objs[0]->windows.size()+1);
  printmatrixf((char*)"emis",objs[0]->emis,tk_l,objs[0]->windows.size()+1);
  printmatrixf((char*)"P",objs[0]->P,7,tk_l);
#endif

#if 1
  double fwllh,bwllh,qval2;
  fwllh=bwllh=qval2=0;
  for(int i=0;i<nChr;i++){
    //    fprintf(stderr,"\t-> hmm.fwllh for chr:%d\n",i);
    fwllh += objs[i]->fwllh;
    bwllh += objs[i]->bwllh;
    qval2 += objs[i]->qval;
  }
  fprintf(stderr,"\t[total llh]  fwllh:%f\n\t[total llh]  bwllh:%f\n\t[total qval] qval:%f\n",fwllh,bwllh,qval2);
#endif
}


void main_analysis_optim(double *tk,int tk_l,double *epsize,double rho){

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




void main_analysis(double *tk,int tk_l,double *epsize,double rho){

  //first make_hmm for all chrs;
  main_analysis_make_hmm(tk,tk_l,epsize,rho);
 
  //  runoptim2(tk,tk_l,epsize,rho);
  for(int i=0;i<20;i++){
    


  }


}

int psmc_wrapper(args *pars,int block) {
  fprintf(stderr,"\t-> we are in file: %s function: %s line:%d\n",__FILE__,__FUNCTION__,__LINE__);
  psmc_par *p=pars->par;
#if 1 //print pars
  fprintf(stderr,"par->n:%d\tpar->n_free:%d\tpar_map:%p\tpar->pattern:%s\tpar->times:%p\tpar->params:%p\n",p->n,p->n_free,p->par_map,p->pattern,p->times,p->params);
  for(int i=0;i<pars->par->n+1;i++)
    fprintf(stderr,"%i)\t%e\t%e\n",i,pars->par->times[i],pars->par->params[i]);
  //  exit(0);
#endif
  int tk_l = pars->par->n+1;
  double *tk = new double [tk_l];
  double *epsize = new double [tk_l];
  setEPSize(epsize,tk_l,p->params);
  //(nelems,array,max_t,alpha,array with values from file, can be NULL)
  setTk(tk_l,tk,15,0.01,p->times);//<- last position will be infinity
  //  fprintf(stderr,"[%s] tk=(%f,%f)\n",__FUNCTION__,tk[0],tk[1]);//exit(0);
#if 0
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
#if 0
  if(0) {
    fprintf(stderr,"\t-> We have now allocated hmm's for: %d chromosomes\n",nChr);
    for(int i=0;i<nobs;i++){
      fprintf(stderr,"\t-> make_hmm for chr:%d\n",i);
      objs[i]->make_hmm(tk,tk_l,epsize,rho);
    }
    double fwllh,bwllh;
    fwllh=bwllh=0;
    for(int i=0;i<nobs;i++){
      //    fprintf(stderr,"\t-> hmm.fwllh for chr:%d\n",i);
      fwllh += objs[i]->fwllh;
      bwllh += objs[i]->bwllh;
    }
    fprintf(stderr,"\t[total llh] fwllh:%f\n\t[total llh] bwllh:%f\n",fwllh,bwllh);
    
    double qs =0;
    for(int i=0;i<nobs;i++){
      fprintf(stderr,"%d/%d\n",i,nobs);
      oPars op;
      op.nP = objs[i]->P;
      op.PP = objs[i]->PP;
      op.tk = tk;
      op.tk_l = tk_l;
      op.pix = objs[i]->pix;
      op.numWind = objs[i]->windows.size();
      op.rho = rho ;
      op.epsize = epsize;
      double tmp = qFunction(NULL,&op);
      fprintf(stderr,"\t -> valQ[%d]: %f\n",i,tmp);
    }
  }else{
#endif
    main_analysis(tk,tk_l,epsize,rho);

  
  
#if 0
  printarrayf("tk",tk,tk_l);
  printmatrixf("fw",objs[0]->fw,tk_l,objs[0]->windows.size()+1);
  printmatrixf("bw",objs[0]->bw,tk_l,objs[0]->windows.size()+1);
  printmatrixf("emis",objs[0]->emis,tk_l,objs[0]->windows.size()+1);
  printmatrixf("P",obj.P,7,tk_l);
#endif
  //printmatrixf("pp",objs[0]->pp,tk_l,objs[0]->windows.size()+1);
  // 
  /*
    printarrayf("stationary",obj.stationary,tk_l);
    
    printarrayf("epsize",epsize,tk_l);

    printmatrixf("emis",obj.emis,tk_l,obj.
windows.size()+1);

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
