#include <ctime>
#include <vector>
#include <cassert>
#include <cmath>
#include "psmcreader.h"
#include "main_psmc.h"
#include "hmm_psmc.h"
#include "compute.h"

int fastPSMC::tot_index=0;

//#define __SHOW_TIME__


void printmatrixf(char *fname,double **m,int x,int y);
void printarrayf(char *fname,double *m,int x);

void printarray(FILE *fp,double *ary,int l){
  for(int i=0;i<l;i++)
    fprintf(fp,"%f\t",ary[i]);

}
void printmatrix(FILE *fp,double **mat,int x,int y){
  fprintf(fp,"#printmatrix with x:%d y:%d\n",x,y);
  for(int i=0;i<y;i++){
    for(int j=0;j<x-1;j++)
      fprintf(fp,"%e\t",mat[j][i]);
    fprintf(fp,"%e\n",mat[x-1][i]);
  }

}

/*
  k: is the state k
  pix: is the perchromsome value of the forward probs
  numwin: is not being used currently
  nP: is the 'new' Pi's (the ones computed by compuateGlobalProbablies)
  PP: is the old PPi's whic correspond to expectation. 

  returnvalue is in logspace
 */

/*
  if something doesnt work check this formula:

  trans<-(exp(read.table("transitions.txt")));
fw<-t(exp(read.table("fw")))
bw<-t(exp(read.table("bw")))
emis<-t(exp(read.table("emis")))
stat<-exp(scan("stationary"))
P<-exp(read.table("P.txt"));colnames(P)<-paste0("P",0:7);
PP<-exp(read.table("PP.txt"));colnames(PP)<-paste0("PP",0:7);
##emis[,1:10]
##fw[,1:10]



 */

double qkFunction(unsigned k, double pix, unsigned numWind,double **nP,double **PP,int tk_l,double &esum){
  //  fprintf(stderr,"k:%u pix:%f numWind:%u tk_l:%d\n",k,pix,numWind,tk_l);
  /*
  //This block is needed if eimission probabilities depend on estimated parameters, e.g. on time disctretisation 
    for (unsigned l = 1; l < numWind + 1; l++) 
    qi += log(emis[K][l])*fw[K][l]*bw[K][l];
    qi /= pix;
  */
  //  printmatrixf((char*)"nP.txt",nP,8,tk_l);
  double qi[7];
  double expec[7];
  double npfac[7];

  
  //i follows PP index.
  for(int i=1;i<8;i++){
    expec[i-1] = exp(lprod(PP[i][k],-pix));
    npfac[i-1] = nP[i][k];
    if(i!=2&&i!=5)
      esum +=expec[i-1];

    if(std::isinf(npfac[i-1]))
      qi[i-1] = 0;
    else
      qi[i-1] = npfac[i-1]*expec[i-1];
    //fprintf(stderr,"qk[%d]:%f\tnpfac:%f\texpec:%f esum:%f\n",i,qi[i-1],npfac[i-1],expec[i-1],esum);
  }
  
  //  fprintf(stderr,"ESUM:%f\n",esum);

#if 0
  qi[0] = nP[1][k]*exp(lprod(PP[1][k],-pix));
  fprintf(stderr,"qk[0]:%f np[1][k]:%f PP[1][k]:%f pix:%f\n",qi[0],nP[1][k],PP[1][k],pix);
  qi[1] = nP[2][k]*exp(lprod(PP[2][k],-pix));
  fprintf(stderr,"qk[1]:%f\n",qi[1]);
  qi[2] = nP[3][k]*exp(lprod(PP[3][k],-pix));
  fprintf(stderr,"qk[2]:%f\n",qi[2]);
  qi[3] = nP[4][k]*exp(lprod(PP[4][k],-pix));
  fprintf(stderr,"qk[3]:%f\n",qi[3]);
  qi[4] = nP[5][k]*exp(lprod(PP[5][k],-pix));
  fprintf(stderr,"qk[4]:%f\n",qi[4]);
  qi[5] = nP[6][k]*exp(lprod(PP[6][k],-pix));
  fprintf(stderr,"qk[5]:%f\n",qi[5]);
  qi[6] = nP[7][k]*exp(lprod(PP[7][k],-pix));
  fprintf(stderr,"qk[6]:%f\n",qi[6]);
  for(int i=0;(k==tk_l-1)&&i<7;i++){//DRAGON
    //nan occurs if expected value is zero, and nP is -Inf
    if(isinf(qi[i])||isnan(qi[i])){
      qi[i]=0.0;
    }
  }
#endif

  double ret = 0;
  for(int i=0;i<7;i++) 
    //    if(i!=2&&i!=)
    ret += qi[i];

  return ret;
}
/*
  tk: intervals, fixed
  tk_l: length of tk
  P: matrix where results will be plugged in
  epsize: the effective populationsize (lambda)
  rho: rho, or theta to rho value.
 */

void ComputeGlobalProbabilities(double *tk,int tk_l,double **P,double *epsize,double rho){
#ifdef __SHOW_TIME__
    clock_t t=clock();
    time_t t2=time(NULL);
#endif
#if 0
  fprintf(stderr,"[%s] rho:%e\n",__FUNCTION__,rho);
  fprintf(stderr,"[%s] tks[0]:%f tks[1]:%f\n",__FUNCTION__,tk[0],tk[1]);
  fprintf(stderr,"[%s] epsize[0]:%f epsize[1]:%f\n",__FUNCTION__,epsize[0],epsize[1]);
#endif
  ComputeP1(tk,tk_l,P[1],epsize,rho);
  ComputeP5(tk,tk_l,P[5],epsize);
  ComputeP6(tk,tk_l,P[6],epsize,rho);
  ComputeP2(tk_l,P[2],P[5]);
  ComputeP3(tk,tk_l,P[3],epsize,rho);
  ComputeP4(tk,tk_l,P[4],epsize,rho);
  ComputeP7(tk,tk_l,P[7],P[3],epsize,rho);
  ComputeP0(tk_l,P[0],P[5]);
  
  for(int p=0;1&&p<8;p++){
    for(int i=0;i<tk_l;i++){
      //      fprintf(stderr,"P[%d][%d]: %f\n",p,i,P[p][i]);
      assert(exp(P[p][i])>=0&&exp(P[p][i])<=1);
    }

  }
#ifdef __SHOW_TIME__
    fprintf(stderr, "\t[TIME] cpu-time used =  %.2f sec for computeglobalprobabilites\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(stderr, "\t[Time] walltime used =  %.2f sec for computeglobalprobablibies\n", (float)(time(NULL) - t2));  
#endif
}

double qFunction_inner(double *tk,int tk_l,double *epsize,double rho,double pix,int numWind,double **nP,double **PP){
#ifdef __SHOW_TIME__
  clock_t t=clock();
  time_t t2=time(NULL);
#endif
  
  ComputeGlobalProbabilities(tk,tk_l,nP,epsize,rho);
  //  printmatrixf((char*)"P_check.txt",nP,8,tk_l);
  double Q = 0;
  double esum =0;
  for (unsigned i = 0; i < tk_l; i++){
    double tmpQ = qkFunction(i, pix,numWind,nP,PP,tk_l,esum);
    Q += tmpQ;
    //    fprintf(stderr,"\t-> Q[%d]:%f\n",i,tmpQ);
  }
  fprintf(stderr,"\t-> ESUM:%f\n",esum);
#ifdef __SHOW_TIME__
   fprintf(stderr, "\t[TIME]:%s cpu-time used =  %.2f sec \n",__func__, (float)(clock() - t) / CLOCKS_PER_SEC);
   fprintf(stderr, "\t[Time]:%s walltime used =  %.2f sec \n",__func__, (float)(time(NULL) - t2)); 
#endif
  return Q;
  
}

/*
  Functions below are template for future ...


 */
void ComputeSW(int maxTime,double W[],double sW[]){
  double tmp = 0;
  for (int i = maxTime; i >=0 ; i--){
    tmp += W[i];
    sW[i] = tmp;
  }
}
void UpdateEPSize(int maxTime, double W[],double sW[],double epSize[],double T[]){
  for (unsigned i = 1; i < maxTime; i++){
    double tmp;
    tmp = W[i]/sW[i];
    epSize[i] = -log(1 - tmp)/(T[i+1]-T[i]);
  }
}




/*
  This functions either set the tk, NOT the intervals.
  n, is the true length of tk. First entry zero, last entry INF
 */
void setTk(int n, double *t, double max_t, double alpha, double *inp_ti){
  //  fprintf(stderr,"[%s] (n,tk,max_t,alpha,inp_ti)=(%d,%p,%f,%f,%p)\n",__FUNCTION__,n,t,max_t,alpha,inp_ti);
  int k;
  if (inp_ti == 0) {
    double beta;
    beta = log(1.0 + max_t / alpha) / n; // beta controls the sizes of intervals
    for (k = 0; k < n; ++k)
      t[k] = alpha * (exp(beta * k) - 1);
    t[n-1] = max_t;
    t[n] = PSMC_T_INF; // the infinity: exp(PSMC_T_INF) > 1e310 = inf
  } else {
    memcpy(t, inp_ti, n * sizeof(double));
  }
}


void setEPSize(double *ary,int l,double *from_infile){
  // fprintf(stderr,"l:%d\n",l);
  if(!from_infile)
    for (int i = 0; i <l; i++)
      ary[l]=1;
  else{
    memcpy(ary,from_infile,l*sizeof(double));
  }
  for(int i=0;0&&i<l;i++)
    fprintf(stderr,"%d) %f\n",i,from_infile[i]);
  //  exit(0);
}


void fastPSMC::calculate_FW_BW_PP_Probs(double *tk,int tk_l,double *epsize,double rho){
#ifdef __SHOW_TIME__
  clock_t t=clock();
  time_t t2=time(NULL);
#endif
    //we first set the initial fwprobs to stationary distribution
  for(int i=0;i<tk_l;i++){
      fw[i][0] = stationary[i];
      //      fprintf(stderr,"fw[%d]:stat[%d]:%f\n",i,i,fw[i][0]);
  }
  //  double bwPrime[tk_l];

    //we now loop over windows.
    //v=0 is above and is the initial distribution, we therefore plug in at v+1
    for(int v=0;v<windows.size();v++){
      ComputeRs(v,fw);//<-prepare R1,R2
#if 0
      printarrayf("r1",R1,tk_l);
      printarrayf("r2",R2,tk_l);
      // exit(0);
#endif
      fw[0][v+1] = addProtect3(lprod(fw[0][v],P[1][0]) , lprod(R1[0],P[3][0]) , lprod(fw[0][v],P[4][0]))+emis[0][v+1] ;
      for (unsigned i = 1; i < tk_l; i++)
	fw[i][v+1]= addProtect4(lprod(fw[i][v],P[1][i]) , lprod(R2[i-1],P[2][i]) , lprod(R1[i],P[3][i]) , lprod(fw[i][v],P[4][i]))+emis[i][v+1];
    }

    //    fprintf(stderr,"tk_l:%d\n",tk_l);
    double tmp[tk_l];
    for(int i=0;i<tk_l;i++){
      tmp[i] = fw[i][windows.size()];
    }
    pix = addProtectN(tmp,tk_l);
    // fprintf(stderr,"pix:%f\n",pix);exit(0);
    assert(!isnan(pix));
    fwllh = pix;
    //fprintf(stderr,"forward(pic) llh:%f\n",pix);


    //now do backward algorithm
    //initialize by stationary
    for(int i=0;i<tk_l;i++)
      bw[i][windows.size()] = lprod(stationary[i],emis[i][windows.size()]);
     

    //we plug in values at v-1, therefore we break at v==1
    for(int v=windows.size();v>0;v--){
      ComputeRs(v,bw);//<-prepare R1,R2
#if 0
      double p1= lprod(bw[0][v],P[1][0]);
      double p2= lprod( R1[0],P[3][0]) ;
      double p3= lprod(bw[0][v],P[4][0]);
#endif
      //      fprintf(stderr,"p1:%f\tp2:%f\tp3:%f emis:%f\n",p1,p2,p3,emis[0][v]);
      bw[0][v-1] = addProtect3(lprod(bw[0][v],P[1][0]) , lprod( R1[0],P[3][0]) , lprod(bw[0][v],P[4][0]))+emis[0][v-1];
      bw[0][v] -= lprod(stationary[0],emis[0][v]);
      //fprintf(stderr,"bw[[0][%d]:%f\n",v,bw[0][v]);
      for (unsigned i = 1; i < tk_l; i++){
	bw[i][v-1] = addProtect4(lprod(bw[i][v],P[1][i]),lprod(R2[i-1],P[2][i]),lprod(R1[i],P[3][i]),lprod(bw[i][v],P[4][i]))+emis[i][v-1];
	bw[i][v] -= lprod(stationary[i],emis[i][v]);
	//	fprintf(stderr,"bw[[%d][%d]:%f stationary[i]:%d emis[i][v]:%f\n",i,v,bw[0][v],emis[i][v]);
      }
    }

    for(int i=0;i<tk_l;i++)
      bw[i][0] -= stationary[i];
     

#if 0
    for(int v=1;v<=windows.size();v++)
      for(int i=0;i<tk_l;i++)
	bw[i][v] =0;
#endif
    
    for(int i=0;i<tk_l;i++)
      tmp[i] = bw[i][1]+stationary[i]+emis[i][1];
    double tmptmp= addProtectN(tmp,tk_l);
    assert(!isnan(tmptmp));
    bwllh = tmptmp;
    //    fprintf(stderr,"backward llh:%f\n",tmptmp);

    //calculate post prob per window per state
#if 0
    for(int v=1;v<windows.size();v++){
      for(int j=0;j<tk_l;j++)
	tmp[j] = fw[j][v]+bw[j][v];
      tmptmp= addProtectN(tmp,tk_l);
      assert(!isnan(tmptmp));
    }
#endif
    //fprintf(stderr,"[%s] stop\n",__FUNCTION__ );

;
#ifdef __SHOW_TIME__
 fprintf(stderr, "\t[TIME]:%s cpu-time used =  %.2f sec \n",__func__, (float)(clock() - t) / CLOCKS_PER_SEC);
 fprintf(stderr, "\t[Time]:%s walltime used =  %.2f sec \n",__func__, (float)(time(NULL) - t2));  
#endif
}



void fastPSMC::allocate(int tk_l_arg){
#ifdef __SHOW_TIME__
  clock_t t=clock();
  time_t t2=time(NULL);
#endif
  
  int numWindows = windows.size();
  //  fprintf(stderr,"\t-> [%s]: will allocate tk with length: %d\n",__FUNCTION__,tk_l_arg);
  tk_l = tk_l_arg;
  //tk = new double[tk_l];
 
  //  printarray(stderr,tk,tk_l);
  stationary = new double[tk_l];
  R1 = new double[tk_l];
  R2 = new double[tk_l];
  fw = new double *[tk_l];
  bw = new double *[tk_l];
  //pp = new double *[tk_l];
  emis = new double *[tk_l];
  for(int i=0;i<tk_l;i++){
    emis[i] = new double[numWindows+1];
#if 0
    for(int j=0;j<numWindows+1;j++)
      emis[i][j] = 0;//1.0;
#endif
    fw[i] = new double[numWindows+1];
    bw[i] = new double[numWindows+1];
    //    pp[i] = new double[numWindows+1];

  }
  //  fprintf(stderr,"\t-> emission allocated with [%d][%d]\n",tk_l,numWindows+1);
  P = new double *[8];
  PP= new double *[8];
  nP = new double*[8];
  for(int i=0;i<8;i++){
    P[i] = new double[tk_l];
    PP[i]= new double[tk_l];
    nP[i]= new double[tk_l];
  }
  for(int i=0;i<tk_l;i++)
    PP[0][i] = -666;
  workspace = new double[windows.size()];
  if(DOTRANS){
    trans = new double *[tk_l];
    for(int i=0;i<tk_l;i++){
      trans[i] = new double[tk_l];
      for(int j=0;j<tk_l;j++)
	trans[i][j] = -888;//placeholder, to spot if something shouldnt be happening;
    }
    //    printmatrixf("TRANS",trans,tk_l,tk_l);
  }
#ifdef __SHOW_TIME__
  fprintf(stderr, "\t[TIME] cpu-time used =  %.2f sec for allocating internal structures\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[Time] walltime used =  %.2f sec for  allocating internal structures\n", (float)(time(NULL) - t2));  
#endif
}
/*
  Function will set the indices for the windows
  first index and last index INCLUSIVE
 */
//void fastPSMC::setWindows(perpsmc *perc,char *chooseChr,int start,int stop,int block){
void fastPSMC::setWindows(double *gls_a,int *pos ,int last,int block){
  gls = new double[last*2];
  memcpy(gls,gls_a,last*sizeof(double)*2);
  //  fprintf(stderr,"\t-> [%s] gls:(%f, %f,%f, %f) \n",__FUNCTION__,gls_a[0],gls_a[1],gls_a[2],gls_a[3] );
  int beginIndex =0;
  int endIndex=0;
  int beginPos = 0;
  int endPos = beginPos+block-1;
  
  while(1) {
    //    fprintf(stderr,"Beginpos:%d EndPos:%d\n",beginPos,endPos);
    wins w;
    if(endPos>pos[last-1])
      break;
    
    while(pos[beginIndex]<beginPos)
      beginIndex++;
    while(pos[endIndex]<endPos)
      endIndex++;
    endIndex--;
#if 0
    //fprintf(stdout,"\t-> endpos:%d\n",pos[endIndex]);
    fprintf(stdout,"\t-> winsize:%d bp:%d,ep:%d bi:%d ei:%d ei-bi:%d\n",block,beginPos,endPos,beginIndex,endIndex,endIndex-beginIndex);
#endif
    w.from = beginIndex;
    w.to = endIndex;
    windows.push_back(w);
    beginPos+=block;
    endPos+=block;
    //    exit(0);
  }
  //  exit(0);
}


/* 
  Calculate emission probabilityes
  tk array of length tk_l
  lambda array effective population sizes
  both has length tk_l
  
  emission probablities will be put in the **emis
  
  stationary(i) = exp(-sum_{j=0}^{i-1}{tau_j/lambda_j}*P2[i])
 */
void calculate_emissions(double *tk,int tk_l,double *gls,std::vector<wins> &windows,double theta,double **emis){
#ifdef __SHOW_TIME__
  clock_t t=clock();
  time_t t2=time(NULL);
#endif 

  //  fprintf(stderr,"\t-> [Calculating emissions with tk_l:%d and windows.size():%lu:%s ] theta:%f gls:(%f,%f,%f,%f) start\n",tk_l,windows.size(),__TIME__,theta,gls[0],gls[1],gls[2],gls[3]);
  //initialize the first:
  for(int j=0;j<tk_l;j++)
    emis[j][0] = log(0);
 
  //  double tmp[windows.size()];
  for(int v=0;v<windows.size();v++){//for each window
    //    fprintf(stderr,"v:%d from:%d to:%d\n",v,windows[v].from,windows[v].to);
    for(int j=0;j<tk_l;j++){//for each interval/state
      //fprintf(stderr,"\tj:%d\n",j);
      emis[j][v+1] = 0;
      double inner = exp(-2.0*tk[j]*theta); // this part relates to issue #1
      //      fprintf(stderr,"\t\t%d from:%d to:%d inner:%f\n",j,windows[v].from,windows[v].to,inner);
      for(int i=windows[v].from;i<=windows[v].to;i++){//for all elements in window
#if 0
	fprintf(stderr,"\t\t\tgls(%d,%d)=",2*i,2*i+1);
	fprintf(stderr,"(%f,%f)\n",gls[2*i],gls[2*i+1]);
#endif
	emis[j][v+1] += log((exp(gls[i*2])/4.0) *inner + (exp(gls[2*i+1])/6)*(1-inner));//<- check
	if(isinf(emis[j][v+1])){
	  fprintf(stderr,"\t-> Huge bug in code contact developer. Emissions evaluates to zero\n");
	  /*
	    This will corrupt the computation of the backward probability.
	  */
	  exit(0);
	}
      }
    }
  }
  //  exit(0);
  //fprintf(stderr,"\t-> [Calculating emissions with tk_l:%d and windows.size():%lu:%s ] stop\n",tk_l,windows.size(),__TIME__);
#if __SHOW_TIME__
  fprintf(stderr, "\t[TIME]:%s cpu-time used =  %.2f sec \n",__func__, (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[Time]:%s walltime used =  %.2f sec \n",__func__, (float)(time(NULL) - t2));  
#endif
}


/*
  Calculate stationary distrubution
  tk array of length tk_l
  lambda array effective population sizes
  both has length tk_l
  
  stationary distribution will be put in results array, also of length tk
  
  stationary(i) = exp(-sum_{j=0}^{i-1}{tau_j/lambda_j}*P2[void])
*/

void fastPSMC::calculate_stationary(double *tk,int tk_l,double *lambda,double *results,double **P){
  results[0] = P[2][0];//fix this
  for(int i=1;i<tk_l;i++){
    results[i]  = P[2][i]+P[0][i-1];
  }
#if 1 //check it sums to one
  double tmp=0;
  for(int i=0;i<tk_l;i++)
    tmp += exp(results[i]);
  //  fprintf(stderr,"\t-> sum of stationary:%f\n",tmp);
  //  exit(0);

#endif

}

void ComputePii(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary,double **emis,double *workspace){
#ifdef __SHOW_TIME__
  clock_t t=clock();
  time_t t2=time(NULL);
#endif
  //  static double *workspace = new double [numWind];
  ComputeP11(numWind,tk_l,P[1],PP[1],fw,bw,workspace,emis);
  ComputeP22(numWind,tk_l,P,PP[2],fw,bw,emis);
  ComputeP33(numWind,tk_l,P[3],PP[3],fw,bw,emis);
  ComputeP44(numWind,tk_l,P[4],PP[4],fw,bw,workspace,emis);
  ComputeP55(numWind,tk_l,P,PP[5],fw,bw,stationary,emis);
  ComputeP66(numWind,tk_l,P,PP[6],fw,bw,stationary,emis);
  ComputeP77(numWind,tk_l,P,PP[7],fw,bw,stationary,emis);
  
  //printmatrixf((char*)"PP.txt",PP,8,tk_l);
  
  for(int p=1;0&&p<8;p++){//CHECK IF THIS SHOULD BE RENAABLED
    for(int i=0;i<tk_l;i++){
      //      fprintf(stderr,"P[%d][%d]: %f\n",p,i,PP[p][i]);
      assert(exp(PP[p][i])>=0&&exp(PP[p][i])<=1);
    }
  }
#ifdef __SHOW_TIME__
 fprintf(stderr, "\t[TIME]:%s cpu-time used =  %.2f sec \n",__func__, (float)(clock() - t) / CLOCKS_PER_SEC);
 fprintf(stderr, "\t[Time]:%s walltime used =  %.2f sec \n",__func__, (float)(time(NULL) - t2)); 
#endif
  //  exit(0);
}


double fastPSMC::make_hmm(double *tk,int tk_l,double *epsize,double rho){

  fprintf(stderr,"\t-> [%s][%d] tk=(%f,%f) gls:(%f, %f,%f, %f) \n",__FUNCTION__,index,tk[0],tk[1],gls[0],gls[1],gls[2],gls[3] );
  //prepare global probs

  ComputeGlobalProbabilities(tk,tk_l,P,epsize,rho);//only the P* ones
  calculate_emissions(tk,tk_l,gls,windows,theta,emis);
  calculate_stationary(tk,tk_l,epsize,stationary,P);
  calculate_FW_BW_PP_Probs(tk,tk_l,epsize,rho);
  ComputePii(windows.size(),tk_l,P,PP,fw,bw,stationary,emis,workspace);

#if 0
  if(index==0)
    printmatrixf((char*)"P_0.txt",P,8,tk_l);
  else if(index==1)
    printmatrixf((char*)"P_1.txt",P,8,tk_l);
  //calculate emissions¯
  if(index==0)
    printmatrixf((char*)"emis_0",emis,tk_l,windows.size()+1); 
  else if(index==1)
    printmatrixf((char*)"emis_1",emis,tk_l,windows.size()+1); 
  if(index==0)
    printarrayf((char*)"stationary_0",stationary,tk_l);
  else if(index==1)
    printarrayf((char*)"stationary_1",stationary,tk_l);
  if(index==0){
    printmatrixf((char*)"fw_0",fw,tk_l,windows.size()+1); 
    printmatrixf((char*)"bw_0",bw,tk_l,windows.size()+1); 
  }else if(index==1){
    printmatrixf((char*)"fw_1",fw,tk_l,windows.size()+1); 
    printmatrixf((char*)"bw_1",bw,tk_l,windows.size()+1); 
  }
  if(index==0)
    printmatrixf((char*)"PP_0.txt",PP,8,tk_l);
  else if(index==1)
    printmatrixf((char*)"PP_1.txt",PP,8,tk_l);
#endif

  qval=qFunction_inner(tk,tk_l,epsize,rho,pix,windows.size(),P,PP);//no need to recompute P. But we fix this later;
  fprintf(stderr,"\t-> hmm[%d]\tqval: %f fwllh: %f bwllh: %f\n",index,qval,fwllh,bwllh);
  //  fprintf(stderr,"\t-> [%s] stop\n",__FUNCTION__ );
  //  exit(0);
  if(DOTRANS){
    for(int i=0;i<tk_l;i++)
      for(int j=0;j<tk_l;j++){
	//	fprintf(stderr,"i:%d j:%d\n",i,j);
	trans[i][j] = calc_trans(i,j,P);
      }
    printmatrixf((char*)"transitions.txt",trans,tk_l,tk_l);
  }
  return qval;
  exit(0);
  
}
/*
double fastPSMC::fwllh(){
 double tmp[tk_l];
 for(int i=0;i<tk_l;i++)
   tmp[i] = fw[i][windows.size()];
 pix = addProtectN(tmp,tk_l);
 assert(!isnan(pix));
 fprintf(stderr,"\t-> forward(pic) llh:%f\n",pix);
 return pix;
}


double fastPSMC::bwllh(){
  double tmp[tk_l];
  for(int i=0;i<tk_l;i++)
    tmp[i] = bw[i][1]+stationary[i]+emis[i][1];
  double tmptmp= addProtectN(tmp,tk_l);
  assert(!isnan(tmptmp));
  fprintf(stderr,"\t <- backward llh:%f\n",tmptmp);
  return tmptmp;
}
*/
