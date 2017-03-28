#include <vector>
#include <cassert>
#include <cmath>
#include "psmcreader.h"
#include "main_psmc.h"
#include "hmm_psmc.h"
#include "compute.c"

void printmatrixf(char *fname,double **m,int x,int y);
void printarrayf(char *fname,double *m,int x);

void printarray(FILE *fp,double *ary,int l){
  for(int i=0;i<l;i++)
    fprintf(fp,"%d)\t%f\n",i,ary[i]);

}
void printmatrix(FILE *fp,double **mat,int x,int y){
  fprintf(fp,"#printmatrix with x:%d y:%d\n",x,y);
  for(int i=0;i<y;i++){
    for(int j=0;j<x-1;j++)
      fprintf(fp,"%f\t",mat[j][i]);
    fprintf(fp,"%f\n",mat[x-1][i]);
  }

}

double qkFunction(unsigned i, double pix, unsigned numWind,double **nP,double **PP){
  /*
    for (unsigned l = 1; l < numWind + 1; l++) //This block is needed if eimission probabilities depend on estimated parameters, e.g. on time disctretisation 
    qi += log(emis[K][l])*fw[K][l]*bw[K][l];
    qi /= pix;
  */
  
  double qi[7];
  qi[0] = nP[1][i]+PP[1][i];
  qi[1] = nP[2][i]+PP[2][i];
  qi[2] = nP[3][i]+PP[3][i];
  qi[3] = nP[4][i]+PP[4][i];
  qi[4] = nP[5][i]+PP[5][i];
  qi[5] = nP[6][i]+PP[6][i];
  qi[6] = nP[7][i]+PP[7][i];

  return addProtectN(qi,7);
}


void ComputeGlobalProbabilities(double *tk,int tk_l,double **P,double *epsize,double rho){
  ComputeP1(tk,tk_l,P[1],epsize,rho);
  ComputeP5(tk,tk_l,P[5],epsize);
  ComputeP6(tk,tk_l,P[6],epsize,rho);
  ComputeP2(tk_l,P[2],P[5]);
  ComputeP3(tk,tk_l,P[3],epsize,rho);
  ComputeP4(tk,tk_l,P[4],epsize,rho);
  ComputeP7(tk,tk_l,P[7],epsize,rho);
  ComputeP0(tk_l,P[0],P[5]);
}



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
  if(!from_infile)
    for (int i = 0; i <l; i++)
      ary[l]=1;
  else{
    memcpy(ary,from_infile,l*sizeof(double));
  }
	
}


void fastPSMC::calculate_FW_BW_PP_Probs(){
    //we first set the initial fwprobs to stationary distribution
    for(int i=0;i<tk_l;i++)
      fw[i][0] = stationary[i];

  

    //we now loop over windows.
    //v=0 is above and is the initial distribution, we therefore plug in at v+1
    for(int v=0;v<windows.size();v++){
      ComputeRs(v,fw);//<-prepare R1,R2
#if 0
      printarrayf("r1",R1,tk_l);
      printarrayf("r2",R2,tk_l);
      // exit(0);
#endif
      fw[0][v+1] = addProtect3(fw[0][v]+P[1][0] , R1[0]+P[3][0] , fw[0][v]+P[4][0])+emis[0][v+1] ;
      fprintf(stderr,"fw[0][1]:%f\n",0,v+1);
      for (unsigned i = 1; i < tk_l; i++){
	fprintf(stderr,"l1:%f l2:%f l3:%f l4:%f\n",fw[i][v]+P[1][i] , R2[i-1]+P[2][i-1] , R1[i]+P[3][i] , fw[i][v]+P[4][i]);
	fw[i][v+1]= addProtect4(fw[i][v]+P[1][i] , R2[i-1]+P[2][i-1] , R1[i]+P[3][i] , fw[i][v]+P[4][i])+emis[i][v+1];
	fprintf(stderr,"fw[%d][%d]:%f\n",i,v+1,fw[i][v+1]);
	//	exit(0);
      }
      //      break;
    }
    fprintf(stderr,"PROGRSM BRESKING IN FW BW\n");
    return ;
    double tmp[tk_l];
    for(int i=0;i<tk_l;i++)
      tmp[i] = fw[i][windows.size()];

    pix = addProtectN(tmp,tk_l);
    fprintf(stderr,"forward(pic) llh:%f\n",pix);


    //now do backward algorithm
    for(int i=0;i<tk_l;i++)
      bw[i][windows.size()] = stationary[i];

    //we plug in values at v-1, therefore we break at v==1
    for(int v=windows.size();v>0;v--){
      ComputeRs(v,bw);//<-prepare R1,R2
      bw[0][v-1] = addProtect3(bw[0][v]+P[1][0] , R1[0]+P[3][0] , bw[0][v]+P[4][0])-emis[0][v] ;
      for (unsigned i = 1; i < tk_l; i++)
	bw[i][v-1] = addProtect4(stationary[i]+bw[i][v]+emis[i][v]+P[1][i],R2[i-1]+P[2][i-1],R1[i]+P[3][i],stationary[i]+bw[i][v]+emis[i][v]+P[4][i])-stationary[i];
      
    }
    
    for(int i=0;i<tk_l;i++)
      tmp[i] = bw[i][windows.size()];
    fprintf(stderr,"backward llh:%f\n",addProtectN(tmp,tk_l));


    //calculate post prob per window per state
    for(int v=1;v<windows.size();v++){
      for(int j=0;j<tk_l;j++)
	pp[j][v] = fw[j][v]+bw[j][v];
	
    }
    
    fprintf(stderr,"[%s] stop\n",__FUNCTION__ );
  }



void fastPSMC::allocate(int tk_l_arg){
  int numWindows = windows.size();
  fprintf(stderr,"\t-> [%s]: will allocate tk with length: %d\n",__FUNCTION__,tk_l_arg);
  tk_l = tk_l_arg;
  //tk = new double[tk_l];
 
  //  printarray(stderr,tk,tk_l);
  stationary = new double[tk_l];
  R1 = new double[tk_l];
  R2 = new double[tk_l];
  fw = new double *[tk_l];
  bw = new double *[tk_l];
  pp = new double *[tk_l];
  emis = new double *[tk_l];
  for(int i=0;i<tk_l;i++){
    emis[i] = new double[numWindows+1];
    fw[i] = new double[numWindows+1];
    bw[i] = new double[numWindows+1];
    pp[i] = new double[numWindows+1];

  }
  fprintf(stderr,"\t-> emission allocated with [%d][%d]\n",tk_l,numWindows+1);
  P = new double *[8];
  PP= new double *[8];
  for(int i=0;i<8;i++){
    P[i] = new double[tk_l];
    PP[i]= new double[tk_l];
  }
 
}
/*
  Function will set the indices for the windows
 */
void fastPSMC::setWindows(perpsmc *perc,char *chooseChr,int start,int stop,int block){
  myMap::const_iterator it = iter_init(perc,chooseChr,start,stop);
  int beginIndex =0;
  int endIndex=0;
  int beginPos = 1;
  int endPos = beginPos+block-1;
  
  while(1){
    wins w;
    if(endPos>perc->pos[perc->last-1])
      break;
    
    while(perc->pos[beginIndex]<beginPos)
      beginIndex++;
    while(perc->pos[endIndex]<endPos)
      endIndex++;
#if 0
    fprintf(stdout,"\t-> endpiadsf:%d\n",perc->pos[endIndex]);
    fprintf(stdout,"\t-> winsize:%d bp:%d,ep:%d bi:%d ei:%d ei-bi:%d\n",block,beginPos,endPos,beginIndex,endIndex,endIndex-beginIndex);
#endif
    w.from = beginIndex;
    w.to = endIndex+1;
    windows.push_back(w);
    beginPos+=block;
    endPos+=block;
  }
  fprintf(stderr,"\t-> [%s] chr: \'%s\' has number of windows:%lu\n",__FUNCTION__,chooseChr,windows.size());
}


/* 
  Calculate emission probabilityes
  tk array of length tk_l
  lambda array effective population sizes
  both has length tk_l
  
  emission probablities will be put in the **emis
  
  stationary(i) = exp(-sum_{j=0}^{i-1}{tau_j/lambda_j}*P2[i])
 */
void calculate_emissions(double *tk,int tk_l,double *gls,std::vector<wins> &windows,double mu,double **emis){
  fprintf(stderr,"\t-> [Calculating emissions with tk_l:%d and windows.size():%lu:%s ] start\n",tk_l,windows.size(),__TIME__);
  //initialize the first:
  for(int j=0;j<tk_l;j++)
    emis[j][0] = log(0);
 
  double tmp[windows.size()];
  for(int v=0;v<windows.size();v++){
    for(int j=0;j<tk_l;j++){
      emis[j][v+1] = 0;
      double inner = exp(-2.0*tk[j]*mu);
      for(int i=windows[v].from;i<windows[v].to;i++)
	emis[j][v+1] += log((exp(gls[i*2])/4.0) *inner + (exp(gls[2*i+1])/6)*(1-inner));//<- check
    }
  }
  fprintf(stderr,"\t-> [Calculating emissions with tk_l:%d and windows.size():%lu:%s ] stop\n",tk_l,windows.size(),__TIME__);
}


/*
  Calculate stationary distrubution
  tk array of length tk_l
  lambda array effective population sizes
  both has length tk_l
  
  stationary distribution will be put in results array, also of length tk
  
  stationary(i) = exp(-sum_{j=0}^{i-1}{tau_j/lambda_j}*P2[i])
 */
void fastPSMC::calculate_stationary(double *tk,int tk_l,double *lambda,double *results,double **P){
  results[0] = P[2][0];//fix this
  for(int i=1;i<tk_l;i++){
    results[i]  = P[2][i]+P[0][i-1];
  }
#if 0 //check it sums to one
  double tmp=0;
  for(int i=0;i<tk_l;i++)
    tmp += exp(results[i]);
  fprintf(stderr,"\t-> sum of stationary:%f\n",tmp);
  exit(0);

#endif

}
void fastPSMC::ComputePii(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary){
  static double *workspace = new double [numWind];
  ComputeP11(numWind,tk_l,P[1],PP[1],fw,bw,stationary,workspace);
  ComputeP22(numWind,tk_l,P,PP[2],fw,bw,stationary);
  ComputeP33(numWind,tk_l,P[3],PP[3],fw,bw,stationary);
  ComputeP44(numWind,tk_l,P[4],PP[4],fw,bw,stationary,workspace);
  ComputeP55(numWind,tk_l,P,PP[5],fw,bw,stationary);
  ComputeP66(numWind,tk_l,P,PP[6],fw,bw,stationary);
  ComputeP77(numWind,tk_l,P,PP[7],fw,bw,stationary);
}


void fastPSMC::make_hmm(double *tk,int tk_l,double *gls,double *epsize){
  fprintf(stderr,"\t-> [%s] start\n",__FUNCTION__ );
  //prepare global probs
  ComputeGlobalProbabilities(tk,tk_l,P,epsize,rho);//only the P* ones
  //calculate emissions
  calculate_emissions(tk,tk_l,gls,windows,mu,emis);
  //  printmatrixf("emis",emis,tk_l,windows.size()+1);exit(0);
  
  calculate_stationary(tk,tk_l,epsize,stationary,P);
  //    printarray(stderr,stationary,tk_l);
  calculate_FW_BW_PP_Probs();
  fprintf(stderr,"\t-> [%s] stop\n",__FUNCTION__ );
}
