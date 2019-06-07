#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "compute.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

double lprod(double a,double b);
double lprod(double a,double b,double c);
double lprod(double a,double b,double c,double d);

#if 0 //this is undef 9feb 2019 and subsituted with the function below
void ComputeBR1(int tk_l, double *bR1, double **P, double **bw,double **emis,int w){
  bR1[tk_l - 1] = log(0);
  double tmp_denom = 0.0;
  for (int i = 0; i < tk_l - 1; i++)
    tmp_denom += P[5][i];
  for (int i = tk_l - 2; i >= 0 ; i--){
    //      bR1[i] =  addProtect2(bR1[i+1] , lprod(bw[i+1][l+1],emis[i+1][l+1],stationary[i+1]));
    bR1[i] =  addProtect3(bR1[i+1] , lprod(bw[i+1][w+1],emis[i+1][w+1]), -tmp_denom);
    tmp_denom -= P[5][i];
  }
}
#endif 

void ComputeBR1(int tk_l, double *bR1, double **P, double *stationary, double **bw,double **emis,int w){
  bR1[tk_l - 1] = log(0);
  double tmp_denom = 0.0;
  for (int i = 0; i < tk_l - 1; i++)
    tmp_denom += P[5][i];
  for (int i = tk_l - 2; i >= 0 ; i--){
    bR1[i] =  addProtect2(lprod(bR1[i+1],P[0][i+1]) , lprod(bw[i+1][w+1],emis[i+1][w+1],stationary[i+1]));
    bR1[i] = lprod(bR1[i], -P[0][i]);
  }
}

double addProtect3(double a,double b, double c){
  if(isinf(a)&&isinf(b)&&isinf(c))
    return log(0.0);
  //function does: log(exp(a)+exp(b)+exp(c)) while protecting for underflow
  double maxVal;// = std::max(a,std::max(b,c));
  if(a>b&&a>c)
    maxVal=a;
  else if(b>c)
    maxVal=b;
  else
    maxVal=c;
  double sumVal = exp(a-maxVal)+exp(b-maxVal)+exp(c-maxVal);
  return log(sumVal) + maxVal;
}

double addProtect4(double a,double b, double c,double d){
  //function does: log(exp(a)+exp(b)+exp(c)) while protecting for underflow
  double maxVal;// = std::max(a,std::max(b,c));
  if(a>b&&a>c&&a>d)
    maxVal=a;
  else if(b>a&&b>c&&b>d)
    maxVal=b;
  else if(c>a&&c>b&&c>d)
    maxVal=c;
  else
    maxVal=d;
  double sumVal = exp(a-maxVal)+exp(b-maxVal)+exp(c-maxVal)+exp(d-maxVal);
  return log(sumVal) + maxVal;
}
double addProtectN(double a[],int len){
  //function does: log(sum(exp(a))) while protecting for underflow
  double maxVal = a[0];
  
  for(int i=1;i<len;i++){
    if(maxVal<a[i])
      maxVal=a[i];
  }

  double sumVal = 0;
  for(int i=0;i<len;i++)
    sumVal += exp(a[i]-maxVal);

  return log(sumVal) + maxVal;
}
double addProtect2(double a,double b){
  //function does: log(exp(a)+exp(b)) while protecting for underflow
  if(!isfinite(a)&&!isfinite(b) )
    return a;
  double maxVal;// = std::max(a,b));
  if(a>b)
    maxVal=a;
  else
    maxVal=b;
  double sumVal = exp(a-maxVal)+exp(b-maxVal);
  return log(sumVal) + maxVal;
}


void ComputeP11(unsigned numWin,int tk_l,double *P1,double *PP1,double **fw,double **bw,double *workspace,double **emis){
  for (unsigned i = 0; i < tk_l; i++){
    workspace[0] = log(0);
    for (unsigned l = 1; l < numWin; l++){
      //      fprintf(stderr,"l:%d\n",l);
      //fprintf(stderr,"fw[][]:%f\n",fw[i][l]);
      //fprintf(stderr,"fw[][]:%f\n",bw[i][l]);
      workspace[l] = lprod(fw[i][l],P1[i],bw[i][l+1],emis[i][l+1]);
      //NOTE: In appendix of K.H. paper it seems to be an extra emission probability for site l+1, it is already inside bw[]
    }
    PP1[i] = addProtectN(workspace,numWin);
  }
}

void ComputeP22(unsigned numWind,int tk_l,double **P,double *PP2,double **fw,double **bw,double **emis){
  double R1[tk_l];
  double R2[tk_l];
  double tmp[tk_l];
  for (unsigned i = 0; i < tk_l; i++)
    PP2[i] = log(0);
  for (unsigned l = 1; l < numWind; l++) {
    R1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = addProtect2(R1[i+1],fw[i+1][l]);
    double tmp = log(0);
    for (unsigned i = 0; i < tk_l ; i++){
      R2[i] = addProtect3(lprod(tmp,P[5][i]) , lprod(fw[i][l],P[6][i]) , lprod(R1[i],P[7][i]));
      tmp = R2[i];
    }

    for (unsigned i = 1; i < tk_l; i++)
      PP2[i] =addProtect2(PP2[i] , lprod(R2[i-1],P[2][i],bw[i][l+1],emis[i][l+1]));//CHECK
  }
}

void ComputeP33(unsigned numWind,int tk_l,double *P3,double *PP3,double **fw,double **bw,double **emis){
  double R1[tk_l];
  for (unsigned i = 0; i < tk_l; i++)
    PP3[i] = log(0);
  for (unsigned l = 1; l < numWind; l++){
    R1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--){
      R1[i] = addProtect2( R1[i+1] , fw[i+1][l]);
      //   fprintf(stderr,"R1[%d]:%f\n",i,R1[i]);
    }
    for (unsigned i = 0; i < tk_l - 1; i++){
      //      fprintf(stderr,"%d) PP[3]:%f lprod:%f\n",i,PP3[i],lprod(R1[i],P3[i],bw[i][l+1],emis[i][l+1]));
      PP3[i] = addProtect2(PP3[i],lprod(R1[i],P3[i],bw[i][l+1],emis[i][l+1]));
      //      fprintf(stderr,"%d) PP[3]:%f lprod:%f\n",i,PP3[i],lprod(R1[i],P3[i],bw[i][l+1],emis[i][l+1]));
    }
  }
  for(int i=0;1&&i<tk_l;i++){//CHECK THIS LASTER 12oct 2017
    assert(!isnan(PP3[i]));
  }

}


void ComputeP44(unsigned numWind,int tk_l,double *P4,double *PP4,double **fw,double **bw,double *workspace,double **emis){

  for (unsigned i = 0; i < tk_l; i++){
    workspace[0] = log(0);
    //    int ntot=0;
    for (unsigned l = 1; l < numWind; l++){
      workspace[l] = lprod(fw[i][l],P4[i],bw[i][l+1],emis[i][l+1]);
      //fprintf(stderr,"fw:%f P4:%f bw:%f emis:%f\n",fw[i][l],P4[i],bw[i][l+1],emis[i][l+1] );
      //ntot++;
    }
    //    fprintf(stderr,"nont:%d\n",ntot);
    PP4[i] = addProtectN(workspace,numWind);
  }

}
void ComputeP55(unsigned numWind,int tk_l,double **P,double *PP5,double **fw,double **bw,double *stationary,double **emis){
  double R1[tk_l];
  double R2[tk_l];
  double bR1[tk_l];

  for (unsigned i = 0; i < tk_l; i++)
    PP5[i] = log(0);
  for (unsigned l = 1; l < numWind; l++){
    R1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = addProtect2(R1[i+1] , fw[i+1][l]);
    //    fprintf(stderr,"ComputeP55_R1[%d]:\t%f\t%f\t%f\n",l,R1[0],R1[1],R1[2]);	
    double tmp = log(0);
    for (unsigned i = 0; i < tk_l ; i++){
      R2[i] = addProtect3(lprod(tmp,P[5][i]),lprod(fw[i][l],P[6][i]),lprod(R1[i],P[7][i]));
      tmp = R2[i];
    }
    //fprintf(stderr,"ComputeP55_R2[%d]:\t%f\t%f\t%f\n",l,R2[0],R2[1],R2[2]);
    ComputeBR1(tk_l,bR1,P,stationary,bw,emis,l);
    for (unsigned i = 1; i < tk_l - 1; i++)
      PP5[i] = addProtect2(PP5[i],lprod(R2[i-1],P[5][i],bR1[i]));//<- CHECK ptgi
  }
}

void ComputeP66(unsigned numWind,int tk_l,double **P,double *PP6,double **fw,double **bw,double *stationary,double **emis){
  double bR1[tk_l];
  for (unsigned i = 0; i < tk_l; i++)
    PP6[i] = log(0);

  for (unsigned l = 1; l < numWind; l++){
    /*
    bR1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      bR1[i] = addProtect2(bR1[i+1] , lprod(bw[i+1][l+1],emis[i+1][l+1],stationary[i+1]));
    */
    ComputeBR1(tk_l,bR1,P,stationary,bw,emis,l);
    for (unsigned i = 0; i < tk_l - 1; i++)
      PP6[i] =addProtect2(PP6[i], lprod(fw[i][l],P[6][i],bR1[i]));//<- CHECK btgi
  }
}

void ComputeP77(unsigned numWind,int tk_l,double **P,double *PP7,double **fw,double **bw,double *stationary,double **emis){
  double R1[tk_l];
  double bR1[tk_l];
  for (unsigned i = 0; i < tk_l; i++)
    PP7[i] = log(0);
  for (unsigned l = 1; l < numWind; l++){
    R1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = addProtect2(R1[i+1] , fw[i+1][l]);
    /*
    bR1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      bR1[i] = addProtect2(bR1[i+1] , lprod(bw[i+1][l+1],emis[i+1][l+1],stationary[i+1]));
    */
    ComputeBR1(tk_l,bR1,P,stationary,bw,emis,l);
    for (unsigned i = 0; i < tk_l - 1; i++)
      PP7[i] = addProtect2(PP7[i],lprod(R1[i],P[7][i],bR1[i]));//<-CHECK ptgi
  }
}

//TODO: tk_l is in fact the number of time intervals
/*
  prob of not recombining given states
*/

void ComputeP1(double *tk,int tk_l,double *P,const double *epsize,double rho){
  //  fprintf(stderr,"[%s] tks=(%f,%f) epssize=(%f,%f)\n",__FUNCTION__,tk[0],tk[1],epsize[0],epsize[1]);
  for (unsigned i = 0; i < tk_l-1; i++){
    P[i] = 1.0/(1.0+epsize[i]*2.0*rho);
    P[i] *= exp( -rho*2.0*tk[i] ) - exp(-rho*2.0*tk[i+1]-(tk[i+1]-tk[i])/epsize[i]);
    P[i] /= 1.0 - exp( -(tk[i+1]-tk[i])/epsize[i] );
    P[i] = log(P[i]);
    //    fprintf(stderr,"P1[%d]:%f\n",i,P[i]);
  }

  //Last interval ends with +infinity
  P[tk_l-1] = log(1.0/(1.0+epsize[tk_l-1]*2.0*rho)* exp( -rho*2.0*tk[tk_l-1] ));
  //  fprintf(stderr,"P1[%d]:%f\n",tk_l-1,P[tk_l-1]);
}
/*
  prob of coalscne happens after timepoint k given that it happens after k-1

 */

void ComputeP5(double *tk,int tk_l,double *P,const double *epsize){
  for (unsigned i = 0; i < tk_l-1; i++)
    P[i] = log(exp( -(tk[i+1] - tk[i])/epsize[i] ));
  P[tk_l-1] = log(0.0);
}


void ComputeP6(double *tk,int tk_l,double *P,const double *epsize,double rho){
    for (unsigned i = 0; i < tk_l-1; i++){
      P[i] = 1/(1-exp(-(tk[i+1]-tk[i])/epsize[i]));
      P[i] *= exp(-(tk[i+1]-tk[i])/epsize[i]);
      double tmp = exp(-2*rho*tk[i]);
      tmp -= 1/(1-2*rho*epsize[i])*exp(-2*rho*tk[i+1]);
      tmp += 2*rho*epsize[i]/(1 - 2*rho*epsize[i])*exp(-2*rho*tk[i]-(tk[i+1]-tk[i])/epsize[i]);
      P[i] *= tmp;
      P[i] = log(P[i]);
    }
    P[tk_l - 1] = log(0.0);
}

/*
void ComputeP6(double *tk,int tk_l,double *P,const double *epsize,double rho){
  double en=2.3;
  double to=3.4;
  double tre=addProtect2(en,to);
  fprintf(stderr,"%f %f %f\n",en,to,tre);
    for (unsigned i = 0; i < tk_l-1; i++){
      double inner = -(tk[i+1]-tk[i])/epsize[i];
      fprintf(stderr,"inner:%f\n",inner);
      P[i] = log(1)-log(1-exp(inner));
      fprintf(stderr,"pi no addprotect: %f\n",P[i]);
      P[i] = log(1)-addProtect2(log(1.0),-log(-exp(-inner)));
      fprintf(stderr,"pi with addprotect: %f\n",P[i]);
      exit(0);
      P[i] += inner;
      double tmp = exp(-2*rho*tk[i]);
      tmp -= 1/(1-2*rho*epsize[i])*exp(-2*rho*tk[i+1]);
      tmp += 2*rho*epsize[i]/(1 - 2*rho*epsize[i])*exp(-2*rho*tk[i]+inner);
      P[i] += log(tmp);
      fprintf(stderr,"pi:%f\n",P[i]);
      exit(0);
    }
    P[tk_l - 1] = log(0.0);
}
*/


void ComputeP2(int tk_l,double *P2,double *P5){
  for (unsigned i = 0; i < tk_l; i++){

    P2[i] = log(1.0 - exp(P5[i]));
    //    fprintf(stderr,"%d):p5:%f value:%f\n",i,P5[i],P2[i]);
  //  P2[tk_l-1]=-0.0;
  }
}



void ComputeP3(double *tk,int tk_l,double *P3,const double *epsize,double rho){
  for (unsigned i = 0; i < tk_l - 1; i++){
    P3[i] = exp(-tk[i]*2.0*rho);
    P3[i] += epsize[i]*2.0*rho/(1.0 - epsize[i]*2.0*rho)*exp(-(tk[i+1]-tk[i])/epsize[i]-tk[i]*2.0*rho);
    P3[i] -= 1.0/(1.0 - epsize[i]*2.0*rho)*exp(-tk[i+1]*2.0*rho);
    P3[i] = log(P3[i]);
  }
  P3[tk_l-1] = log(0.0);
}

/*
  Check P4

 */

void ComputeP4(double *tk,int tk_l,double *P4,const double *epsize,double rho){
  //  fprintf(stderr,"\t->[%s] rho: %f\n",__FUNCTION__,rho);
  //  fprintf(stderr,"\t->[%s] epsizes=(%f,%f) \n",__FUNCTION__,epsize[0],epsize[1]);
  for (unsigned i = 0; i < tk_l-1; i++){

    //    double fact1 = (exp(-(2.0*rho*tk[i])))/(1.0 - exp(-(tk[i+1]-tk[i])/epsize[i]) );
    double fact1 =(exp(-(2.0*rho*tk[i])))/((1.0 - exp(-(tk[i+1]-tk[i])/epsize[i]) ));

    double part1 = 2.0/(1-2.0*epsize[i]*rho)/(1+2.0*epsize[i]*rho);
    //    fprintf(stderr,"part1:%f epszie:%f\n",part1,epsize[i]);
    //    fprintf(stderr,"\tparst1fact1:%f\n",part1);
    double part1exp1 = -(tk[i+1]-tk[i])/epsize[i];
    double part1exp2 = -(tk[i+1]-tk[i])*2.0*rho;
    part1 *= exp(part1exp1+part1exp2);

    double part2 = (2.0*rho*epsize[i])/(1+2*epsize[i]*rho);

    double part3 = (2*epsize[i]*rho)/(1-2*epsize[i]*rho);
    double part3exp = -2*(tk[i+1]-tk[i])/epsize[i];
    part3 *= exp(part3exp);
    

    double part4exp = -(tk[i+1]-tk[i])/epsize[i];
    double part4 = 2*exp(part4exp);
    
    double fact2= part1+part2-part3-part4;
    if(fact2<0){
      fprintf(stderr,"\t-> This fix shouldnt happen that often: compute.cpp computep4 epsize[%u]:%f\n",i,epsize[i]);
      fact2=fabs(fact2);
    }
    P4[i] = log(fact1)+log(fact2);//log(fact1)+log(fact2);
    
#if 0
    fprintf(stderr,"\t-> fact1: %e fact2: %e fact1*fact2: %f\n",fact1,fact2,fact1*fact2);
    fprintf(stderr,"\t-> part1: %f part2: %f part3: %f part4: %f\n",part1,part2,part3,part4);
    fprintf(stderr,"\t-> part1+part2: %f -part3-part4: %f \n",part1+part2,-part3-part4);
    fprintf(stderr,"\t-> fact2: part1+part2-part3-part4: %f \n",part1+part2-part3-part4);
    fprintf(stderr,"\t-> P4[%d]: %f \n",i,P4[i]);
#endif
    if(isnan(P4[i])){
      fprintf(stderr,"P[4][%d]: %f never happens cmputep4\n",i,P4[i]);
      fprintf(stderr,"\t-> fact1: %e fact2: %e fact1*fact2: %f\n",fact1,fact2,fact1*fact2);
      fprintf(stderr,"\t-> part1: %f part2: %f part3: %f part4: %f\n",part1,part2,part3,part4);
      fprintf(stderr,"\t-> part1+part2: %f -part3-part4: %f \n",part1+part2,-part3-part4);
      fprintf(stderr,"\t-> fact2: part1+part2-part3-part4: %f \n",part1+part2-part3-part4);
      fprintf(stderr,"\t-> P4[%d]: %f \n",i,P4[i]);

      for(int i=0;i<tk_l;i++)
	fprintf(stderr,"epsize[%d]: (%f,%f)\n",i,tk[i],epsize[i]);
      exit(0);
      //exit(0);
      //    assert(P4[i]>=0&&P4[i]<=1);
    }
   
  }

  double top = 2.0*rho*epsize[tk_l-1];
  double bottom = (1.0 + 2.0*rho*epsize[tk_l-1]);
  double expot = exp(-2.0*rho*tk[tk_l-1]);
  //  fprintf(stderr,"TOPTOP top:%f bot:%f exp:%f rho:%f epSize[tk_l-1]:%e\n",top,bottom,expot,rho,epsize[tk_l-1]);
  P4[tk_l-1] = log(2.0*rho*epsize[tk_l-1]/(1.0 + 2.0*rho*epsize[tk_l-1])*exp(-2.0*rho*tk[tk_l-1]));
  // fprintf(stderr,"P4 last:%f\n",exp(P4[tk_l-1]));exit(0);
  assert(exp(P4[tk_l-1])>=0&&exp(P4[tk_l-1])<=1);

}

  
void ComputeP7(double *tk,int tk_l,double *P7,double *P3,const double *epsize,double rho){
  for (unsigned i = 0; i < tk_l - 1; i++)
    P7[i] = log(exp(-2.0*rho*tk[i])-exp(-2.0*rho*tk[i+1]) -exp(P3[i]));
  
  unsigned i = tk_l - 1;
  
  P7[i] = log(0.0);
  //  fprintf(stderr,"P7 (%f,%f)\n",exp(P7[0]),exp(P7[1]));
}
  
void ComputeP0(int tk_l,double *P0,double *P5){ //probability P(T > i)
  P0[0] = P5[0];
  for (unsigned i = 1; i < tk_l; i++)
    P0[i] = P0[i-1]+P5[i];
}

double ComputeXXX(int i,int j,double **P){
  assert(i<j);//check if needed
  double sum =0;
  for(int l=i+1;l<=j-1;l++)
    sum += P[5][l];
  double returnVal = P[7][i]+P[2][j]+sum;

  return returnVal;
}

void ComputeExpectedCoalTime(double *tk,int tk_l,double *expectCoalT,const double *epsize){
  for (unsigned i = 0; i < tk_l-1; i++){
    double expterm = exp(-1/epsize[i]*(tk[i+1] - tk[i]));
    double fact1 = epsize[i]/(1 - expterm);
    double fact2 = 1 + tk[i]/epsize[i] - (1 + tk[i+1]/epsize[i])*expterm;
    expectCoalT[i] = fact1*fact2;
    
#if 0
    fprintf(stderr,"\t-> P4[%d]: %f \n",i,P4[i]);
#endif
	assert(expectCoalT[i] >= tk[i] && expectCoalT[i] <= tk[i+1]);
    if(isnan(expectCoalT[i])){
      fprintf(stderr,"expectCoalT[%d]: %f, aborted\n", i, expectCoalT[i]);
      exit(0);
    }
  }

  //  fprintf(stderr,"TOPTOP top:%f bot:%f exp:%f rho:%f epSize[tk_l-1]:%e\n",top,bottom,expot,rho,epsize[tk_l-1]);
    expectCoalT[tk_l-1] = epsize[tk_l - 1] + tk[tk_l - 1];
    assert(expectCoalT[tk_l-1] >= tk[tk_l - 1]);

}


double calc_trans(int k, int j,double **P){
  double ret;
  if(k<j){
    double sumL=0;
    for(int l=k+1;l<=j-1;l++)
      sumL += P[5][l];

    ret = exp(P[6][k]+P[2][j]+sumL);
    double sum =0;
    for(int i=0;i<k;i++)
      sum += exp(ComputeXXX(i,j,P));//underflow stuff
    ret += sum;//underflow stuff
    ret = log(ret);
  }else if(j==k){
    double sum =0;
    for(int i=0;i<k;i++){
      sum += exp(ComputeXXX(i,k,P));//underflow stuff
    }

    ret = log(exp(P[1][k])+exp(P[4][k])+sum);
  }else if(k>j){
    double sum =0;
    for(int i=0;i<j;i++){
      sum += exp(ComputeXXX(i,j,P));//underflow stuff
    }
    ret = log(exp(P[3][j])+sum);//addProtect2(P[3][j],log(sum));
  }else{
    assert(0==1);
    ret=0;//<- is never set, just to silence compiler
  }

  assert(!isnan(ret));

  return ret;
}

#define NUM_LIN 1

void ComputeU1(double *tk,int tk_l,double **U,const double *epsize,double rho){
  for (unsigned i = 1; i < tk_l-1; i++){
    double exponent = (tk[i+1]-tk[i])/epsize[i];
    double fact1 = 1/(1-exp(-exponent));
    double fact2 = -(NUM_LIN+1)*exponent;
    fact2 = exp(fact2)-(NUM_LIN+1)*exp(-exponent);
    fact2 = 1+fact2/NUM_LIN;
    U[1][i] = fact1*fact2;
    //		P[i] = log(P[i]);
  }
  U[1][0] = 0;
  U[1][tk_l-1] = 1;
}

//wrong
void ComputeU2(int tk_l,double **U){
  U[2][0] = U[4][0];
  for (unsigned i = 1; i < tk_l; i++){
    U[2][i] = U[2][i-1]*U[5][i]+U[4][i];
  }
}

//Check that ComputeU3 = ComputeP4 for #define NUM_LIN 1
void ComputeU3(double *tk,int tk_l,double **U,const double *epsize,double rho){
  for (unsigned i = 0; i < tk_l-1; i++){
    double exponent = (tk[i+1]-tk[i])/epsize[i];
    double term1 = -exp(-2*rho*tk[i] - exponent)*(1 + NUM_LIN)/NUM_LIN;
    double term2 = 2*rho/(1/epsize[i] + 2*rho)*exp(-2*rho*tk[i]);
    double term3 = exp(-2*rho*tk[i] - (1 + NUM_LIN)*exponent);
    term3 *= 2*rho/NUM_LIN/(-NUM_LIN/epsize[i] + 2*rho);
    double term4 = exp(-2*rho*tk[i+1] - exponent);
    term4 *= (1/(NUM_LIN/epsize[i] - 2*rho) + 1/(1/epsize[i] + 2*rho))/epsize[i];
    U[3][i] = term1 + term2 + term3 + term4;
    U[3][i] /= 1 - exp(-exponent);
    //		P[i] = log(P[i]);
  }
  U[3][tk_l-1] = 2*rho/(1/epsize[tk_l-1] + 2*rho)*exp(-2*rho*tk[tk_l-1]);
  //  fprintf(stderr,"U[3]: %f\n",U[3][tk_l-1]);exit(0);
}

//Check that ComputeU4 = ComputeP7 for #define NUM_LIN 1
void ComputeU4(double *tk,int tk_l,double **U,const double *epsize,double rho){
  for (unsigned i = 0; i < tk_l-1; i++){
    double exponent = (tk[i+1]-tk[i])/epsize[i];
    U[4][i] = -exp(-2*rho*tk[i] - NUM_LIN*exponent) + exp(-2*rho*tk[i+1]);
    U[4][i] = U[4][i]*2*rho/(NUM_LIN/epsize[i] - 2*rho);
    //		P[i] = log(P[i]);
  }
  U[4][tk_l-1]=0;
}

//Check that ComputeU5 = ComputeP5 for #define NUM_LIN 1
void ComputeU5(double *tk,int tk_l,double **U,const double *epsize,double rho){
  for (unsigned i = 0; i < tk_l-1; i++){
    U[5][i] = exp( - NUM_LIN*(tk[i+1] - tk[i])/epsize[i] );
    //		P[i] = log(P[i]);
  }
  U[5][tk_l-1] = 0.0;
  //	P[tk_l-1] = log(0.0);
}

//Check that ComputeU6 = ComputeP6 for #define NUM_LIN 1
void ComputeU6(double *tk,int tk_l,double **U,const double *epsize,double rho){
  for (unsigned i = 0; i < tk_l-1; i++){
    double exponent = (tk[i+1]-tk[i])/epsize[i];
    double fact1 = 1/(1-exp(-exponent));
    double term1 = -exp(-2*rho*tk[i])*(-1+exp(-exponent));
    double term2 = -exp(-2*rho*tk[i]) + exp(-exponent - 2*rho*tk[i+1]);
    term2 = term2/epsize[i]/(1/epsize[i] + 2*rho);
    U[6][i] = fact1*(term1 + term2)-U[3][i];
  }
  U[6][tk_l-1] = 0;
}

//wrong
void ComputeU7(double *tk,int tk_l,double **U,const double *epsize,double rho){
  for (unsigned i = 0; i < tk_l-1; i++){
    U[7][i] = 1 - U[5][i];
  }
}

//wrong
void ComputeU8(double *tk,int tk_l,double **U,const double *epsize,double rho){
  U[8][0] = U[6][0];
  for (unsigned i = 1; i < tk_l; i++){
    U[8][i] = U[8][i-1]*U[7][i]+U[4][i];
  }
}

//Check that ComputeU9 = ComputeP3 for #define NUM_LIN 1
void ComputeU9(double *tk,int tk_l,double **U,const double *epsize,double rho){
  for (unsigned i = 0; i < tk_l-1; i++){
    U[9][i] = exp(-2*rho*tk[i]) - exp(-2*rho*tk[i+1]) - U[4][i];
  }
  U[9][tk_l-1] = 0.0;
}

void ComputeU10(double *tk,int tk_l,double **U,const double *epsize,double rho){
  U[10][0] = U[6][0];
  for (unsigned i = 1; i < tk_l; i++){
    U[10][i] = U[10][i-1]*U[5][i]+U[6][i];
  }
}

void ComputeR3_original(int tk_l,double **fw,double **P,double **U,double *R3,int v){
  R3[0] = fw[0][v]*U[10][0];
  fprintf(stderr,"ComputeR3: fw[0][%d] U[10][0]:%f R3[0]:%f\n",v,fw[0][v],U[10][0],R3[0]);
  for (unsigned i = 1; i < tk_l; i++){
    R3[i] = R3[i-1]*exp(P[5][i])+fw[i][v]*U[10][i];
    fprintf(stderr,"ComputeR3: R3[%d]:%f R3[%d]:%f P[2][%d]:%f fw[%d][%d]:%f U[10][%d]:%f\n",i,R3[i],i-1,R3[i-1],i,P[2][i],i,v,fw[i][v],i,U[10][i]);
  }
}

void ComputeR3(int tk_l,double **fw,double **P,double **U,double *R3,int v){
  for (int k = 0; k < tk_l; k++){

    double r3 = 0;
    for (int i = 0; i < k; i++){
      double part1=exp(P[6][i]), part2=1.0;
      for (int l = 0; l < i-1; l++) {
	double prod = 1.0;
	for (int j = l+1; j < i+1; j++){
	  prod *= exp(P[5][j]);
	}
	prod *= exp(P[7][l]);
	part1 += prod;
      }
      for (int j = i+1; j < k+1; j++){
	part2 *= exp(P[5][j]);
      }
      // fprintf(stderr,"r3:%f fw:%f part1:%f part2:%f fw*part1*part2:%f\n",r3,fw[i][v],part1,part2,fw[i][v]*part1*part2);
      r3 += fw[i][v]*part1*part2;
    }
    //    fprintf(stderr,"r3[%d] = %f k:%d\n",k, r3,k);   
  }

  
  R3[0] = fw[0][v]*U[10][0];
  //fprintf(stderr,"ComputeR3: fw[0][%d]:%f U[10][0]:%f R3[0]:%f\n",v,fw[0][v],U[10][0],R3[0]);
  for (unsigned i = 1; i < tk_l; i++){
    R3[i] = R3[i-1]*exp(P[5][i])+fw[i][v]*U[10][i];
    //fprintf(stderr,"ComputeR3: R3[%d]:%f R3[%d]:%f P[2][%d]:%f fw[%d][%d]:%f U[10][%d]:%f\n",i,R3[i],i-1,R3[i-1],i,P[2][i],i,v,fw[i][v],i,U[10][i]);
  }
  //  exit(0);
}


//original
void NextFW_original(int tk_l,double **P,double **U,double **fw,int v,double **emis,double *R1, double *R3){
  ComputeR3(tk_l,fw,P,U,R3,v);
  fw[0][v+1] = fw[0][v]*(P[1][0]+U[1][0]+U[3][0])+R3[0]*P[2][0]+R1[0]*U[9][0];
  fw[0][v+1] = fw[0][v+1]*exp(emis[0][v]);
  for (unsigned i = 1; i < tk_l; i++){
    fw[i][v+1] = fw[i][v]*(P[1][i]+U[1][i]+U[3][i])+R3[i]*P[2][i]+exp(R1[i])*(U[8][i-1]*U[7][i]+U[9][i]);
    fw[i][v+1] = fw[i][v+1]*exp(emis[i][v]);
  }
}
//p in log
//fw in normal
//u in normal
//emis in log
//r1 and r3 in normal
void NextFW(int tk_l,double **P,double **U,double **fw,int v,double **emis,double *R1, double *R3){
  ComputeR3(tk_l,fw,P,U,R3,v);


for (int k = 0; k < tk_l; k++){
	double r3 = 0.0;
	for (int i = 0; i < k+1; i++){
		double part1=exp(P[6][i]);
		for (int l = 0; l <= i-1; l++){
			double prod = 1.0;
			for (int j = l+1; j < i+1; j++)
				prod *= exp(P[5][j]);
			prod *= exp(P[7][l]);
			part1 += prod;
		}
		double part2=1.0;
		for (int j = i+1; j < k+1; j++)
			part2 *= exp(P[5][j]);
		r3 += fw[i][v]*part1*part2;
	}
	fprintf(stderr, "r3[%d] = %f\n", k, r3);
}

 fprintf(stderr, "special case r3[0] = %f\n", fw[0][v]*exp(P[6][0]));
 fprintf(stderr, "special case r3[0] = %f\n", fw[0][v]*U[10][0]);
 double tmp = fw[1][v]*(exp(P[6][1])+exp(P[7][0]+P[5][1]))+fw[0][v]*exp(P[6][0]+P[5][1]);
 fprintf(stderr, "special case r3[0] = %f\n", tmp);

  //  fw[0][v+1] =   fw[0][v]*(exp(P[1][0])+U[3][0])+R3[0]*exp(P[2][0])+R1[0]*U[9][0];
  fw[0][v+1] =   fw[0][v]*(exp(P[1][0])+U[3][0])+R1[0]*U[9][0];
  fw[0][v+1] = fw[0][v+1]*exp(emis[0][v+1]);

fprintf(stderr,"R3 test[0]: %f \n", R3[0]);
  //  fprintf(stderr,"U12 test[0]: %f\n",U[1][0]*U[2][0]);
  //  fprintf(stderr  ,"fw[0][%d]:%f lastfw[0][%d]:%f p10:%f U10:%f u30:%f r30:%f p20:%f R1:%f u90:%f emis[0][%d]:%f\n",v+1,log(fw[0][v+1]),v,fw[0][v],exp(P[1][0]),U[1][0],U[3][0],R3[0],exp(P[2][0]),R1[0],U[9][0],v+1,emis[0][v+1]);
  for (unsigned i = 1; i < tk_l; i++){
    fprintf(stderr,"R3 test[%d]: %f \n",i , R3[i]);
    fw[i][v+1] = fw[i][v]*(exp(P[1][i])+U[1][i]*U[2][i-1]+U[3][i])+R3[i-1]*exp(P[2][i])+R1[i]*(U[8][i-1]*U[7][i]+U[9][i]);

    fw[i][v+1] = fw[i][v+1]*exp(emis[i][v+1]);
    //    fprintf(stderr  ,"fw[%d][%d]:%f lastfw[%d][%d]:%f p1i:%f U1i:%f u3i:%f r3i:%f p2i:%f R1:%f u9i:%f emis[0][%d]:%f U8[i-1]:%f U7i:%f\n",i,v+1,log(fw[i][v+1]),i,v,fw[i][v],exp(P[1][i]),U[1][i],U[3][i],R3[i],exp(P[2][i]),R1[i],U[9][i],v+1,emis[i][v+1],U[8][i-1],U[7][i]);
  }
}
