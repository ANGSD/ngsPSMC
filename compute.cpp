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

  for (unsigned l = 1; l < numWind; l++){
    R1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = addProtect2(R1[i+1] , fw[i+1][l]);
  }
  for (unsigned i = 0; i < tk_l; i++)
    PP5[i] = log(0);
  for (unsigned l = 1; l < numWind; l++){
    double tmp = log(0);
    for (unsigned i = 0; i < tk_l ; i++){
      R2[i] = addProtect3(lprod(tmp,P[2][i]),lprod(fw[i][l],P[6][i]),lprod(R1[i],P[7][i]));
      tmp = R2[i];
    }
    bR1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      bR1[i] =  addProtect2(bR1[i+1] , lprod(bw[i+1][l+1],emis[i+1][l+1],stationary[i+1]));
    
    for (unsigned i = 0; i < tk_l - 1; i++)
      PP5[i] = addProtect2(PP5[i],lprod(R2[i],P[5][i],bR1[i],-P[0][i]));//<- CHECK ptgi
  }
}

void ComputeP66(unsigned numWind,int tk_l,double **P,double *PP6,double **fw,double **bw,double *stationary,double **emis){
  double bR1[tk_l];
  for (unsigned i = 0; i < tk_l; i++)
    PP6[i] = log(0);
  for (unsigned l = 1; l < numWind; l++){
    bR1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      bR1[i] = addProtect2(bR1[i+1] , lprod(bw[i+1][l+1],emis[i+1][l+1],stationary[i+1]));
    for (unsigned i = 0; i < tk_l - 1; i++)
      PP6[i] =addProtect2(PP6[i], lprod(fw[i][l],P[6][i],bR1[i],-P[0][i]));//<- CHECK btgi
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
    
    bR1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      bR1[i] = addProtect2(bR1[i+1] , lprod(bw[i+1][l+1],emis[i+1][l+1],stationary[i+1]));

    for (unsigned i = 0; i < tk_l - 1; i++)
      PP7[i] = addProtect2(PP7[i],lprod(R1[i],P[7][i],bR1[i],-P[0][i]));//<-CHECK ptgi
  }
  //DRAGON missing PP7[tk_l-1]
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
    double fact1 =(exp(-(2.0*rho*tk[i])))/(addProtect2(1.0, - exp(-(tk[i+1]-tk[i])/epsize[i]) ));

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
