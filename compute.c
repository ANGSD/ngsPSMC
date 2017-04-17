#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

double addProtect3(double a,double b, double c){
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

  for(int i=1;i<10;i++)
    if(maxVal<a[i])
      maxVal=a[i];

  double sumVal = 0;
  for(int i=1;i<10;i++)
    sumVal += exp(a[i]-maxVal);

  return log(sumVal) + maxVal;
}
double addProtect2(double a,double b){
  //function does: log(exp(a)+exp(b)) while protecting for underflow
  double maxVal;// = std::max(a,b));
  if(a>b)
    maxVal=a;
  else
    maxVal=b;
  double sumVal = exp(a-maxVal)+exp(b-maxVal);
  return log(sumVal) + maxVal;
}


void ComputeP11(unsigned numWin,int tk_l,double *P1,double *PP1,double **fw,double **bw,double *stationary,double *workspace){
  
  for (unsigned i = 0; i < tk_l; i++){
    workspace[i] = log(0);
    for (unsigned l = 1; l < numWin; l++)
      workspace[i] = fw[i][l]+P1[i]+bw[l+1][i]-stationary[i];//NOTE: In appendix of K.H. paper it seems to be an extra emission probability for site l+1, it is already inside bw[]
    PP1[i] = addProtectN(workspace,numWin);
  }
}

void ComputeP22(unsigned numWind,int tk_l,double **P,double *PP2,double **fw,double **bw,double *stationary){
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
      R2[i] = addProtect3(tmp+P[2][i],fw[l][i]+P[6][i],R1[i]+P[7][i]);
      tmp = R2[i];
    }
    for (unsigned i = 1; i < tk_l; i++)
      PP2[i] = R2[i-1]+P[2][i]+bw[i][l+1]-stationary[i];//CHECK
  }
}

void ComputeP33(unsigned numWind,int tk_l,double *P3,double *PP3,double **fw,double **bw,double *stationary){
  double R1[tk_l];
  for (unsigned i = 0; i < tk_l; i++)
    PP3[i] = log(0);
  for (unsigned l = 1; l < numWind; l++){
    R1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = addProtect2( R1[i+1] , fw[i+1][l]);
    for (unsigned i = 0; i < tk_l - 1; i++)
      PP3[i] = addProtect2(PP3[i],R1[i]+P3[i]+bw[i][l]-stationary[i]);
  }
}

void ComputeP44(unsigned numWind,int tk_l,double *P4,double *PP4,double **fw,double **bw,double *stationary,double *workspace){
  for (unsigned i = 0; i < tk_l; i++){
    workspace[i] = log(0);
    for (unsigned l = 1; l < numWind; l++)
      workspace[i] = fw[i][l]+P4[i]+bw[i][l+1]-stationary[i];
    PP4[i] = addProtectN(workspace,numWind);
  }
}

void ComputeP55(unsigned numWind,int tk_l,double **P,double *PP5,double **fw,double **bw,double *stationary){
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
      R2[i] = addProtect3(tmp+P[2][i],fw[l][i]+P[6][i],R1[i]+P[7][i]);
      tmp = R2[i];
    }
    bR1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      bR1[i] = addProtect2(bR1[i+1] , bw[i+1][l+1]);
    for (unsigned i = 0; i < tk_l - 1; i++)
      PP5[i] = addProtect2(PP5[i],R2[i]+P[5][i]+bR1[i]-P[0][i]);//<- CHECK ptgi
    //dragon missing tk_l-1 entry
  }
}

void ComputeP66(unsigned numWind,int tk_l,double **P,double *PP6,double **fw,double **bw,double *stationary){
  double bR1[tk_l];
  for (unsigned i = 0; i < tk_l; i++)
    PP6[i] = log(0);
  for (unsigned l = 1; l < numWind; l++){
    bR1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--)
      bR1[i] = addProtect2(bR1[i+1] , bw[i+1][l+1]);
    for (unsigned i = 0; i < tk_l - 1; i++)
      PP6[i] =addProtect2(PP6[i], fw[i][l]+P[6][i]+bR1[i]-P[0][i]);//<- CHECK btgi
    //dragon missing tk_l-1 entry
  }
}

void ComputeP77(unsigned numWind,int tk_l,double **P,double *PP7,double **fw,double **bw,double *stationary){
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
      bR1[i] = addProtect2(bR1[i+1] , bw[i+1][l+1]);

    for (unsigned i = 0; i < tk_l - 1; i++)
      PP7[i] = addProtect2(PP7[i],R1[i]+P[7][i]+bR1[i]-P[0][i]);//<-CHECK ptgi
  }
  //DRAGON missing PP7[tk_l-1]
}

//TODO: tk_l is in fact the number of time intervals
/*
  prob of not recombining given states
*/

void ComputeP1(double *tk,int tk_l,double *P,double *epsize,double rho){ 
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

void ComputeP5(double *tk,int tk_l,double *P,double *epsize){
  for (unsigned i = 0; i < tk_l-1; i++)
    P[i] = log(exp( -(tk[i+1] - tk[i])/epsize[i] ));
  P[tk_l-1] = log(0.0);
}


void ComputeP6(double *tk,int tk_l,double *P,double *epsize,double rho){
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
  for (unsigned i = 0; i < tk_l-1; i++)
    P2[i] = log(1.0 - exp(P5[i]));
  P2[tk_l-1]=-0.0;
}



void ComputeP3(double *tk,int tk_l,double *P3,double *epsize,double rho){
  for (unsigned i = 0; i < tk_l - 1; i++){
    P3[i] = exp(-tk[i]*2.0*rho);
    P3[i] += epsize[i]*2.0*rho/(1.0 - epsize[i]*2.0*rho)*exp(-(tk[i+1]-tk[i])/epsize[i]-tk[i]*2.0*rho);
    P3[i] -= 1.0/(1.0 - epsize[i]*2.0*rho)*exp(-tk[i+1]*2.0*rho);
    P3[i] = log(P3[i]);
  }
  P3[tk_l-1] = exp(-tk[tk_l-1]*2.0*rho);
}

/*
  Check P4

 */

void ComputeP4_old(double *tk,int tk_l,double *P4,double *epsize,double rho){
  //  fprintf(stderr,"\t->[%s] rho: %f\n",__FUNCTION__,rho);
  for (unsigned i = 0; i < tk_l-1; i++){
    //fprintf(stderr,"\t->[%s] epsize[%d]: %f tk[%d+1]: %f tk[%d]: %f tk[i+1]-tk[i]: %f\n",__FUNCTION__,i,epsize[i],i,tk[i+1],i,tk[i],tk[i+1]-tk[i]);
    double fact1 = 1.0/(1.0 - exp(-(tk[i+1]-tk[i])/epsize[i]) );

    double part1 = 2.0/(1-4.0*epsize[i]*epsize[i]*rho*rho);
    double part1exp = -(tk[i+1]-tk[i])/epsize[i]-2*rho*tk[i+1];
    part1 *= exp(part1exp);

    double part2 = (2.0*rho*epsize[i])/(1+2*epsize[i]*rho)*exp(-2*rho*tk[i]);

    double part3 = (2*epsize[i]*rho)/(1-2*epsize[i]*rho);
    double part3exp = -2*(tk[i+1]-tk[i])/epsize[i]-2*rho*tk[i];
    part3 *= exp(-part3exp);
    

    double part4exp = -(tk[i+1]-tk[i])/epsize[i]-2*rho*tk[i];
    double part4 = 2*exp(part4exp);
    
    double fact2= part1+part2-part3-part4;
    
    //fprintf(stderr,"\t-> fact1: %f fact2: %f fact1*fact2: %f\n",fact1,fact2,fact1*fact2);
    //fprintf(stderr,"\t-> part1: %f part2: %f part3: %f part4: %f\n",part1,part2,part3,part4);
    P4[i] = log(fact1)+log(fact2);
    fprintf(stderr,"P[4][%d]: %f\n",i,P4[i]);

    //exit(0);
    assert(P4[i]>=0&&P4[i]<=1);
  }
  P4[tk_l-1] = log(2.0*rho/(1.0 + 2.0*rho*epsize[tk_l-1])*exp(-2.0*rho*tk[tk_l-1]));
  assert(P4[tk_l-1]>=0&&P4[tk_l-1]<=1);
}

//this is the appendix version, ordering has been swapped to match, with rho=rho*2, it gives similar results
void ComputeP4(double *tk,int tk_l,double *P4,double *epsize,double rho){
  fprintf(stderr,"\t->[%s] rho: %f\n",__FUNCTION__,rho);
  //rho = rho*2.0;
  for (unsigned i = 0; i < tk_l-1; i++){
    fprintf(stderr,"\t->[%s] epsize[%d]: %f tk[%d+1]: %f tk[%d]: %f tk[i+1]-tk[i]: %f\n",__FUNCTION__,i,epsize[i],i,tk[i+1],i,tk[i],tk[i+1]-tk[i]);
    double fact1 = 1.0/(1.0 - exp(-(tk[i+1]-tk[i])/epsize[i]) );

    double part2 = (rho*epsize[i])/(1+epsize[i]*rho)*exp(rho*tk[i]);

    double part4exp = -(tk[i+1]-tk[i])/epsize[i]-rho*tk[i];
    double part4 = 2*exp(part4exp);

    double part3 = (epsize[i]*rho)/(1-epsize[i]*rho);
    double part3exp = -2*(tk[i+1]-tk[i])/epsize[i]-rho*tk[i];
    part3 *= exp(part3exp);

    
    double part1 = 2.0/((1-epsize[i]*rho)*(1+epsize[i]*rho));
    double part1exp = -(tk[i+1]-tk[i])/epsize[i]-rho*tk[i+1];
    part1 *= exp(part1exp);

    
    
    double fact2= part1+part2-part3-part4;
    
    fprintf(stderr,"\t-> fact1: %f fact2: %f fact1*fact2: %f\n",fact1,fact2,fact1*fact2);
    fprintf(stderr,"\t-> part1: %f part2: %f part3: %f part4: %f\n",part1,part2,part3,part4);
    P4[i] = log(fact1)+log(fact2);
    fprintf(stderr,"P[4][%d]: %f\n",i,P4[i]);

    exit(0);
  }
  P4[tk_l-1] = log(2.0*rho/(1.0 + 2.0*rho*epsize[tk_l-1])*exp(-2.0*rho*tk[tk_l-1]));
}

  
void ComputeP7(double *tk,int tk_l,double *P7,double *epsize,double rho){
  for (unsigned i = 0; i < tk_l - 1; i++){
    P7[i] = 1.0 - exp(-(tk[i+1]-tk[i])*2.0*rho) - exp(-tk[i]*2.0*rho);
    P7[i] -= epsize[i]*2*rho/(1 - epsize[i]*2.0*rho)*exp(-(tk[i+1]-tk[i])/epsize[i]-tk[i]*2.0*rho);
    P7[i] += 1.0/(1.0 - epsize[i]*2.0*rho)*exp(-tk[i]*2.0*rho);
    P7[i] = log(P7[i]);
  }
  unsigned i = tk_l - 1;
  P7[i] = log(1.0 - exp(-2.0*rho*tk[i]));
}
  
void ComputeP0(int tk_l,double *P0,double *P5){ //probability P(T > i)
  P0[0] = P5[0];
  for (unsigned i = 1; i < tk_l; i++)
    P0[i] = P0[i-1]+P5[i];
}
