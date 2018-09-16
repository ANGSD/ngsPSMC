#include <cmath>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "splineEPsize.h"

#ifdef __WITH_MAIN__
//remember to subtract by one
void setTk(int n, double *t, double max_t, double alpha, double *inp_ti){
  //  assert(inp_ti!=NULL);
  fprintf(stderr,"[%s] (n,tk,max_t,alpha,inp_ti)=(%d,%p,%f,%f,%p)\n",__FUNCTION__,n,t,max_t,alpha,inp_ti);
  int k;
  if (inp_ti == 0) {
    double beta;
    beta = log(1.0 + max_t / alpha) / n; // beta controls the sizes of intervals
    for (k = 0; k <= n; ++k)
      t[k] = alpha * (exp(beta * k) - 1);
    t[n+1] = max_t;
    //    t[n] = PSMC_T_INF; // the infinity: exp(PSMC_T_INF) > 1e310 = inf
  } else {
    memcpy(t, inp_ti, n * sizeof(double));
  }
}
#endif


double splineEPSize::Poly(int degree, double *coef, double x){
  double val = 0;
  for(int i = 0; i < degree+1; i++){
    val += coef[i]*pow(x, degree - i);
  }
  return val;
}

void splineEPSize::printAll(FILE *fp,double *epsize){
  fprintf(fp,"(tk_l,numInt,degree)\t%d %d %d\n",tk_l,nsplines,degree);
  for(int i=0;i<nsplines;i++){
    fprintf(fp,"sp\t%f\t%f",tk[Tk[i]],tk[Tk[i+1]]);
    for(int j=0;j<degree+1;j++)
      fprintf(fp,"\t%f",spline[i][j]);
    fprintf(fp,"\n");
  }
  int at=0;
  for(int i=0;i<tk_l;i++)
    fprintf(fp,"xy[%d][%d]\t%f\t%f\n",at++,i,tk[i],epsize[i]);
  for(int i=0;i<nsplines+1;i++)
    fprintf(fp,"fd[%d]\t%f\t%f\t%f\n",at++,tk[Tk[i]],fv[i],dv[i]);
}

void splineEPSize::computeSpline(){
  for(int i = 0; i < nsplines; i++){
    double t0 = tk[Tk[i]];
    //    fprintf(stderr,"t0:%f Tk[i]:%d\n",t0,Tk[i]);
    double t1 = tk[Tk[i+1]];
    double v0 = fv[i];
    double v1 = fv[i+1];
    double g0 = dv[i];
    double g1 = dv[i+1];
    double a0, a1, a2, a3;

    double tmpval = pow(t0 - t1, 3);
    //    fprintf(stderr,"tmpval: %f t0:%f t1:%f\n",tmpval,t0,t1);
    a0 = ((g0 + g1)*(t0 - t1) - 2*v0 + 2*v1)/tmpval;
    a1 = (-g0*(t0 - t1)*(t0 + 2*t1) + g1*(-2*t0*t0 + t0*t1 + t1*t1) + 3*(t0 + t1)*(v0 - v1))/pow(t0 - t1, 3);
    a2 = (g1*t0*(t0 - t1)*(t0 + 2*t1) - t1*(g0*(-2*t0*t0 + t0*t1 + t1*t1) + 6*t0*(v0 - v1)))/pow(t0 - t1, 3);
    a3 = (t1*(t0*(-t0 + t1)*(g1*t0 + g0*t1) - t1*(-3*t0 + t1)*v0) + t0*t0*(t0 - 3*t1)*v1)/pow(t0 - t1, 3);
    //fprintf(stderr,"a0:%f, a1:%f, a2:%f, a3:%f\n",a0,a1,a2,a3);
    double tmp1[4]={a0, a1, a2, a3};
    double tmp2[3]={3*a0, 2*a1, a2};


    // fprintf(stderr,"ls: %f rs: %f diff:%e\n",Poly(3,tmp1,t1),v1,Poly(3,tmp1,t1)-v1);exit(0);
    assert(fabs(Poly(3, tmp1, t0) - v0)<1e-6);
    assert(fabs(Poly(3, tmp1, t1) - v1)<1e-6);
    assert(fabs(Poly(2, tmp2, t0) - g0)<1e-6);
    assert(fabs(Poly(2, tmp2, t1) - g1)<1e-6);


    spline[i][0] = a0;
    spline[i][1] = a1;
    spline[i][2] = a2;
    spline[i][3] = a3;
    
  }
  //  fprintf(stderr,"[%s] done\n",__FUNCTION__);fflush(stderr);
}

void splineEPSize::computeEPSize(double *epsize){//FIXME epsize should be of length tk_l+1?
  //  fprintf(stderr,"[%s]\n",__FUNCTION__);fflush(stderr);
  double t0,t1; //start and end of integration interval
  int si;//<-  splinenumber
  si = 0;
  for (int i = 0; i < tk_l-1 ; i++){
    t0 = tk[i];
    t1 = tk[i+1];
    if (t0 >= tk[Tk[si + 1]]){
      si += 1;
      if(si==nsplines)
	break;
    }
    //    fprintf(stderr,"i:%d si:%d \n",i,si);
    double tmp1[5]={spline[si][0]/4, spline[si][1]/3, spline[si][2]/2, spline[si][3], 0};//check
    epsize[i] = (Poly(4,tmp1, t1) - Poly(4,tmp1, t0))/(t1-t0);
    
  }
  for(int i=Tk[nsplines];i<tk_l;i++)
    epsize[i] = fv[nsplines];

#if 0
  for(int i=0;i<tk_l;i++){
    if(epsize[i]<0){
      printAll(stderr,epsize);
      exit(0);
    }
  }

#endif
}



//tofull zero,moves from pars to pars2 assuming 
//tofull one, takes the
void splineEPSize::convert(const double *pars,double *pars2,int tofull){
  if(tofull==0){
    setfd(pars);
    computeSpline();
    computeEPSize(pars2);
  }else{
    getfd(pars2);
  }
}


#ifdef __WITH_MAIN__

int main(){
  splineEPSize obj(14,3,7,23.889078);
  //  obj.fillit();
  double pars[obj.n_free];
  for(int i=0;i<obj.n_free;i++)
    if(i<obj.nsplines+1)
      pars[i] = i+drand48()*5;
    else
      pars[i] = drand48()*5-2.5;
  obj.setfd(pars);
  obj.computeSpline();
  
  double *epsize= new double[obj.tk_l];
  obj.computeEPSize(epsize);
  obj.printAll(stdout,epsize);
  delete [] epsize;
  /*  
      splineEPSize obj(8,2);
      
      obj.tk[0]=0.0;
      for(int i=1;i<8;i++)
      obj.tk[i]=i+drand48();
      obj.spline[1][0] = 888;
      
      obj.Tk[0]=0;
      obj.Tk[1]=3;
      obj.Tk[2]=7;
      
      obj.fv[0] = 4.4;
      obj.fv[1] = 0.1;
      obj.fv[2] = 10.7;
      obj.dv[0] = 2.1;
      obj.dv[1] = -3;
      obj.dv[2] = 6;
      
      
      obj.computeSpline();
      double epsize[8+1];
      obj.computeEPSize(epsize);
      obj.printAll(stdout,epsize);
      return 0;
  */
}

#endif
