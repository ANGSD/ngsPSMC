
class splineEPSize{
public:
  int tk_l;//number of time points in *tk
  int nsplines;//number of splines
  int n_free;
  int degree;//spline degree
  double *tk;//time points
  int *Tk;//time points in *tk where spline values are set
  double *fv, *dv;//time points, splineepsize value, derivative value of slpineepsize value
  double **spline;//cubic spline: four coefficients for each interval
  double alpha;
  double max_t;
  splineEPSize(int nsplines_arg,int pointsperinterval,  int ppl_arg){
    //    fprintf(stderr,"tk:%d intNum:%d nlast:%d\n",tk_l_arg,intNum_arg,nlast_arg);
    //constructor
    alpha=0.01;
    max_t=22;
    nsplines=nsplines_arg;

    tk_l = ppl_arg+nsplines*pointsperinterval+nsplines+1;
    degree=3;
    tk = new double[tk_l];
    Tk = new int[nsplines+1];
    fv = new double[nsplines+1];
    dv = new double[nsplines+1];
    spline = new double *[nsplines];
    n_free = 2*(nsplines+1);
    for(int i = 0; i < nsplines; i++){
      //      fprintf(stderr,"allcing spline[%d][%d]\n",i,degree+1);
      spline[i] = new double[degree + 1];
    }
    void setTk(int n, double *t, double max_t, double alpha, double *inp_ti);
    setTk(tk_l,tk,max_t,alpha,NULL);
    int at=0;
    for(int i=0;i<nsplines+1;i++){
      Tk[i]=i*(pointsperinterval+1);
    }
  }
  void printAll(FILE *fp,double *epsize);
  double Poly(int degree, double *coef, double x);
  void computeEPSize(double *epsize);
  void computeSpline();
  void setfd(double *ary){
    for(int i=0;i<nsplines+1;i++)
      fv[i]=ary[i];
    for(int i=0;i<nsplines+1;i++)
      dv[i]=ary[nsplines+1+i];
  }
  ~splineEPSize(){
    for(int i=0;i<nsplines;i++)
      delete [] spline[i];
    delete[] tk;
    delete [] Tk;
    delete [] fv;
    delete [] dv;
    delete [] spline;

  }
  void fillit(){
    srand48(100);
    for(int i=0;i<nsplines+1;i++){
      fv[i]=drand48()*1+1;
      dv[i]=drand48()*1-0.5;
    }

  }
private:
};
