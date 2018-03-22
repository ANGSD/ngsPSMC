/*
  Class that contains the perChr hmm
 */
double lprod(double a,double b);
double lprod(double a,double b,double c,double d);
#define DOTRANS 1

double qkFunction(unsigned i, double pix, unsigned numWind,double **nP,double **PP,int);
void setTk(int n, double *t, double max_t, double alpha, double *inp_ti);
void setEPSize(double *ary,int l,double *from_infile);
#define PSMC_T_INF 1000.0

struct wins{
  int from;//inclusive
  int to;//inclusive
};

//void calculate_emissions(double *tk,int tk_l,double *gls,std::vector<wins> &windows,double **emis);
/*
  baumwelch,PP is not in logspace

 */
class fastPSMC {
public:
  int index;
  static int tot_index;
  static char *outnames;
  double pix;
  int tk_l;
  double max_t;
  double **P;
  double **PP;
  double **nP;
  double *stationary,*R1,*R2;//tk_l long
  double **fw;//tk_l x nWindows+1
  double **bw;//tk_l x nWindows+1
  //  double **pp;//tk_l x nWindows+1
  double **emis;//tk_l x nWindows+1
  double *gls;//deep copy of the gls for a chr
  std::vector<wins> windows;
  double **trans;
  double *workspace;
  double fwllh;
  double bwllh;

  double **baumwelch;
  double qval;
  fastPSMC(){
    trans = NULL;
    pix = -666;
    max_t = 15;
    //rho = 0.207;
    index=tot_index++;
  }
  /*
    double fwllh();
    double bwllh();
  */
  void setWindows(double *gls_a,int *pos,int last,int block);
  void printWindows(FILE *fp){
    //print indices for endpoint of windows
    for(int w=0;w<windows.size();w++)
      fprintf(stdout,"win[%d]=(%d,%d)\n",w,windows[w].from,windows[w].to);
  //print out data:
#if 0
  for(int w=0;0&&w<windows.size();w++)
    for(int p=windows[w].from;p<windows[w].to;p++)
      fprintf(stdout,"%d\t%d\t%f\t%f\n",p,w,gls[2*p],gls[2*p+1]);
#endif
  }
  void allocate(int tk_l);
  //  void ComputePii(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary);
  //  void calculate_stationary(int tk_l,double *results,double **P);
  void calculate_FW_BW_PP_Probs(double *tk,int tk_l,double *epsize,double rho);
  //  void calculate_FW_BW_PP_Probs();
  double make_hmm(double *tk,int tk_l,double *epsize,double theta,double rho);
  
private:
  void ComputeR2(int v,double **mat){
    double addProtect3(double,double,double);
    double tmp = log(0);
    for (unsigned i = 0; i < tk_l ; i++){
      double p1=  lprod(tmp,P[5][i]);
      double p2= lprod(mat[i][v],P[6][i]);
      double p3 = lprod(R1[i],P[7][i]);
      //      fprintf(stderr,"p1:%f\tp2:%f\tp3:%f\n",p1,p2,p3);
      R2[i] = addProtect3(p1,p2,p3);
      if(std::isnan(R2[i])){
	fprintf(stderr,"[hmm_psmc.h:computeR2] R2[%d] evaluates to NaN p5:%f tmp:%f p1:%f p2:%f p3:%f\n",i,P[5][i],tmp,p1,p2,p3);
	exit(0);
      }
      tmp = R2[i];
    }
  }
  
  void ComputeR1(int v,double **mat){
    double addProtect2(double,double);
    R1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--){
      //      fprintf(stderr,"R1[%d]:%f mat[%d][%d]:%f\n",i+1,R1[i+1],i+1,v,mat[i+1][v]);
      R1[i] = addProtect2(R1[i+1] , mat[i+1][v]);
      //fprintf(stderr,"R1[%d]:%f\n",i,R1[i]);
    }
    //exit(0);
  }
  
  void ComputeRs(int v,double **mat){
    ComputeR1(v,mat);
    ComputeR2(v,mat);
  }

};
