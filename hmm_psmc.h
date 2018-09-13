/*
  Class that contains the perChr hmm
 */
double lprod(double a,double b);
double lprod(double a,double b,double c,double d);
#define DOTRANS 1

#define PSMC_T_INF 1000.0
struct wins{
  int from;//inclusive
  int to;//inclusive
};
/*
  baumwelch,PP is not in logspace

 */
class fastPSMC {
public:
  //shared between all threads
  double **trans;
  double **P;

  int index;
  static int tot_index;
  static char *outnames;
  double pix;
  int tk_l;
  double max_t;
 
  double **PP;
  double **nP;
  double *stationary,*R1,*R2;//tk_l long
  double **fw;//tk_l x nWindows+1
  double **bw;//tk_l x nWindows+1
  double **emis;//tk_l x nWindows+1
  double *gls;//deep copy of the gls for a chr
  std::vector<wins> windows;

  double *workspace;
  double fwllh;
  double bwllh;
  double **baumwelch;
  double qval;
  fastPSMC(){
    trans = NULL;
    pix = -666;
    max_t = 15;
    index=tot_index++;
  }
  ~fastPSMC();
  void setWindows(double *gls_a,int *pos,int last,int block);
  void printWindows(FILE *fp){
    for(int w=0;w<windows.size();w++)
      fprintf(fp,"win[%d]=(%d,%d)\n",w,windows[w].from,windows[w].to);
  }
  void allocate(int tk_l);
  void calculate_FW_BW_Probs(double *tk,int tk_l,double *epsize,double rho);

  double make_hmm(double *tk,int tk_l,double *epsize,double theta,double rho);
  
private:
  void ComputeR2(int v,double **mat){
    double addProtect3(double,double,double);
    double tmp = log(0);
    for (unsigned i = 0; i < tk_l ; i++){
      double p1=  lprod(tmp,P[5][i]);
      double p2= lprod(mat[i][v],P[6][i]);
      double p3 = lprod(R1[i],P[7][i]);
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
      R1[i] = addProtect2(R1[i+1] , mat[i+1][v]);
    }
    
  }
  
  void ComputeRs(int v,double **mat){
    ComputeR1(v,mat);
    ComputeR2(v,mat);
  }

};
