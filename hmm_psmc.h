/*
  Class that contains the perChr hmm
 */
double lprod(double a,double b);
double lprod(double a,double b,double c,double d);

void ComputeGlobalProbabilities(double *tk,int tk_l,double **P,const double *epsize,double rho);

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
  static double **trans;// tk_l x tk_l
  static double **P;//8xtk_l
  static double *stationary;//tk_l
  static double **nP;//8xtk_l
  int index;
  static int tot_index;
  double pix;
  int tk_l;
  double max_t;
 
  double **PP;//8xtk_l
  
  double *R1,*R2;//tk_l long
  double **fw;//tk_l x nWindows+1
  double **bw;//tk_l x nWindows+1
  double **emis;//tk_l x nWindows+1
  mygltype *gls;//deep copy of the gls for a chr
  std::vector<wins> windows;

  double *workspace;
  double fwllh;
  double bwllh;
  double **baumwelch;
  double qval;
  fastPSMC(){
    pix = -666;
    max_t = 15;
    index=tot_index++;
    has_calc_emissions = 0;
  }
  ~fastPSMC();
  void setWindows(int *pos,int last,int block);
  void printWindows(FILE *fp){
    for(int w=0;w<windows.size();w++)
      fprintf(fp,"win[%d]=(%d,%d)\n",w,windows[w].from,windows[w].to);
  }
  void allocate(int tk_l);
  void calculate_FW_BW_Probs(double *tk,int tk_l,double *epsize,double rho);
  void make_hmm_pre(double *tk,int tk_l,double *epsize,double theta,double rho);
  double make_hmm(double *tk,int tk_l,double *epsize,double theta,double rho);
  void print_emission(const char *fname){
    FILE *fp=NULL;
    fp=fopen(fname,"w");
    assert(fp);
    for(int w=0;w<windows.size();w++){
      for(int i=0;i<tk_l-1;i++){
	fprintf(fp,"%f\t",emis[i][w]);
      }
      fprintf(fp,"%f\n",emis[tk_l-1][w]);
    }
    fclose(fp);
    exit(0);
  }
private:
  int has_calc_emissions;
  void ComputeR2(int v,double **mat,int direction){
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
    if(direction==0)
      fprintf(stderr,"ComputeRs_R2[%d]:\t%f\t%f\t%f\n",v,R2[0],R2[1],R2[2]);
  }
  
  void ComputeR1(int v,double **mat,int direction){
    double addProtect2(double,double);
    R1[tk_l - 1] = log(0);
    for (int i = tk_l - 2; i >= 0 ; i--){
      R1[i] = addProtect2(R1[i+1] , mat[i+1][v]);
    }
    if(0&&direction==0)//0=from start to end, 1=from end to start
      fprintf(stderr,"ComputeRs_R1[%d]:\t%f\t%f\t%f\n",v,R1[0],R1[1],R1[2]);
  }
  
  void ComputeRs(int v,double **mat,int direction){
    ComputeR1(v,mat,direction);
    ComputeR2(v,mat,direction);
  }

};
