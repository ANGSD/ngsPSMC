int psmc_wrapper(args *pars,int block);

typedef struct{
  clock_t t;
  time_t t2;
  double tids[2];
}timer;

timer starttimer();
void stoptimer(timer &t);
