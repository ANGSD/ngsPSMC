#include <vector>		// 
#include <ctime>
#include  <cmath>
#include <ctype.h>
#include "main_psmc.h"
#include "fpsmc.h"
#include <htslib/kstring.h>
extern int nThreads;


double lprod(double a,double b,double c,double d){
  if(isinf(a)||isinf(b)||isinf(c)||isinf(d)){
    //    fprintf(stderr,"Something is inf\n");
    return log(0.0);
  }else{
    // fprintf(stderr,"Something is not inf\n");
    return a+b+c+d;

  }

}


double lprod(double a,double b){
  if(isinf(a)||isinf(b)){
    //    fprintf(stderr,"Something is inf\n");
    return log(0.0);
  }else{
    // fprintf(stderr,"Something is not inf\n");
    return a+b;

  }

}


double lprod(double a,double b,double c){
  if(isinf(a)||isinf(b)||isinf(c)){
    //    fprintf(stderr,"Something is inf\n");
    return log(0.0);
  }else{
    // fprintf(stderr,"Something is not inf\n");
    return a+b+c;

  }

}


#define DEFAULT_PATTERN "4+5*3+4"
// parse a pattern like "4+5*3+4"
// the returned array holds which parameters are linked together
// number of parameters and number of free parameters will be also returned

int *psmc_parse_pattern(const char *pattern, int *n_free, int *n_pars)
{
  fprintf(stderr,"\t-> parsing pattern :\"%s\"\n",pattern);
	char *q, *p, *tmp;
	int top = 0, *stack = (int*)malloc(sizeof(int) * 0x100);
	int *pars_map, k, l, i;
	p = q = tmp = strdup(pattern);
	k = 1;
	while (1) {
		assert(isdigit(*p) || *p == '*' || *p == '+' || *p == '\0'); // allowed characters
		if (*p == '+' || *p == '\0') {
			int is_end = (*p == 0)? 1 : 0;
			*p++ = '\0';
			l = atoi(q); q = p;
			for (i = 0; i < k; ++i) {
				stack[top++] = l;
				assert(top <= 0xff);
			}
			k = 1;
			if (is_end) break;
		} else if (*p == '*') {
			*p = '\0';
			k = atoi(q); // number of repeats
			*p++ = '*'; q = p;
		} else ++p;
	}
	for (k = l = 0; k != top; ++k) l += stack[k];
	*n_pars = l - 1; *n_free = top;
	pars_map = (int*)malloc(sizeof(int) * (*n_pars + 1));
	for (k = i = 0; k != top; ++k)
		for (l = 0; l < stack[k]; ++l)
			pars_map[i++] = k;
	free(tmp); free(stack);
	return pars_map;
}

void setpars( char *fname,psmc_par *pp,int which) {
  fprintf(stderr,"\t-> [%s]:%s which:%d\n",__FUNCTION__,fname,which);
  FILE *fp = NULL;
  fp=fopen(fname,"r");
  if(!fp){
    fprintf(stderr,"\t-> Problem opening file:%s\n",fname);
    exit(0);
  }
  char *buf = new char[fsize(fname)+10];
  memset(buf,0,fsize(fname)+10);
  fread(buf,sizeof(char),fsize(fname),fp);
  fclose(fp);
  char *slashslash[100];
  
  //stupid loop below....
  int n=0;

  //catch first case seperately

  for(int i=0;1&i<strlen(buf)-1;i++){//offset with one so we dont get the last empty output from PSCMC
    if(strncmp(buf+i,"\nRD\t",4)==0){
      assert(i>8);//check that we can plug in dummy values.
      slashslash[n] = buf+i;
      buf[i-1] ='0';
      buf[i-2] ='\t';
      buf[i-3] ='T';
      buf[i-4] ='I';
      buf[i-5] ='\n';
      buf[i-6] = buf[i-7]= '/';
      buf[i-8] ='\n';
      //      slashslash[n++] = buf+i-8;
      // fprintf(stderr,"buf:\'%s\'",slashslash[n]);
      break;
    }
  }
  //  exit(0);

  for(int i=0;i<strlen(buf)-1;i++){//offset with one so we dont get the last empty output from PSCMC
    if(strncmp(buf+i,"\n//\n",4)==0)
      slashslash[n++] = buf+i;
  }
  fprintf(stderr,"\t-> Number of rounds:%d\n",n-2);
  if(which>=n-1){
    fprintf(stderr,"\t-> which higher than number of rounds, setting which to last element");
    which=-1;
  }
  //  fprintf(stderr,"whcich:%d slah:%s\n",which,slashslash[which]);
  char *last = which!=-1?slashslash[which]:slashslash[n-2];
  //  fprintf(stderr,"slashslash[%d]:%s\n",which,last);exit(0);
  char *line = NULL;
  strtok(last,"\n");

  line=strtok(NULL,"\n");
  int IT=-1;
  sscanf(line,"IT\t%d",&IT);
  int RD=-1;
  line=strtok(NULL,"\n"); sscanf(line,"RD\t%d",&RD);
  fprintf(stderr,"\t-> Using round: %d\n",RD);
  double LK=-1;
  line=strtok(NULL,"\n"); sscanf(line,"LK\t%lf",&LK);
  double QD[2]={-1,-1};
  line=strtok(NULL,"\n"); sscanf(line,"QD\t%lf -> %lf",&QD[0],&QD[1]);
  double RI=-1;
  line=strtok(NULL,"\n"); sscanf(line,"RI\t%lf",&RI);
  double TR[2]={-1,-1};
  line=strtok(NULL,"\n"); sscanf(line,"TR\t%lf\t%lf",&TR[0],&TR[1]);
  double MT={-1};
  line=strtok(NULL,"\n"); sscanf(line,"MT\t%lf",&MT);
  double C_pi=-1;
  double n_recomb=-1;
  line=strtok(NULL,"\n"); sscanf(line,"MM\tC_pi: %lf, n_recomb: %lf",&C_pi,&n_recomb);
  std::vector<char *> RS;
  fprintf(stderr,"\t-> IT:%d RD:%d lk:%f qd[0]:%f qd[1]:%f ri:%f tr[0]:%f tr[1]:%f mt:%f c_pi:%f n_rebomc:%f\n",IT,RD,LK,QD[0],QD[1],RI,TR[0],TR[1],MT,C_pi,n_recomb);
  while(((line=strtok(NULL,"\n")))){
    if(line[0]=='R'&&line[1]=='S')
      RS.push_back(line);
    else{
      break;
    }
  }
  fprintf(stderr,"\t-> Number of lines with RS:%lu\n",RS.size());
  char *nline = strdup(line);
  char *tok = strtok(nline,"\n\t ");
  tok = strtok(NULL,"\n\t ");
  pp->pattern=strdup(tok);

  pp->par_map= psmc_parse_pattern(pp->pattern,&pp->n_free,&pp->n);
  assert(RS.size()-1==pp->n);
  pp->params = new double[RS.size()];
  pp->times = new double[RS.size()];
  pp->TR[0] = TR[0];
  pp->TR[1] = TR[1];
  //  fprintf(stderr,"RS:%lu\n",RS.size());
  for(int i=0;i<RS.size();i++){
    int val;
    sscanf(RS[i],"RS\t%d\t%lf\t%lf\t",&val,&pp->times[i],&pp->params[i]);
    //    fprintf(stderr,"PP->params[%d]:%e\n",i,pp->params[i]);
    assert(val==i);
  }
  fprintf(stderr,"\t-> Done reading parameters from file: \'%s\'\n",fname);
  //  exit(0);
}


void readtkfile(psmc_par *pp,const char *fname){
  fprintf(stderr,"\t[%s]\n",__FUNCTION__);
  FILE *fp = NULL;
  if(!(fp=fopen(fname,"rb"))){
    fprintf(stderr,"\t-> Problem writing file: \'%s\'\n",fname);
    exit(0);
  }
  char buf[4096];
  std::vector<double> tk;
  std::vector<double> lambda;
  while(fgets(buf,4096,fp)){
    if(buf[0]=='#')
      continue;
    tk.push_back(atof(strtok(buf,"\n\t ")));
    lambda.push_back(atof(strtok(NULL,"\n\t ")));
  }
  fprintf(stderr,"\t-> Number of tks read from file:%lu\n",tk.size());
  pp->TR[0] = -666.0;
  pp->TR[1] = -666.0;
  pp->params= new double[tk.size()-1];
  pp->times = new double[tk.size()-1];
  for(int i=0;i<tk.size()-1;i++){
    pp->times[i] = tk[i+1];
    pp->params[i] = lambda[i+1];
  }
  std::vector<int> mult;
  int at=0;
  while(at<tk.size()-1){
    int counter =0;
    //  fprintf(stderr,"at:%d params[%d]:%f\n",at,at,pp->params[at]);
    while(counter<tk.size()-1){
      //  fprintf(stderr,"at:%d counter:%d params[%d]:%f params[%d]:%f\n",at,counter,counter,pp->params[at+counter],at,pp->params[at]);
      if(pp->params[at+counter]!=pp->params[at])
	break;
      counter++;
    }
    fprintf(stderr,"counter:%d at:%d\n",counter,at);
    at+=counter;
    mult.push_back(counter);
  }
  for(int i=0;i<mult.size();i++)
    fprintf(stderr,"mult: %d %d\n",i,mult[i]); 
  //  exit(0)  ;
  std::vector<int> mult2;
  kstring_t kstr;kstr.s=NULL;kstr.l=0;kstr.m=0;
  for(int at=0;at<mult.size();at++){
    int counter =0;
    //    fprintf(stderr,"at:%d params[%d]:%f\n",at,at,pp->params[at]);
    while(counter+at<mult.size()){
      //fprintf(stderr,"params[%d]:%f params[%d]:%f\n",counter,pp->params[at+counter],at,pp->params[at]);
      if(mult[at+counter]!=mult[at])
	break;
      counter++;
    }
    assert(counter>0);
    if(counter>1)
      ksprintf(&kstr,"%d*%d",counter,mult[at]);
    else
      ksprintf(&kstr,"%d",mult[at]);
    at+=counter-1;
    mult2.push_back(counter);
    if(at<mult.size()-1)
      ksprintf(&kstr,"+");
  }
  
  for(int i=0;0&i<mult2.size();i++)
    fprintf(stderr,"mult2: %d %d\n",i,mult2[i]);
  pp->pattern=strdup(kstr.s);
  free(kstr.s);
  //  fprintf(stderr,"buf:%s\n",kstr.s);
  //exit(0);
  fprintf(stderr,"\t[%s] done\n",__FUNCTION__);
}

int doGlStyle =0;

args * getArgs(int argc,char **argv){
  args *p = new args;
  p->chooseChr=NULL;
  p->start=p->stop=-1;
  p->maxIter=1e2;
  p->tole=1e-6;
  p->nSites =0;
  p->fname = NULL;
  p->onlyOnce = 0;
  p->seed =1;
  p->blocksize = 100;//default 100bp
  p->par =(psmc_par*) calloc(1,sizeof(psmc_par));
  p->tkfile = NULL;
  p->RD = -1;
  p->nChr = -1;
  p->nThreads =1;
  p->doQuad =0;
  p->smartsize =0;
  p->outname = strdup("output");
  p->fres = p->flog = NULL;
  p->init = -1;
  char *inffilename=NULL;
  if(argc==0)
    return p;

  while(*argv){
    //    fprintf(stderr,"%s\n",*argv);
    if(!strcasecmp(*argv,"-tole"))
      p->tole = atof(*(++argv));
    else if(!strcasecmp(*argv,"-smartsize"))
      p->smartsize = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-maxIter"))
      p->maxIter = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-winSize"))
      p->blocksize = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-RD"))
      p->RD = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-nThreads"))
      p->nThreads = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-doGlStyle"))
      doGlStyle = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-nIter"))
      p->nIter = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-p"))
      p->par->pattern =  strdup(*(++argv));
    else  if(!strcasecmp(*argv,"-out"))
      p->outname =  strdup(*(++argv));
    else  if(!strcasecmp(*argv,"-tkfile"))
      p->tkfile =  strdup(*(++argv));
    else  if(!strcasecmp(*argv,"-nSites"))
      p->nSites = atol(*(++argv));
    else  if(!strcasecmp(*argv,"-seed"))
      p->seed = atol(*(++argv));
    else  if(!strcasecmp(*argv,"-infile"))
      inffilename = strdup(*++argv);
    else  if(!strcasecmp(*argv,"-init"))
      p->init = atof(*++argv);
   else  if(!strcasecmp(*argv,"-nChr"))
      p->nChr = atoi(*++argv);
    else  if(!strcasecmp(*argv,"-doLinear")){
      int mytmp = atoi(*++argv);
      if(mytmp ==1)
	p->doQuad = 0;
      else
	p->doQuad =1;
    }else  if(!strcasecmp(*argv,"-r")){
      p->chooseChr = get_region(*(++argv),p->start,p->stop);
      if(!p->chooseChr)
	return NULL;
    }
    else{

      p->fname = *argv;

    }
    argv++;
  }
  p->perc = perpsmc_init(p->fname,p->nChr);
  extern int doQuadratic;
  doQuadratic = p->doQuad;
  if(inffilename)
    setpars(inffilename,p->par,p->RD);

  //adjust theta:
  p->par->TR[0] = p->par->TR[0]/2.0;
  fprintf(stderr,"\t-> p->perc->version:%d (one is gls, otherwise fasta)\n",p->perc->version);
  if(p->perc->version==1){
    fprintf(stderr,"\t-> Adjusing theta with blocksize: %d\n",p->blocksize);
    p->par->TR[0] = p->par->TR[0]/(1.0*p->blocksize);
  }
  
  if(p->seed==0)
    p->seed = time(NULL);
  srand48(p->seed);
  fprintf(stderr,"\t-> args: tole:%f maxiter:%d chr:%s start:%d stop:%d fname:%s seed:%ld winsize:%d RD:%d nThreads:%d doLinear:%d doGlStyle:%d\n",p->tole,p->maxIter,p->chooseChr,p->start,p->stop,p->fname,p->seed,p->blocksize,p->RD,p->nThreads,p->doQuad,doGlStyle);
  //  fprintf(stderr,"par:%p par->pattern:%p DEFAULT_PATTERN:%s\n",p->par,p->par->pattern,DEFAULT_PATTERN);
  if(p->tkfile)
    readtkfile(p->par,p->tkfile);
  if(p->par->pattern==NULL)
    p->par->pattern = strdup(DEFAULT_PATTERN);
  //  fprintf(stderr,"par:%p par->pattern:%p DEFAULT_PATTERN:%s\n",p->par,p->par->pattern,DEFAULT_PATTERN);
  if(p->par->pattern!=NULL)
    p->par->par_map = psmc_parse_pattern(p->par->pattern, &p->par->n_free, &p->par->n);
  
  char tmp[1024];
  snprintf(tmp,1024,"%s.log",p->outname);
  fprintf(stderr,"\t-> Writing file: \'%s\'\n",tmp);
  if(fexists(tmp)){
    fprintf(stderr,"\t-> File exists, will exit\n");
    exit(0);
  }
  p->flog = fopen(tmp,"w");
  assert(p->flog!=NULL);
  snprintf(tmp,1024,"%s.res",p->outname);
  fprintf(stderr,"\t-> Writing file: \'%s\'\n",tmp);
  if(fexists(tmp)){
    fprintf(stderr,"\t-> File exists, will exit\n");
    exit(0);
  }
  p->fres = fopen(tmp,"w");
  
  nThreads = p->nThreads;
  if(p->init!=-1)
    for(int i=0;i<p->par->n+1;i++)
      p->par->params[i] = p->init;
      
  return p;
}

//made a seperate function for this. Im assuming our args will contain allocated data at some point.
void destroy_args(args *p){
  perpsmc_destroy(p->perc);
  fclose(p->flog);
  fclose(p->fres);
  delete p;
}


//simple function 
int main_psmc(int argc, char **argv){
  fprintf(stderr,"\t-> we are in file: %s function: %s line:%d\n",__FILE__,__FUNCTION__,__LINE__);

  //we loop over the single chromosomes
#ifdef __SHOW_TIME__
  clock_t t=clock();
  time_t t2=time(NULL);
#endif
  args *pars = getArgs(argc,argv);
#ifdef __SHOW_TIME__
  fprintf(stderr, "\t[TIME] cpu-time used =  %.2f sec for reading arguments\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[Time] walltime used =  %.2f sec for reading arguments\n", (float)(time(NULL) - t2));  
#endif
  for(int i=0;0&&i<pars->par->n+1;i++)
    fprintf(stderr,"%d) tk:%f lambda:%f\n",i,pars->par->times[i],pars->par->params[i]);
  if(!pars)
    return 0;
  //this will printout the header
  writepsmc_header(stderr,pars->perc);

  if(1){
    assert(pars->flog!=NULL);
    psmc_wrapper(pars,pars->blocksize);
  }else{
    //below is old printout, keeping for reference
    for(myMap::iterator it=pars->perc->mm.begin();it!=pars->perc->mm.end();++it){
      //set perchr iterator, if pars->chooseChr, then we have only use a single chr
      it = pars->chooseChr?iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop):iter_init(pars->perc,it->first,pars->start,pars->stop);
      
      //print out the chromosome position and the two gls
      for(size_t s=pars->perc->first;0&&s<pars->perc->last;s++)
	fprintf(stdout,"%s\t%d\t%e\t%e\n",it->first,pars->perc->pos[s]+1,pars->perc->gls[2*s],pars->perc->gls[2*s+1]);
      
      if(pars->chooseChr!=NULL)
	break;
    }
  }
  destroy_args(pars);
  return 0;
}
