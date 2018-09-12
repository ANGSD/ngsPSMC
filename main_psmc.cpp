#include <vector>		// 
#include <ctime>
#include  <cmath>
#include <ctype.h>
#include "msArg_toPars.h"
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
	//	fprintf(stderr,"psmc_parse_pattern: n_pars:%d\n",*n_pars);
	pars_map = (int*)malloc(sizeof(int) * (*n_pars + 1));//<-bug
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

  for(int i=0;1&&i<strlen(buf)-1;i++){//offset with one so we dont get the last empty output from PSCMC
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
  if(pp->pattern)
    free(pp->pattern);
  pp->pattern=strdup(tok);
  if(pp->par_map)
    free(pp->par_map);
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
  free(nline);
  delete [] buf;
}


int doGlStyle =0;

args * getArgs(int argc,char **argv,int dontprint){
  args *p = new args;
  p->flog = NULL;
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
  p->RD = -1;
  p->nChr = -1;
  p->nThreads =1;
  p->doQuad =0;
  p->smartsize =0;
  p->outname = strdup("output");
  p->flog = NULL;
  p->init = -1;
  p->msstr = NULL;
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
    else  if(!strcasecmp(*argv,"-out")){
      if(p->outname)
	free(p->outname);
      p->outname =  strdup(*(++argv));
    }
    else  if(!strcasecmp(*argv,"-ms"))
      p->msstr =  strdup(*(++argv));
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
  if(inffilename){
    setpars(inffilename,p->par,p->RD);
    free(inffilename);
  }
  
  if(p->seed==0)
    p->seed = time(NULL);
  srand48(p->seed);
  fprintf(stderr,"\t-> args: tole:%f maxiter:%d chr:%s start:%d stop:%d\n\t-> fname:\'%s\' seed:%ld winsize:%d RD:%d nThreads:%d doLinear:%d doGlStyle:%d -nChr:%d -ms:\'%s\'\n",p->tole,p->maxIter,p->chooseChr,p->start,p->stop,p->fname,p->seed,p->blocksize,p->RD,p->nThreads,p->doQuad,doGlStyle,p->nChr,p->msstr);
  //  fprintf(stderr,"par:%p par->pattern:%p DEFAULT_PATTERN:%s\n",p->par,p->par->pattern,DEFAULT_PATTERN);
 
  if(p->par->pattern==NULL)
    p->par->pattern = strdup(DEFAULT_PATTERN);

  if(p->par->pattern!=NULL){
    if(p->par->par_map)
      free(p->par->par_map);
    p->par->par_map = psmc_parse_pattern(p->par->pattern, &p->par->n_free, &p->par->n);
  }
  if(dontprint!=1){
    char tmp[1024];
    snprintf(tmp,1024,"%s.log",p->outname);
    fprintf(stderr,"\t-> Writing file: \'%s\'\n",tmp);
    if(0&&fexists(tmp)){
      fprintf(stderr,"\t-> File exists, will exit\n");
      destroy_args(p);
      return NULL;
    }
    p->flog = fopen(tmp,"w");
    assert(p->flog!=NULL);
  }
  nThreads = p->nThreads;
  if(p->init!=-1)
    for(int i=0;i<p->par->n+1;i++)
      p->par->params[i] = p->init;
      
  return p;
}

//made a seperate function for this. Im assuming our args will contain allocated data at some point.
void destroy_args(args *p){
  if(p->msstr)
    free(p->msstr);
  perpsmc_destroy(p->perc);
  if(p->flog)
    fclose(p->flog);
  if(p->outname)
    free(p->outname);
  if(p->par->par_map)
    free(p->par->par_map);
  if(p->par->pattern)
    free(p->par->pattern);
  if(p->par->params)
    delete [] p->par->params;
  if(p->par->times)
    delete [] p->par->times;
  if(p->par)
    free(p->par);
  delete p;
}

 
//simple function 
int main_psmc(int argc, char **argv){
  fprintf(stderr,"\t-> we are in file: %s function: %s line:%d\n",__FILE__,__FUNCTION__,__LINE__);

  timer t = starttimer();
  //we loop over the single chromosomes

  args *pars = getArgs(argc,argv,0);
  if(!pars)
    return 0;

  if(!pars)
    return 0;
  //this will printout the header
  writepsmc_header(stderr,pars->perc,1);

  assert(pars->flog!=NULL);
  psmc_wrapper(pars,pars->blocksize);
#if 0
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

#endif
  destroy_args(pars);
  stoptimer(t);
  fprintf(stdout,"MM\ttotaltime(wall(min),cpu(min)):(%f,%f) \n",t.tids[1],t.tids[0]);
  return 0;
}
