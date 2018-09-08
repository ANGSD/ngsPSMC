#include <cstdio>
#include <cstring>
#include <vector>
#include <cstdlib>
#include "main_psmc.h"
#include "msArg_toPars.h"
#define MSARG "./ms 2 100 -t 8196 -r 1355 3000000 -l -eN 0.01 0.05 -eN 0.0375 0.5 -eN 1.25 1.0 "

//simple fucntion that only READS msatyle input similar to MSARG
msarg function(char *str){
  std::vector<double> eN1; //times
  std::vector<double> eN2; //sizes
  double theta=-1.0;
  double rho[2];// recombination regionlength
  
  //fprintf(stderr,"*str:%s\n",str);
  char *tok = strtok(str," \n");
  while(tok){
    //    fprintf(stderr,"*tok: %s\n",tok);
    if(strcmp("-t",tok)==0)
      theta=atof(strtok(NULL," \n"));
    else if (strcmp("-r",tok)==0){
      rho[0] = atof(strtok(NULL," \n"));
      rho[1] = atof(strtok(NULL," \n"));
    }
    else if (strcmp("-eN",tok)==0){
      double val[2]; 
      val[0] =  atof(strtok(NULL," \n"));
      val[1] =  atof(strtok(NULL," \n"));
      if(eN1.size()==0&&val[0]!=0.0){
	eN1.push_back(0);
	eN2.push_back(1);
      }
      eN1.push_back(val[0]);
      eN2.push_back(val[1]);

    }

    tok = strtok(NULL," \n");
  }
#if 0
  fprintf(stderr,"theta: %f rho:(%f,%f) eN: ",theta,rho[0],rho[1]);
  for(int i=0;i<eN1.size();i++)
    fprintf(stderr," (%.3f,%.3f)",eN1[i],eN2[i]);
  fprintf(stderr,"\n");
  exit(0);
#endif
  msarg ret;
  ret.eN1 = eN1;
  ret.eN2 = eN2;
  ret.rho[0]=rho[0];   ret.rho[1]=rho[1];
  ret.theta=theta;
  return ret;
}

//plug ms into par. 
//parstruct is in main_psmc.h
//msarg struct is in main_psmc.h
void transform(msarg &ms,psmc_par *par,int winsize){
  fprintf(stderr,"\t->[%s] theta: %f rho:(%f,%f) eN: par->pattern:%s par->par_map:%p winsize:%d\n",__FUNCTION__,ms.theta,ms.rho[0],ms.rho[1],par->pattern,par->par_map,winsize);
  for(int i=0;0&&i<ms.eN1.size();i++)
    fprintf(stderr,"(%.3f,%.3f)\n",ms.eN1[i],ms.eN2[i]);

#if 1
  int *psmc_parse_pattern(const char *pattern, int *n_free, int *n_pars);
  if(par->par_map)
    free(par->par_map);
  par->par_map= psmc_parse_pattern(par->pattern,&par->n_free,&par->n);
   extern int remap_l;
   extern int * remap;
  void make_remapper(psmc_par *pp);
  make_remapper(par);
  //fprintf(stderr,"remap_l: %d eN1.size():%d\n",remap_l,ms.eN1.size());
  assert(ms.eN1.size()==remap_l);

#endif
 
#if 1
  par->times=new double[par->n+1];
  par->params=new double[par->n+1];
  int at =0;
  //  fprintf(stderr,"remap_l:%d eN1.size():%d par->n_free:%d\n",remap_l,ms.eN1.size(),par->n_free);

  for(int i=0;i<remap_l;i++){
    double diff;
    for(int r=0;r<remap[i];r++){
      if((r==0)&&i<remap_l-1){
	diff = (ms.eN1[i+1]-ms.eN1[i])/((double) remap[i]);
	//	fprintf(stderr,"diff:%f\n",diff);
      }
      par->times[at] = ms.eN1[i]+diff*r;
     
      par->params[at++] = ms.eN2[i];
      // fprintf(stderr,"times[%d]:%f\n",at-1,par->times[at-1]);     
      //fprintf(stderr,"params[%d]:%f\n",at,par->params[at]);
      //  fprintf(stderr,"(times:pars)[%d]: (%f,%f)\n",at-1,par->times[at-1],par->params[at-1]);
    }
  }

#endif 

#if 1
  fprintf(stderr,"ms.thteta:%f ms.rho[0]:%f ms.rho[1]:%f\n",ms.theta,ms.rho[0],ms.rho[1]);
  double ngsPsmcTheta=ms.theta/ms.rho[1]/2.0*winsize;
  //rho[0] recombination; rho[1] regionlength
  double ngsPsmcRho = ms.rho[0]/ms.rho[1]/2.0*winsize;

  par->TR[0]=ngsPsmcTheta;
  par->TR[1]=ngsPsmcRho;
  fprintf(stderr,"\t-> ms calculated theta:%f rho:%f\n",par->TR[0],par->TR[1]);
#endif
  //  exit(0);
}



#ifdef __WITH_MAIN__
int main(int argc,char** argv){
  char *msarg= strdup(MSARG);
  function(msarg);
  return 0;
}
#endif
