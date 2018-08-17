#include <cstdio>
#include <cstring>
#include <vector>
#include <cstdlib>
#include "msArg_toPars.h"
#define MSARG "./ms 2 100 -t 8196 -r 1355 3000000 -l -eN 0.01 0.05 -eN 0.0375 0.5 -eN 1.25 1.0 "


msarg function(char *str){
  std::vector<double> eN1; //times
  std::vector<double> eN2; //sizes
  eN1.push_back(0);
  eN2.push_back(1);
  double theta=-1.0;
  double rho[2];
  
  //fprintf(stderr,"*str:%s\n",str);
  char *tok = strtok(str," \n");
  while(tok){
    //    fprintf(stderr,"*tok: %s\n",tok);
    if(strcmp("-t",tok)==0)
      theta=atof(strtok(NULL," \n"));
    else if (strcmp("-t",tok)==0){
      rho[0] = atof(strtok(NULL," \n"));
      rho[1] = atof(strtok(NULL," \n"));
    }
    else if (strcmp("-eN",tok)==0){
      eN1.push_back( atof(strtok(NULL," \n")));
      eN2.push_back( atof(strtok(NULL," \n")));
    }

    tok = strtok(NULL," \n");
  }
#if 0
  fprintf(stderr,"theta: %f rho:(%f,%f) eN: ");
  for(int i=0;i<eN1.size();i++)
    fprintf(stderr," (%.3f,%.3f)",eN1[i],eN2[i]);
  fprintf(stderr,"\n");
#endif
  msarg ret;
  ret.eN1 = eN1;
  ret.eN2 = eN2;
  ret.rho[0]=rho[0];   ret.rho[1]=rho[1];
  ret.theta=theta;
  return ret;
}





#ifdef __WITH_MAIN__
int main(int argc,char** argv){
  char *msarg= strdup(MSARG);
  function(msarg);
  return 0;
}
#endif
