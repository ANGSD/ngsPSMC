#include <ctime>
#include <vector>
#include <cassert>
#include <cmath>
#include <htslib/kstring.h> //<- included directly in this source due to some string errors in latest htslib 17sep2019
#include "psmcreader.h"
#include "main_psmc.h"
#include "hmm_psmc.h"
#include "compute.h"

int fastPSMC::tot_index=0;
double **fastPSMC::trans=NULL;
double **fastPSMC::P=NULL;
double *fastPSMC::stationary=NULL;
double **fastPSMC::nP = NULL;
//#define __SHOW_TIME__

extern int doQuadratic;

/*
  Calculate stationary distrubution
  tk array of length tk_l
  lambda array effective population sizes
  both has length tk_l
  
  stationary distribution will be put in results array, also of length tk
  
  stationary(i) = exp(-sum_{j=0}^{i-1}{tau_j/lambda_j}*P2[void])
*/

void calculate_stationary(int tk_l,double *results,double **P){
  results[0] = P[2][0];//fix this
  for(int i=1;i<tk_l;i++){
    results[i]  = P[2][i]+P[0][i-1];
  }

  double tmp=0;
  for(int i=0;i<tk_l;i++)
    tmp += exp(results[i]);

  assert(fabs(1-tmp)<1e-6);//check that it sums to one
}

/*
  k: is the state k
  pix: is the perchromsome value of the forward probs
  numwin: is not being used currently
  nP: is the 'new' Pi's (the ones computed by compuateGlobalProbablies)
  PP: is the old PPi's whic correspond to expectation. 

  returnvalue is in logspace
 */

/*
  if something doesnt work check this formula:

  trans<-(exp(read.table("transitions.txt")));
fw<-t(exp(read.table("fw")))
bw<-t(exp(read.table("bw")))
emis<-t(exp(read.table("emis")))
stat<-exp(scan("stationary"))
P<-exp(read.table("P.txt"));colnames(P)<-paste0("P",0:7);
PP<-exp(read.table("PP.txt"));colnames(PP)<-paste0("PP",0:7);
##emis[,1:10]
##fw[,1:10]



 */

double qkFunction(unsigned k, double pix, unsigned numWind,double **nP,double **PP,int tk_l,double &esum){
  //  fprintf(stderr,"k:%u pix:%f numWind:%u tk_l:%d\n",k,pix,numWind,tk_l);
  /*
  //This block is needed if eimission probabilities depend on estimated parameters, e.g. on time disctretisation 
    for (unsigned l = 1; l < numWind + 1; l++) 
    qi += log(emis[K][l])*fw[K][l]*bw[K][l];
    qi /= pix;
  */
  double qi[7];
  double expec[7];
  double npfac[7];

  
  //i follows PP index.
  for(int i=1;i<8;i++){
    expec[i-1] = exp(lprod(PP[i][k],-pix));
    npfac[i-1] = nP[i][k];
    if(i!=2&&i!=5)
      esum +=expec[i-1];

    if(std::isinf(npfac[i-1]))
      qi[i-1] = 0;
    else
      qi[i-1] = npfac[i-1]*expec[i-1];
  
  }
  
  double ret = 0;
  for(int i=0;i<7;i++) 
    if(i==6)
      ret += qi[i];

  return ret;
}
/*
  tk: intervals, fixed
  tk_l: length of tk
  P: matrix where results will be plugged in
  epsize: the effective populationsize (lambda)
  rho: rho, or theta to rho value.
 */

void ComputeGlobalProbabilities(double *tk,int tk_l,double **P,const double *epsize,double rho){
  ComputeP1(tk,tk_l,P[1],epsize,rho);
  ComputeP5(tk,tk_l,P[5],epsize);
  ComputeP6(tk,tk_l,P[6],epsize,rho);
  ComputeP2(tk_l,P[2],P[5]);
  ComputeP3(tk,tk_l,P[3],epsize,rho);
  ComputeP4(tk,tk_l,P[4],epsize,rho);
  ComputeP7(tk,tk_l,P[7],P[3],epsize,rho);
  ComputeP0(tk_l,P[0],P[5]);

  for(int p=0;1&&p<8;p++){
    for(int i=0;i<tk_l;i++){
      //      fprintf(stderr,"P[%d][%d]: %f\n",p,i,P[p][i]);
      if(!(P[p][i]<=0)){
	fprintf(stderr,"\t->[%s] p:%d i:%d val:%f\n",__FUNCTION__,p,i,P[p][i]);
	exit(0);
      }
    }

  }

}
//linear
double qFunction_inner(double *tk,int tk_l,const double *epsize,double rho,double pix,int numWind,double **nP,double **PP){
  assert(0!=1);
  //  ComputeGlobalProbabilities(tk,tk_l,nP,epsize,rho);
  double Q = 0;
  double esum =0;
  for (unsigned i = 0; i < tk_l; i++){
    double tmpQ = qkFunction(i, pix,numWind,nP,PP,tk_l,esum);
    Q += tmpQ;
  }
  if(fabs(numWind-1-esum)>0.5)
    fprintf(stderr,"\t-> POTENTIAL PROBLEM ESUM:%f Q:%f numWind:%d\n",esum,Q,numWind);
  return Q;
  
}

//quadratic
double qFunction_inner2(double *tk,int tk_l,const double *epsize,double rho,double pix,int numWind,double **nP,double **baumwelch,double **trans){
  
  double newstationary[tk_l];
  calculate_stationary(tk_l,newstationary,nP);

  
  double Q = 0;
  for (unsigned i = 0; i < tk_l; i++){
    if(baumwelch[tk_l][i]!=0.0){
      double tmpQ=newstationary[i]*baumwelch[tk_l][i];
      Q += tmpQ;
      if(1&&std::isnan(tmpQ)){
	fprintf(stderr,"\t->[qFunction_inner2] Q[%d]:%e newstat:%e baumwel:%e\n",i,tmpQ,newstationary[i],baumwelch[tk_l][i]);
      exit(0);
      }
    }
  }
  
  for (unsigned i = 0; i < tk_l; i++){
    for (unsigned j = 0; j < tk_l; j++){
      if(baumwelch[i][j]!=0.0){
	double tmpQ = trans[i][j]*baumwelch[i][j];
	Q += tmpQ;
	if(1&&std::isnan(tmpQ)){
	  fprintf(stderr,"\t->[qFunction_inner2]i:%d j:%d tmpQ:%e trans::%e baumwel:%e\n",i,j,tmpQ,trans[i][j],baumwelch[i][j]);
	  exit(0);
	}
      }

    }
  }
   if(std::isnan(Q)){
     fprintf(stderr,"Q is nan will exit\n");
     exit(0);
   }
   return Q;
  
}

void fastPSMC::calculate_FW_BW_Probs(double *tk,int tk_l,double *epsize,double rho){
    //we first set the initial fwprobs to stationary distribution
  for(int i=0;i<tk_l;i++)
      fw[i][0] = stationary[i];

    //we now loop over windows.
    //v=0 is above and is the initial distribution, we therefore plug in at v+1
    for(int v=0;v<windows.size();v++){
      ComputeRs(v,fw);//<-prepare R1,R2

      fw[0][v+1] = addProtect3(lprod(fw[0][v],P[1][0]) , lprod(R1[0],P[3][0]) , lprod(fw[0][v],P[4][0]))+emis[0][v+1] ;
      for (unsigned i = 1; i < tk_l; i++)
	fw[i][v+1]= addProtect4(lprod(fw[i][v],P[1][i]) , lprod(R2[i-1],P[2][i]) , lprod(R1[i],P[3][i]) , lprod(fw[i][v],P[4][i]))+emis[i][v+1];
    }

    double tmp[tk_l];
    for(int i=0;i<tk_l;i++){
      tmp[i] = fw[i][windows.size()];
    }
    pix = addProtectN(tmp,tk_l);

    assert(!std::isnan(pix));
    fwllh = pix;

    //now do backward algorithm
    //initialize by stationary
    for(int i=0;i<tk_l;i++)
      bw[i][windows.size()] = lprod(stationary[i],emis[i][windows.size()]);
     

    //we plug in values at v-1, therefore we break at v==1
    for(int v=windows.size();v>0;v--){
      ComputeRs(v,bw);//<-prepare R1,R2

      bw[0][v-1] = addProtect3(lprod(bw[0][v],P[1][0]) , lprod( R1[0],P[3][0]) , lprod(bw[0][v],P[4][0]))+emis[0][v-1];
      bw[0][v] -= lprod(stationary[0],emis[0][v]);
      for (unsigned i = 1; i < tk_l; i++){
	bw[i][v-1] = addProtect4(lprod(bw[i][v],P[1][i]),lprod(R2[i-1],P[2][i]),lprod(R1[i],P[3][i]),lprod(bw[i][v],P[4][i]))+emis[i][v-1];
	bw[i][v] -= lprod(stationary[i],emis[i][v]);
      }
    }

    for(int i=0;i<tk_l;i++)
      bw[i][0] -= stationary[i];
     

    for(int i=0;i<tk_l;i++)
      tmp[i] = bw[i][1]+stationary[i]+emis[i][1];
    double tmptmp= addProtectN(tmp,tk_l);
    assert(!std::isnan(tmptmp));
    bwllh = tmptmp;
 
}

void fastPSMC::allocate(int tk_l_arg){
  int numWindows = windows.size();
  tk_l = tk_l_arg;
  if(index==0)
    stationary  = new double[tk_l];
  R1 = new double[tk_l];
  R2 = new double[tk_l];
  fw = new double *[tk_l];
  bw = new double *[tk_l];
  //pp = new double *[tk_l];
  emis = new double *[tk_l];
  baumwelch = new double *[tk_l+1];
  for(int i=0;i<tk_l;i++){
    emis[i] = new double[numWindows+1];
    baumwelch[i] = new double[tk_l];
    fw[i] = new double[numWindows+1];
    bw[i] = new double[numWindows+1];
  }
  baumwelch[tk_l] = new double[tk_l];
  for(int i=0;i<tk_l+1;i++)
    for(int j=0;j<tk_l;j++)
      baumwelch[i][j] = -777;
  if(index==0)
  P = new double *[8];
  PP= new double *[8];
  nP = new double*[8];
  for(int i=0;i<8;i++){
    if(index==0)
    P[i] = new double[tk_l];
    PP[i]= new double[tk_l];
    nP[i]= new double[tk_l];
  }
  for(int i=0;i<tk_l;i++)
    PP[0][i] = -666;
  workspace = new double[windows.size()];
  if(DOTRANS&&index==0){
    //    fprintf(stderr,"allocating\n");
    trans = new double *[tk_l];
    for(int i=0;i<tk_l;i++){
      trans[i] = new double[tk_l];
      for(int j=0;j<tk_l;j++)
	trans[i][j] = -888;//placeholder, to spot if something shouldnt be happening;
    }
  }

}
/*
  Function will set the indices for the windows
  first index and last index INCLUSIVE
 */
void fastPSMC::setWindows(int *pos ,int last,int block){
  int beginIndex =0;
  int endIndex=0;
  int beginPos = 0;
  int endPos = beginPos+block-1;
  
  while(1) {
    //fprintf(stdout,"Beginpos:%d EndPos:%d\n",beginPos,endPos);
    wins w;
    if(endPos>pos[last-1])
      break;
    
    while(pos[beginIndex]<beginPos)
      beginIndex++;
    endIndex=beginIndex;
    if(pos[endIndex]<endPos){
      while(pos[endIndex]<endPos){
	//	fprintf(stdout,"runs?\n");
	endIndex++;
      }
      endIndex--;
    }
#if 0
    fprintf(stdout,"\t-> winsize:%d bp:%d,ep:%d bi:%d ei:%d ei-bi:%d\n",block,beginPos,endPos,beginIndex,endIndex,endIndex-beginIndex);
#endif
    w.from = beginIndex;
    w.to = endIndex;
    windows.push_back(w);
    beginPos+=block;
    endPos+=block;
  }

}


/* 
  Calculate emission probabilityes
  tk array of length tk_l
  lambda array effective population sizes
  both has length tk_l
  
  emission probablities will be put in the **emis
  
  stationary(i) = exp(-sum_{j=0}^{i-1}{tau_j/lambda_j}*P2[i])
 */

void calculate_emissions(double *tk,int tk_l,mygltype *gls,std::vector<wins> &windows,double theta,double **emis,double *epsize){

  //initialize the first:
  for(int j=0;j<tk_l;j++)
    emis[j][0] = log(0);
 
  //  double tmp[windows.size()];
  double nontmpdir[tk_l];
  double expectCoalT[tk_l];
  int emis_approx = 0;
  if (emis_approx == 0)
    ComputeP1(tk,tk_l,nontmpdir,epsize,theta);
  else if (emis_approx == 1){
    ComputeExpectedCoalTime(tk, tk_l, expectCoalT, epsize);
    for (int i = 0; i < tk_l; i++)
      nontmpdir[i] = -2*theta*expectCoalT[i];
  }
  else
    assert(0!=0);

 for(int v=0;v<windows.size();v++){//for each window
   assert(windows[v].from<=windows[v].to);
   //   fprintf(stderr,"v:%d from:%d to:%d\n",v,windows[v].from,windows[v].to);
    for(int j=0;j<tk_l;j++){//for each interval/state
      emis[j][v+1] = 0;

      double inner = exp(nontmpdir[j]);///exp(-2.0*tk[j]*theta); // this part relates to issue #1
#if 0
      //      double inner;
      if (j< tk_l-1)
	inner = exp(-2.0*(tk[j]+tk[j+1])/2.0*theta); // this part relates to issue #1
      else
	inner = exp(-2.0*(tk[j]+tk[j])*theta); // this part relates to issue #1
      //      fprintf(stderr,"\t\t%d from:%d to:%d inner:%f\n",j,windows[v].from,windows[v].to,inner);
#endif
      for(int i=windows[v].from;i<=windows[v].to;i++) {//for all elements in window
	double igl[2]={log(0),log(0)};
	double val = gls[i];
	if(!std::isinf(gls[i])){
	  if(sizeof(mygltype)==1){
	    // fprintf(stderr,"gl: %d ",(int)val);
	    val = log(pow(10.0,val/-100.0));
	    //fprintf(stderr,"gl: %f \n",val);
	    
	  }
	  if(gls[i]<0){
	    igl[0]=0;
	    igl[1]=val;
	  }else{
	    igl[0]=-val;
	    igl[1]=0;
	  }
	  //	  fprintf(stdout,"dims\t%f\t%f\n",igl[0],igl[1]);
	}

	if(igl[0]!=igl[1])
	  emis[j][v+1] += log((exp(igl[0])/4.0) *inner + (exp(igl[1])/6.0)*(1.0-inner));//<- check

	if(std::isinf(emis[j][v+1])){
	  fprintf(stderr,"\t-> Huge bug in code contact developer. Emissions evaluates to zero\n");
	  /*
	    This will corrupt the computation of the backward probability.
	  */
	  exit(0);
	}
      }
    }
  }
}

void ComputePii(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary,double **emis,double *workspace){

  ComputeP11(numWind,tk_l,P[1],PP[1],fw,bw,workspace,emis);
  ComputeP22(numWind,tk_l,P,PP[2],fw,bw,emis);
  ComputeP33(numWind,tk_l,P[3],PP[3],fw,bw,emis);
  ComputeP44(numWind,tk_l,P[4],PP[4],fw,bw,workspace,emis);
  ComputeP55(numWind,tk_l,P,PP[5],fw,bw,stationary,emis);
  ComputeP66(numWind,tk_l,P,PP[6],fw,bw,stationary,emis);
  ComputeP77(numWind,tk_l,P,PP[7],fw,bw,stationary,emis);
  
  
  for(int p=1;0&&p<8;p++)//CHECK IF THIS SHOULD BE RENAABLED
    for(int i=0;i<tk_l;i++)
      assert(exp(PP[p][i])>=0&&exp(PP[p][i])<=1);


}


//results is in baumwelch matrix the content is not in log but normal space
void ComputeBaumWelch(unsigned numWind,int tk_l,double **fw,double **bw,double **emis,double **trans,double **baumwelch,double pix){
  for(int i=0;i<tk_l;i++){
    for(int j=0;j<tk_l;j++){
      double tmp = log(0);
      for(int w=1;w<numWind;w++){
	tmp = addProtect2(tmp,fw[i][w]+trans[i][j]+emis[j][w+1]+bw[j][w+1]);
	if(0&&w>20)
	  exit(0);
      }
      baumwelch[i][j] = exp(tmp-pix);
    }
  }
  for(int i=0; i < tk_l; i++)
    baumwelch[tk_l][i] = exp(fw[i][1]+bw[i][1]-pix);

}

//should only be run once, since calculates 
void fastPSMC::make_hmm_pre(double *tk,int tk_l,double *epsize,double theta,double rho){
  assert(trans);
  assert(index==0);
  ComputeGlobalProbabilities(tk,tk_l,P,epsize,rho);//only the P* ones
  calculate_stationary(tk_l,stationary,P);
  if(DOTRANS){
    for(int i=0;i<tk_l;i++)
      for(int j=0;j<tk_l;j++){
	//	fprintf(stderr,"trans:%p\n",trans);
	trans[i][j] = calc_trans(i,j,P);
      }
  }

}


double fastPSMC::make_hmm(double *tk,int tk_l,double *epsize,double theta,double rho){
  //prepare probs

  calculate_emissions(tk,tk_l,gls,windows,theta,emis,epsize);
  // print_emission("emis.txt");
  calculate_FW_BW_Probs(tk,tk_l,epsize,rho);


  if(doQuadratic==0)
    ComputePii(windows.size(),tk_l,P,PP,fw,bw,stationary,emis,workspace);
  else
    ComputeBaumWelch(windows.size(),tk_l,fw,bw,emis,trans,baumwelch,pix);

  //no need to recompute P. But we fix this later;
  if(doQuadratic==0)
    qval=qFunction_inner(tk,tk_l,epsize,rho,pix,windows.size(),P,PP);
  else
    qval=qFunction_inner2(tk,tk_l,epsize,rho,pix,windows.size(),P,baumwelch,trans);
  return qval;
 }

fastPSMC::~fastPSMC(){
  delete [] gls;
  delete [] R1;
  delete [] R2;
  for(int i=0;i<tk_l;i++){
    delete [] emis[i];
    delete [] fw[i];
    delete [] bw[i];
    delete [] baumwelch[i];
  }
  delete [] baumwelch[tk_l];
  for(int i=0;i<8;i++){
    if(index==0)
      delete [] P[i];
    delete [] PP[i];
    if(index==0)
    delete [] nP[i];
  }
  delete [] emis;
  delete [] fw;
  delete [] bw;
  delete [] baumwelch;
  if(index==0)
  delete [] P;
  delete [] PP;
  if(index==0)
  delete [] nP;
  if(DOTRANS&&index==0){
    for(int i=0;i<tk_l;i++)
      delete [] trans[i];
    delete [] trans;
  }
  delete [] workspace;
  if(index==0)
  delete [] stationary;
}
