#include <cstdio>
#include <htslib/faidx.h>
#include "header.h"

typedef struct{
  char *fastaname;
  faidx_t *fai;//contains the faidx structure
  char *seqs;//contains the reference for the current chr;
  int curChr;//the exact chromosome name for the seqs above
  int chrLen;//length of chromosome
}perFasta;


//checks that newer is newer than older
int isNewer(const char *newer,const char *older){
   if (strstr(older, "ftp://") == older || strstr(older, "http://") == older)
     return 0;
  //  fprintf(stderr,"newer:%s older:%s\n",newer,older);
  // return 0;
  struct stat one;
  struct stat two;
  stat(newer, &one );
  stat(older, &two );
  
  return one.st_mtime>=two.st_mtime;
}

void perFasta_destroy(perFasta *pf){
  if(pf->fastaname)
    free(pf->fastaname);
  if(pf->fai)
    fai_destroy(pf->fai);
  if(pf->seqs)
    free(pf->seqs);
  delete pf;
}

//this will initialize our data
perFasta *perFasta_init(const char *fname){

  fprintf(stderr,"\t-> Reading fasta: %s\n",fname);
  if(fexists(fname)==0){
    fprintf(stderr,"\t-> fastafile: \'%s\' doesn't exists, will exit\n",fname);
    exit(0);
  }
  //check that fa hasn't been updated
  char *strtsk=NULL;
  strtsk = (char*)calloc(strlen(fname) + 5, 1);
  sprintf(strtsk, "%s.fai", fname);
  if(isNewer(fname,strtsk)){
    fprintf(stderr,"\t-> fai index file: \'%s\' looks older than corresponding fastafile: \'%s\'.\n\t-> Please reindex fasta file\n",strtsk,fname);
    exit(0);
  }
  free(strtsk);
    
  perFasta *r= new perFasta;
  r->fastaname=NULL;
  r->fai = NULL;
  r->seqs = NULL;
  r->curChr = -1;

  if(NULL==(r->fai = fai_load(fname))){
    fprintf(stderr,"[%s:%s] error reading fai file:%s\n",__FILE__,__FUNCTION__,fname);
    exit(0);
  }
  r->fastaname=strdup(fname);
  return r;
}



char *loadChr(perFasta *f, char*chrName,int chrId){
  //  fprintf(stderr,"[%s] \t->loading chr:%s from faidx=%p curchr=%d\n",__FUNCTION__,chrName,f,f->curChr);
  free(f->seqs);
  f->seqs=NULL;
  //fprintf(stderr,"f->curChr=%d chrId=%d\n",f->curChr,chrId);
  f->curChr = chrId;
  f->seqs = faidx_fetch_seq(f->fai, chrName, 0, 0x7fffffff, &f->chrLen);
  if(f->seqs==NULL){
    fprintf(stderr,"\n[%s] Error loading fasta info from chr:\'%s\' \n",__FUNCTION__,chrName);
    f->chrLen=0;
  }
  //  fprintf(stderr,"[%s] done\n",__FUNCTION__);

  return f->seqs;
}

