// #define TEST

#include "vcfreader.h"
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <iostream>
/*
For name of vcf file and chromosome name(???ID) returns struct rawdata of gls scores

PARAMETERS
@ inf: char*
file name string
@ chr: char*
chrom name

RETURNS 
@ret: rawdata
structure




*/
rawdata readvcf(char* inf,char* chr){//infile

    int *pl = NULL;
    int npl_arr = 0;
    int npl = 0;

    int ngt_arr = 0;
    int ngt     = 0;
    int *gt     = NULL;
    int nad_arr = 0;
    int nad     = 0;
    int *ad     = NULL;

    rawdata ret;
    htsFile *fp = bcf_open(inf, "r");
    bcf1_t *rec = bcf_init();
    int *tmpgls;
    
    if(fp == NULL) {
        throw std::runtime_error("Unable to open file.");
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);

    if(hdr == NULL) {
        throw std::runtime_error("Unable to read header.");
    }
    fprintf(stderr,"@%s",chr);
    int rid = bcf_hdr_name2id(hdr,chr);
    
    bcf_idpair_t *ctg =hdr->id[BCF_DT_CTG];
    printf("file %s is successfully read\n",inf);
    int length;
    printf("Information from header:\n");
    fprintf(stderr,"@file: %s\n",inf);
    fprintf(stderr,"@chr: %s\n",chr);
    
    for (int i = 0; i < hdr->n[BCF_DT_CTG]; ++i){
        // printf("%s\t%ld\n", ctg[i].key, ctg[i].val->info[0]);
        if (strcmp(ctg[i].key,chr)==0){
            length = ctg[i].val->info[0];
        }



    }

    

    tmpgls = (int*)calloc(2*length,sizeof(int));
    int k = 0;
    int pos = -1;
    
    while(bcf_read(fp, hdr, rec) == 0){
        // fprintf(stderr,"@%d@%d\n",length,2*k);
        if (rec->rid!=rid)continue;
        if (pos == -1)pos = k;
        bcf_unpack(rec, BCF_UN_ALL);
        npl = bcf_get_format_int32(hdr, rec, "PL", &pl, &npl_arr);
        //homozygotic sites are on 1,1+1,1+1+2,1+1+2+3
        //TODO:  !!! CONSIDER <*> IN CODE !!! 
        int u = 0;
        for (int s = 0; s < rec->n_allele;s++){
            tmpgls[2*k] += pl[u];//Homozygotic sum
        
            for(int v = 0; v < s; v++){
                tmpgls[2*k+1] += pl[u+v];//Heterozygotic sum
            }    
            u += s+1;
           
        }
        k+=1;   
    }

    
    ret.pos = new int[1];
    ret.len = k+1; 
    ret.gls = new double[2*k+2];
    ret.pos[0] = pos;

    for (int i = 0; i<k+1;i++)
    {
        ret.gls[2*i] = (mygltype)tmpgls[i]/4;
        ret.gls[2*i+1] = (mygltype)tmpgls[i]/6;
        fprintf(stderr,"@gls[%d] = %lf\n@gls[%d] = %lf\n",2*i,ret.gls[2*i],2*i+1,ret.gls[2*i+1]);
    }
    return ret;

}

// int main(){
//     rawdata check;
//     check = readvcf("out4.vcf.gz","chrM");
// }