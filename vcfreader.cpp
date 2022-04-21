// #define TEST

#include "vcfreader.h"
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <iostream>
#include <cmath>
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
        printf("%s\t%ld\n", ctg[i].key, ctg[i].val->info[0]);
        if (strcmp(ctg[i].key,chr)==0){
            length = ctg[i].val->info[0];
        }



    }

    

    tmpgls = (int*)calloc(2*length,sizeof(int));
    int k = 0;
    int pos = -1;
    
    while(bcf_read(fp, hdr, rec) == 0){
        if (rec->rid != rid)continue;
        rec ->
        if (pos == -1)pos = k;
        bcf_unpack(rec, BCF_UN_ALL);
        npl = bcf_get_format_int32(hdr, rec, "PL", &pl, &npl_arr);
        //homozygotic sites are on 1,1+1,1+1+2,1+1+2+3

        //TODO:  !!! CONSIDER <*> IN CODE !!! 
        int u = 0;
        int iter = 1;
//        fprintf(stderr,"@%d\n",rec->n_allele);
        for (int s = 0; s < rec->n_allele;s++){
            tmpgls[2*k] += pl[u];//Homozygotic sum

            for(int v = u+1; v < u+iter+1; v++){
                tmpgls[2*k+1] += pl[u+v];//Heterozygotic sum
            }    
            u += iter+1;
            iter += 1;
           
        }
//        fprintf(stderr,"@tmpgls: %d,%d\n",tmpgls[2*k],tmpgls[2*k+1]);
        k+=1;

    }

    
    ret.pos = new int[1];
    ret.len = k+1; 
    ret.gls = new double[k];
    ret.pos[0] = pos;

    for (int i = 0; i<k;i++) {
        ret.gls[i] = log(0.0);// = -infinity
        if (tmpgls[2 * i] != tmpgls[2 * i + 1]) {
            fprintf(stderr,"@tmps: %d,%d\n",tmpgls[2 * i], tmpgls[2 * i + 1]);
            double mmax = std::max(tmpgls[2 * i]/4, tmpgls[2 * i + 1]/6);
            double val =
                    std::min(tmpgls[2 * i]/4, tmpgls[2 * i + 1]/6) - mmax; //??? Why fill it min - max ANS:PREVENT UNDERFLOW

            ret.gls[i] = val;
            if (tmpgls[2 * i]/4 < tmpgls[2 * i + 1]/6)
                ret.gls[i] = -ret.gls[i];//????Why reverse

//        }

            fprintf(stderr, "@ret.gls[%d] = %f\n",i, ret.gls[i]);
        }
    }
    ret.lastp = k;
    return ret;

}

// int main(){
//     rawdata check;
//     check = readvcf("NA12878.chr22.vcf.gz","22");
// }