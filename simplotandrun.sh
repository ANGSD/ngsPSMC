DIR=test4/
mkdir -p ${DIR}
HOT=../foreign/msHOT-lite/msHOT-lite
TOFA=../psmc/utils/ms2psmcfa.pl
PSMC=../psmc/psmc

${HOT} 2 20 -t 30000 -r 6000 30000000 -eN 0.01 0.1 -eN 0.06 1 -eN 0.2 0.5 -eN 1 1 -eN 2 2 -l >${DIR}/msraw
${TOFA} <${DIR}/msraw|gzip -c  >${DIR}/ms.lh3.fa.gz
##cleanup name in firstline
../angsd/misc/msHOT2glf -in ${DIR}/msraw -out ${DIR}/angsd -psmcfa 1 -psmc2 1 -do_seq_glf 1

${PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o ${DIR}/ms.lh3.fa.gz.psmc ${DIR}/ms.lh3.fa.gz 
${PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o ${DIR}/angsd.psmc ${DIR}/angsd.fa.gz
bgzip -r test4/angsd.fa.gz
samtools faidx test4/angsd.fa.gz
../angsd/misc/plot_psmc.py ${DIR}/msraw ${DIR}/ms.lh3.fa.gz.psmc ${DIR}/angsd.psmc

 ./ngsPSMC mshot.fa.gz -infile newdelme2 -r 1







DIR=test5/
mkdir -p ${DIR}
HOT=../foreign/msHOT-lite/msHOT-lite
TOFA=../psmc/utils/ms2psmcfa.pl
PSMC=../psmc/psmc

${HOT} 2 1 -t 300000 -r 60000 300000000 -eN 0.01 0.1 -eN 0.06 1 -eN 0.2 0.5 -eN 1 1 -eN 2 2 -l >${DIR}/msraw
${TOFA} <${DIR}/msraw|gzip -c  >${DIR}/ms.lh3.fa.gz
##cleanup name in firstline
../angsd/misc/msHOT2glf -in ${DIR}/msraw -out ${DIR}/angsd -psmcfa 1 -psmc2 1 -do_seq_glf 1 &

${PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o ${DIR}/ms.lh3.fa.gz.psmc ${DIR}/ms.lh3.fa.gz 
${PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o ${DIR}/angsd.psmc ${DIR}/angsd.fa.gz
bgzip -r test4/angsd.fa.gz
samtools faidx test4/angsd.fa.gz
../angsd/misc/plot_psmc.py ${DIR}/msraw ${DIR}/ms.lh3.fa.gz.psmc ${DIR}/angsd.psmc

 ./ngsPSMC mshot.fa.gz -infile newdelme2 -r 1





double ary[1000000];//<- stack
double *ary = new double[10];;<- heap
double *ary = malloc(sizeof(double)*10)

