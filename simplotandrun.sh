DIR=test3/
mkdir -p ${DIR}
HOT=../psmc_project/foreign/msHOT-lite/msHOT-lite
TOFA=../psmc_project/psmc/utils/ms2psmcfa.pl
PSMC=../psmc_project/psmc/psmc

${HOT} 2 20 -t 30000 -r 6000 30000000 -eN 0.01 0.1 -eN 0.06 1 -eN 0.2 0.5 -eN 1 1 -eN 2 2 -l >${DIR}/msraw
${TOFA} <${DIR}/msraw|gzip -c  >${DIR}/ms.lh3.fa.gz
${PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o ${DIR}/ms.lh3.fa.gz.psmc ${DIR}/ms.lh3.fa.gz 
../angsd/misc/msHOT2glf -in ${DIR}/msraw -out ${DIR}/angsd -psmcfa 1 -psmc2 1 -do_seq_glf 1
../angsd/misc/plot_psmc.py test3/msraw test3/ms.lh3.fa.










