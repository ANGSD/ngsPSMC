getVec <- function(x){
file<-x
d<-as.numeric(unlist(strsplit(system(paste0("grep PA ",file," | tail -n1|cut -f2|cut -f1-4 --complement -d\" \""|tr -d \"/\"),intern=T)," ")))
d
}

psmc <- getVec("../test25/angsd.psmc")
ngspsmc1.fa <- getVec("../dongdong1.stdout")
ngspsmc2.fa <- getVec("../dongdong2.stdout")
ngspsmc3.gl <- getVec("../test25/delme3.stdout")
ngspsmc4.gl <- getVec("../test25/delme4.stdout")

pdf("test.pdf")
plot(1:28,psmc,ylim=c(-1,10),type='b',lwd=2,col=1)
lines(1:28,ngspsmc1.fa,type='b',lwd=2,col=2)
lines(1:28,ngspsmc2.fa,type='b',lwd=2,col=3)
lines(1:28,ngspsmc3.gl,type='b',lwd=2,col=4)
lines(1:28,ngspsmc4.gl,type='b',lwd=2,col=5)

legend("topright",c("PSMC","ngsPSMC.fa (init from PSMC)","ngsPSMC.fa init=1","ngsPSMC.gl init=1 (div 4/10)","ngsPSMC.gl init=1 (no div 4/10)"),fill=1:5)
dev.off()


psmc <- getVec("../test25/angsd.psmc")
ngspsmc1.gl <- getVec("../test25/delme3.stdout")	
sub.ngspsmc1.gl <- getVec("../test25/sub.delme3.stdout")	
ngspsmc2.gl <- getVec("../test25/delme4.stdout")	
sub.ngspsmc2.gl <- getVec("../test25/sub.delme4.stdout")	
dongdong1  <- getVec("../test25/dongdong1.stdout")	
dongdong2  <- getVec("../test25/dongdong2.stdout")	

ngspsmc1.fa <- getVec("../test25/dongdong1.stdout")
ngspsmc2.fa <- getVec("../test25/dongdong2.stdout")
ngspsmc3.gl <- getVec("../test25/delme3.stdout")
ngspsmc4.gl <- getVec("../test25/delme4.stdout")

pdf("test.pdf")
plot(1:28,psmc,ylim=c(-1,10),type='b',lwd=2,col=1)
lines(1:28,dongdong1,type='b',lwd=1,col=2)
lines(1:28,dongdong2,type='b',lwd=1,col=3)
lines(1:28,ngspsmc3.gl,type='b',lwd=2,col=4)
lines(1:28,ngspsmc4.gl,type='b',lwd=2,col=5)

legend("topright",c("PSMC","ngsPSMC.fa (init from PSMC)","ngsPSMC.fa init=1","ngsPSMC.gl init=1 (div 4/10)","ngsPSMC.gl init=1 (no div 4/10)"),fill=1:5)
dev.off()


pdf("fromfasta.pdf")
plot(1:28,psmc,ylim=c(-1,10),type='b',lwd=2,col=1,main='fromfasta')
lines(1:28,dongdong1,type='b',lwd=1,col=2)
lines(1:28,dongdong2,type='b',lwd=1,col=3)
legend("topright",c("PSMC","ngsPSMC.fasta (optim from PSMC)","ngsPSMC.fasta init=1"),fill=1:3)
dev.off()

pdf("subset.pdf",width=14)
par(mfrow=c(1,2))

plot(1:28,psmc,ylim=c(-1,10),type='b',lwd=3,col=1,main='full genome, half genome comparison')
lines(1:28,ngspsmc1.gl,type='b',lwd=2,col=2)
lines(1:28,sub.ngspsmc1.gl,type='b',lwd=1,col=3)
legend("topright",c("PSMC","ngsPSMC (div 4/10)","0.5 ngsPSMC (div 4/10)"),fill=1:3)
plot(1:28,psmc,ylim=c(-1,10),type='b',lwd=3,col=1,main='full genome, half genome comparison')
lines(1:28,ngspsmc2.gl,type='b',lwd=2,col=4)
lines(1:28,sub.ngspsmc2.gl,type='b',lwd=1,col=5)
legend("topright",c("PSMC","ngsPSMC (nodiv 4/10)","0.5 ngsPSMC (nodiv 4/10)"),fill=1:3)

dev.off()

if(FALSE){
    pdf("plots.pdf")
    source("read.psmc.R")
    fl <- list.files(".",pattern=".psmc$")
    plot.psmc(read.psmc(fl[1]),lwd=2,col=1,ylim=c(0,5),main="Only psmc")
    plot.psmc(read.psmc(fl[2]),lwd=2,add=T,col=2)
    plot.psmc(read.psmc(fl[3]),lwd=2,add=T,col=3)
    legend("top",c("CEU","YRI","YRI"),fill=1:3)

    plot.psmc(read.psmc("NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.fa.psmcfa.psmc"),lwd=2,col=1,ylim=c(0,5),main="CEU NA12878 ")
    plot.psmc(read.psmc("NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.psmc.idx.1.stdout"),lwd=2,col=2,add=TRUE)
    plot.psmc(read.psmc("NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.psmc.idx.2.stdout"),lwd=2,col=3,add=TRUE)
    legend("topleft",c("PSMC","ngs init 1","ngs init -1"),fill=1:3)

     plot.psmc(read.psmc("NA18507.mapped.ILLUMINA.bwa.YRI.low_coverage.20120522.fa.psmcfa.psmc"),lwd=2,col=1,ylim=c(0,10),main="YRI NA18507 ")
    plot.psmc(read.psmc("NA18507.mapped.ILLUMINA.bwa.YRI.low_coverage.20120522.psmc.idx.1.stdout"),lwd=2,col=2,add=TRUE)
    plot.psmc(read.psmc("NA18507.mapped.ILLUMINA.bwa.YRI.low_coverage.20120522.psmc.idx.2.stdout"),lwd=2,col=3,add=TRUE)
    legend("topleft",c("PSMC","ngs init 1","ngs init -1"),fill=1:3)

    plot.psmc(read.psmc("NA19239.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.fa.psmcfa.psmc"),lwd=2,col=1,ylim=c(0,5),main="YRI NA19239 ")
    plot.psmc(read.psmc("NA19239.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.psmc.idx.1.stdout"),lwd=2,col=2,add=TRUE)
    plot.psmc(read.psmc("NA19239.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.psmc.idx.2.stdout"),lwd=2,col=3,add=TRUE)
    legend("topleft",c("PSMC","ngs init 1","ngs init -1"),fill=1:3)

    fl <- c(list.files(".",pattern=".psmc$"),list.files(".",pattern=".stdout$"))
    for(i in fl)
        plot_chain(read.psmc(i),main=i)
dev.off()
    

}
