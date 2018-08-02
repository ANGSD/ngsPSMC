getVec <- function(x){
file<-x
d<-as.numeric(unlist(strsplit(system(paste0("grep PA ",file," | tail -n1|cut -f2|cut -f1-4 --complement -d\" \""),intern=T)," ")))
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