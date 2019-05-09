

read.psmc1 <- function(x,rd=-1){
    d<-readLines(x)
    if(length(d)<10)
        stop("file looks empty")
    d <- grep("CC",val=T,invert=T,d)
    rds <- grep("RD",d)
    rds <- c(rds,length(d))

    ofs<-rds[1:2]
    for(i in 2:(length(rds)-1))
        ofs<-rbind(ofs,rds[i:(i+1)])
    if(rd==-1){
        rd=nrow(ofs)
    }
    if(rd>nrow(ofs))
        stop("rd is out of bounds")
    sd <- d[ofs[rd,1]:ofs[rd,2]]
    
  
    TR <-as.numeric(unlist(strsplit(grep("TR",sd,val=T),"\t"))[-1])
    LK <-as.numeric(unlist(strsplit(grep("LK",sd,val=T),"\t"))[-1])
    MT <- as.numeric(unlist(strsplit(grep("MT",sd,val=T),"\t"))[-1])
    nrs <- length(grep("RS",sd))
    RS <- matrix(as.numeric(matrix(unlist(strsplit(grep("RS",sd,val=T),"\t")),nrow=nrs,byrow=T)[,-1]),nrow=nrs,byrow=F)
    RS[,2]<-log(RS[,2])
    PAT <-unlist(strsplit(unlist(strsplit(grep("PA",sd,val=T),"\t"))," "))[[2]] 
    PA <-  as.numeric(unlist(strsplit(unlist(strsplit(grep("PA",sd,val=T),"\t"))," "))[-c(1,2)])
    ret <- list(TR=TR,MT=MT,RS=RS,PAT=PAT,PA=PA,LK=LK)
    return(ret)
    
}

read.psmc <- function(x){
    d <- grep("CC",val=T,invert=T,readLines(x))
    nrd <- length(grep("RD",d))
    cat("\t-> Number of RD from file: ",x," is ",nrd,"\n")
    ret <- lapply(1:nrd,function(rd) read.psmc1(x,rd))
    return(ret)
}


if(FALSE){
    source("read.psmc.R")
    x<-"data_std_nomigr2/ms2g2.psmc"
    ret <- read.psmc(x)

}

if(FALSE){
    pdf("forvlad.pdf")
    d0 <- read.psmc1("test27/ngs.ms.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs.ms.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="PSMC optimized pars,win=1",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])) )
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)
    d0 <- read.psmc1("test27/ngs.ms2.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs.ms2.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="PSMC optimized pars,win=100",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)

    d0 <- read.psmc1("test27/ngs2.ms.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs2.ms.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="MS  pars,win=1",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)
    d0 <- read.psmc1("test27/ngs2.ms2.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs2.ms2.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="MS  pars,win=100 version1)",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="MS  pars,win=100 version2)",ylim=c(0,2))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)

    dev.off()
}


if(FALSE){
    pdf("forvlad3.pdf")
    d0 <- read.psmc1("test27/ngs5.ms.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs5.ms.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="MS GL 5+5+5+5,win=100",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])) )
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)

    d0 <- read.psmc1("test27/ngs6.ms.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs6.ms.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="MS GL 10+10+10+10,win=100",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)

    d0 <- read.psmc1("test27/ngs7.ms.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs7.ms.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="MS GL 10+10+10+10 win=1",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)
    d0 <- read.psmc1("test27/ngs8.ms.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs8.ms.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="MS GL 5+5+5+5,win=1 )",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)
    dev.off()
}


if(FALSE){
	   pdf("forvlad4.pdf")
    d0 <- read.psmc1("test27/ngs.ms.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs.ms.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="fasta input,psmc pars,win=1",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])) )
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)

    d0 <- read.psmc1("test27/ngs.ms2.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs.ms2.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="fasta input psmc pars, win=100",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)

    d0 <- read.psmc1("test27/ngs2.gl.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs2.gl.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="gl input psmc pars win=1",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)

    d0 <- read.psmc1("test27/ngs3.ms.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs3.ms.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="gl input psmc pars win=100 ",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)

    d0 <- read.psmc1("test27/ngs4.fa.ms.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs4.fa.ms.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="fa input ms pars win=1 (4x5) ",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)


       d0 <- read.psmc1("test27/ngs5.fa.ms.stdout",rd=0)
    d1 <- read.psmc1("test27/ngs5.fa.ms.stdout",rd=1)
    plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,main="fa input ms pars win=100 (4x5) ",ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])))
    lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
    legend("topright",paste0(c("init.llh","optim1.llh"),"=",c(d0$LK,d1$LK)),fill=1:2)


    


    dev.off()	


}

plot_chain<-function(x,ylim,type,...){
    colfunc<-colorRampPalette(c("black","red"))
    colo=colfunc(length(x))
    if(missing(ylim))
        ylim=c(0,max(sapply(x,function(x) x$RS[,3])))
    if(missing(type))
        type='s'
    plot(x[[1]]$RS[,2],x[[1]]$RS[,3],type=type,ylim=ylim,col=colo[1],lwd=6,...)
    for (i in 2:length(x)){
        d<-x[[i]]
        lines(d$RS[,2],d$RS[,3],col=colo[i],lwd=2,type=type)
    }
    legend("topright",paste0(c("init.llh","optim.llh"),"=",c(x[[1]]$LK,tail(x,1)[[1]]$LK)),fill=1:2)
}

plot.psmc<-function(x,ylim,type,col,lwd,main,add=FALSE,...){
    if(missing(ylim))
        ylim=c(0,max(sapply(x,function(x) x$RS[,3])))
    if(missing(type))
        type='s'
    if(missing(col))
        col='black'
    if(missing(lwd))
        lwd=2
    if(missing(main))
        main=NULL
    x <- tail(x,1)
    if(add==FALSE)
        plot(x[[1]]$RS[,2],x[[1]]$RS[,3],type=type,ylim=ylim,col=col,lwd=lwd,main=main,...)
    else
        lines(x[[1]]$RS[,2],x[[1]]$RS[,3],type=type,ylim=ylim,col=col,lwd=lwd,...)
   # legend("topright",paste0(c("optim.llh"),"=",c(x[[1]]$LK),fill=1))
}

plot_ngs <-function(x,...){
	 d0<-x[[1]]
	 d1<-x[[2]]
  plot(d0$RS[,2],d0$RS[,3],type='l',lwd=6,col=1,ylim=c(0,max(rbind(d0$RS,d1$RS)[,3])),...)
 lines(d1$RS[,2],d1$RS[,3],type='l',lwd=3,col=2)
}


if(FALSE){
    #   pdf("forvlad4.pdf")
    lh3.w100.fa <- read.psmc("test28/lh3.w100.fa.gz.psmc")
    tsk.w100.fa <- read.psmc("test28/tsk.w100.fa.gz.psmc")	
    lh3.w1.fa <- read.psmc("test28/lh3.w1.fa.gz.psmc")
     tsk.w1.fa <- read.psmc("test28/tsk.w1.fa.gz.psmc")	

    

  plot_chain(lh3.w100.fa)  
  plot_chain(tsk.w100.fa)				
   plot_chain(lh3.w1.fa)
    plot_chain(lh3.w1.fa)
    
#    plot_ngs(ngs.tsk.w100.fa)
   
    plot(lh3.w1.fa[[26]]$RS[,2],lh3.w1.fa[[26]]$RS[,3],type='l',lwd=3,col=1,main='fasta (win1 vs win100)(lh3)')
    lines(lh3.w100.fa[[26]]$RS[,2],lh3.w100.fa[[26]]$RS[,3],lwd=1,col=2)	

    plot(tsk.w1.fa[[26]]$RS[,2],tsk.w1.fa[[26]]$RS[,3],type='l',lwd=4,col=1,main='fasta (win1 vs win100)(tsk)')
    lines(tsk.w100.fa[[26]]$RS[,2],tsk.w100.fa[[26]]$RS[,3],lwd=2,col=2)	

    plot(lh3.w1.fa[[26]]$PA,tsk.w1.fa[[26]]$PA,main='win1, lh3fasta vs tskfasta')
    plot(lh3.w100.fa[[26]]$PA,tsk.w100.fa[[26]]$PA,main='win100, lh3fasta vs tskfasta')

    ngs.tsk.w100.fa <- read.psmc("test28/ngs.tsk.w100.fa.gz.stdout")


    #fa
    pdf("ngsPSMC1.pdf",width=21,height=21)
    par(mfrow=c(2,2))
    plot_ngs(read.psmc("test28/ngs.tsk.w1.fa.gz.stdout"),main="Fasta w1 for one optim")
    plot_ngs(read.psmc("test28/ngs.lh3.w100.fa.gz.stdout"),main="Fasta w100 for one optim")	
#gl
  plot_ngs(read.psmc("test28/ngs.tsk.w1.gl.gz.stdout"),main="gl w1 for one optim(based on psmcfasta)")	
  plot_ngs(read.psmc("test28/ngs.tsk.w100o.gl.gz.stdout"),main="gl w100 for one optim(based on psmcfasta)")	
  dev.off()
}



if(FALSE){
plot_chain(read.psmc("ngs.run2.stdout"))

}


if(FALSE){
    source("read.psmc.R")

 pdf('plot.pdf',width=14,height=14)
 par(mfrow=c(2,2))
 plot_chain(read.psmc("data1.res16/l0.tsk.d32.psmc.idx.stdout"),main='quad sharedstruc')
 plot_chain(read.psmc("data1.res16/l1.tsk.d32.psmc.idx.stdout"),main='lin sharedstruc')
 plot_chain(read.psmc("data1.res15/l0.tsk.d32.psmc.idx.stdout"),main='quad nosharedstruc')
 plot_chain(read.psmc("data1.res15/l1.tsk.d32.psmc.idx.stdout"),main='quad nosharedstruc')
 dev.off()

}
