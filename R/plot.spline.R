fun<-function(x){
d<-readLines(x)
sps<-	grep("sp",d,val=T)
splinevals<-matrix(as.numeric(matrix(unlist(strsplit(grep("sp",d,val=T),"\t")),nrow=length(sps),byrow=T)[,-1]),nrow=2,byrow=F)
	d2<-strsplit(grep("fd",d,val=T),"\t")
fd<-matrix(as.numeric(sapply(d2,function(x) x[-1])),byrow=F,ncol=length(d2))
		return(list(sp=splinevals,fd=t(fd)))
		
}

spline<-function(x, spl){
   sp<-spl[x>spl[,1]&x<=spl[,2],-c(1:2)]
   return(   sum(x^(3:0)*sp))
}

plot.spline<-function(x,step=0.1){
	d<-fun(x)
	xs<-seq(0.00002,max(d$sp[,2]),step)
	ys<-sapply(xs,function(x) spline(x,splinevals))
	plot(xs,ys,type='l',lwd=3,col=3)
	points(d$fd[,1],d$fd[,2])
	return(list(xs=xs,ys=ys))
}

plot.spline(x)




