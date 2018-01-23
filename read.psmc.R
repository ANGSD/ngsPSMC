

read.psmc1 <- function(x,rd=-1){
    v<-system(paste0("grep -n RD ",x,"|grep -v C|tr \":\" \"\t\"|cut -f2 --complement"),intern=T)
    rdsplit <-  matrix(as.numeric(unlist(strsplit(v,"\t"))),ncol=2,byrow=T)[,2:1]
    rdsplit <- rbind(rdsplit,c(-1,as.numeric(system(paste0("wc -l ",x,"| cut -f1 -d\" \""),intern=T))))
    if(rd==-1)
        rd <- tail(rdsplit,2)[1,1]
##    cat("\t-> Numaber of RD from file: ",x," is ",nrow(rdsplit)-2,"\n")
    wh <-which(rdsplit[,1]==rd)
    until <- rdsplit[wh+1,2]
    howMany <- until-rdsplit[wh,2]
    cmd <- paste("head -n",until,x ,"|tail -n",howMany)
    TR <-c(NA,NA);
    TR[1] <- as.numeric(system(paste0(cmd,"|grep TR|cut -f2"),intern=T))
    TR[2] <- as.numeric(system(paste0(cmd,"|grep TR|cut -f3"),intern=T))
    MT <- as.numeric(system(paste0(cmd,"|grep MT|cut -f2"),intern=T))
    RS <- t(matrix(as.numeric(unlist(strsplit(system(paste0(cmd,"|grep RS|cut -f2,3,4"),intern=T),"\t"))),3))
    PAT <- system(paste0(cmd,"|grep PA|cut -f2|cut -f1 -d\" \""),intern=T)
    PA <-  as.numeric(unlist(strsplit(system(paste0(cmd,"|grep PA|cut -f2|cut -f1 -d\" \" --complement"),intern=T)," ")))
    ret <- list(TR=TR,MT=MT,RS=RS,PAT=PAT,PA=PA)
    return(ret)
    
}

read.psmc <- function(x){
   v<-system(paste0("grep -n RD ",x,"|grep -v C|tr \":\" \"\t\"|cut -f2 --complement"),intern=T)
    rdsplit <-  matrix(as.numeric(unlist(strsplit(v,"\t"))),ncol=2,byrow=T)[,2:1]
    rdsplit <- rbind(rdsplit,c(-1,as.numeric(system(paste0("wc -l ",x,"| cut -f1 -d\" \""),intern=T))))
   cat("\t-> Number of RD from file: ",x," is ",nrow(rdsplit)-2,"\n")
   ret <- lapply(0:(nrow(rdsplit)-2),function(rd) read.psmc1(x,rd))
   return(ret)
}


if(FALSE){
    source("read.psmc.R")
    x<-"data_std_nomigr2/ms2g2.psmc"
    ret <- read.psmc(x)

}
