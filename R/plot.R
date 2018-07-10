getVec <- function(x){
file<-x
d<-as.numeric(unlist(strsplit(system(paste0("grep PA ",file," | tail -n1|cut -f2|cut -f1-4 --complement -d\" \""),intern=T)," ")))
d
}
