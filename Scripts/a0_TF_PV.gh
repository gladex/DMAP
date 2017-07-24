
rm(list=ls())

CLine="K562"
Folder=paste("./",CLine,"/",sep="")
File<-list.files(Folder,pattern=".rda",full.names=FALSE)

load("All_BKG.rda")
s0=s; rm(s);
t0=hist(s0,prob=TRUE,breaks=40);
LIST=NULL
for (i in 1:length(File)) {
  # load("./A549/a1_ATF3_100.rda")
  load(paste(Folder,File[i],sep=""))
  TF=strsplit(File[i],"_")[[1]][2];
  
  s=x[x[,4]>0,5]; rm(S,x);
  t=hist(s,prob=TRUE,breaks=40);
  
  kstden=ks.test(t$density,t0$density);
  pv=signif(kstden$p.value,5);
  
  LIST=rbind(LIST,c(TF,pv)); # rbind in a Column
    # append(LIST,c(TF,pv))::in a Line
}
colnames(LIST)=c("TF","p-value")
write.csv(LIST,file=paste(CLine,"_PV.csv",sep=""),
          row.names=FALSE)
save(LIST,file=paste(CLine,"_PV.rda",sep=""))

#ref: write.csv(Data,file="Data.csv",row.names=FALSE)



rm(list=ls())

#par(mfrow=c(1,2))
t=hist(s,prob=TRUE,col="light blue",breaks=40,xlim=c(0,1),ylim=c(0,13),
		   xlab="Normalized methylation level",ylab="Density",main="",
			 cex.lab=1.3, cex.axis=1.3) # ,cex.main=1.5,main="ATF2 (H1-hESC cell line)"
kstden=ks.test(t$density,t0$density)
#NA: ks.test(t$counts,t0$counts)
#lines(t$breaks,append(t$density,0),lwd=3,col="red")
# col="red",density=1,breaks=10,
# breaks="Sturges",xlim = range(breaks),include.lowest=TRUE,right=TRUE,
# plot(t$breaks[-1],t$density,type="b",pch=19,lwd=2,col="red")
# lines(density(s,bw=0.05),lwd=3,col="red")
lines(t0$breaks,append(t0$density,0),lwd=2,lty=1,col="blue")
#lines(density(s0, bw=0.05),lwd=3,lty=2,col="blue")
#rug(s,ticksize=0.02,col="light blue")  
legend("topright",c("Methylation level","Background"),
			 lwd=2,lty=1,col=c('light blue','blue'),bty="n",cex=1.3)
legend("right",paste("K-S test p-vaule: ",signif(kstden$p.value,4),sep=""),bty="n",cex=1.1)


rm(list=ls())
s0=s; rm(s)
t0=hist(s0,prob=TRUE,col=NULL, border=NULL,breaks=40,xlim=c(0,1),ylim=c(0,12))
load("F:/Ubuntu.Desktop/DMR/H1.hESC/a1_NRSF_100.rda")
s=x[x[,4]>0,5]; rm(S,x)
pdf("130807_NRSF.pdf", width=5, height=5)
t=hist(s,prob=TRUE,col="light blue",breaks=40,xlim=c(0,1),ylim=c(0,12),
			 xlab="Normalized methylation level",ylab="",main="")#,ylab="Density",main="NRSF (H1-hESC cell line)"
kstden=ks.test(t$density,t0$density)
#ks.test(t$counts,t0$counts)
#lines(t$breaks,append(t$density,0),lwd=3,col="red")
# col="red",density=1,breaks=10,
# breaks="Sturges",xlim = range(breaks),include.lowest=TRUE,right=TRUE,
# plot(t$breaks[-1],t$density,type="b",pch=19,lwd=2,col="red")

#lines(density(s,bw=0.05),lwd=3,col="red")
lines(t0$breaks,append(t0$density,0),lwd=2,lty=1,col="blue")
#lines(density(s0, bw=0.05),lwd=3,lty=2,col="blue")
#rug(s,ticksize=0.02,col="light blue")  
legend("topright",c("Methylation level","Background"),
			 lwd=2,lty=1,col=c('light blue','blue'),bty="n")
legend("right",paste("K-S test p-vaule: ",signif(kstden$p.value,4),sep=""),bty="n")
