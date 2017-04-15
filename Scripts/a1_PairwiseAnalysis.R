
## OK: load HAIB/RRBS/hg19chr
rm(list=ls())

load("Haib.T47D.Rep1.rda")
x=a.sort[,1:5]; rm(a.sort) # each Methy probe=50 bp, deleate to save memory
x[,4]=0;    

load("rrbs.T47D.r1.rda")
a.sort[,4]=a.sort[,6]/100  # convert DNAm Percentage
y=a.sort[,1:5]; rm(a.sort) # each Methy probe=2 bp,
y[,5]=0; 

load("hg19chr.rda")
itv=2000
Sp=array(0,dim=c(2,itv*22))
Sr=array(0,dim=c(2,itv*22))

for(i in 1:nrow(x)) {  # Probes: [1,482421]
  idx=match(x[i,1],hg19chr[1:22,1])  # probe chrN matching peak's
  if(!is.na(idx)) { # matching case
  	probe=x[i,2] # end point [C+G]
  	j=itv*(idx-1)+ceiling(probe/(hg19chr[idx,2]/itv)) # j=index
  	Sp[1,j]=Sp[1,j]+1      # count
  	Sp[2,j]=Sp[2,j]+x[i,5] # score / x[i,5]
  	x[i,4]=j
  }
}

for (i in 1:nrow(y)) {  # Probes: [1,482421]
  idx=match(y[i,1],hg19chr[1:22,1])  # probe chrN matching peak's
	if(!is.na(idx)) { # matching case
		probe=y[i,2] # end point [C+G]
		j=itv*(idx-1)+ceiling(probe/(hg19chr[idx,2]/itv)) # j=index
		Sr[1,j]=Sr[1,j]+1      # count
		Sr[2,j]=Sr[2,j]+y[i,4] # score / x[i,5]
		y[i,5]=j
	}
}

save(Sr,Sp,file="a4_Pairwise_T47D.rda")

NZ=array(0,dim=c(2,ncol(Sp)))
for(i in 1:ncol(Sp)){
  if(Sp[1,i]>0){
    NZ[1,i]=Sp[2,i]/Sp[1,i]
  }else {
    NZ[1,i]=0 
  }
  if(Sr[1,i]>0){
    NZ[2,i]=Sr[2,i]/Sr[1,i]
  }else {
    NZ[2,i]=0
  }
}

## None-Zero Entry
p=Sp[1,]>0; r=Sr[1,]>0; int=(p & r)
#cor(NZ[1,int],NZ[2,int])

X1=NZ[1,int]; X2=NZ[2,int];dat=data.frame(X1,X2)
library(ggplot2)
tiff(file="a4_Pairwise_T47D_nz.tiff",width=4,height=4,units='in',res=300)
ggplot(dat,aes(x=X1,y=X2))+
  theme_bw() +
  xlab("450K methylation (score)")+ylab("RRBS methylation (100%)")+
  geom_point(shape=20,alpha=0.05,colour="navyblue")+ ##1177FF")+ #,color="red")+
  scale_colour_brewer()+#palette="Set2")+
  geom_density2d(color="red",width=3)+
  geom_text(color="black",size=4,  #hjust=1.0, 
            aes(x=0.1,y=0.9,label=paste('R =',signif(cor(X1,X2),5))))+
  theme(axis.title=element_text(face="bold",colour="black",size=12))+
  geom_abline(intercept=0,slope=1,colour="dodgerblue",size=1)
dev.off()


plot(NZ[1,int],NZ[2,int],pch='.',col='lightblue',
     xlab="450K methylation (score)",
     ylab="RRBS methylation (100%)",
     main="Genome-wide methylation pairwise correlation",
     xlim=c(0,1),ylim=c(0,1))
legend("topleft",paste('R =',signif(cor(NZ[1,int],NZ[2,int]),5)),bty="n")
abline(a=0,b=signif(cor(NZ[1,int],NZ[2,int]),5),col=4)
grid()



plot(1:44000,NZ[1,],pch='.',col='violet',axes=FALSE,
     xlim=c(1,44000),ylim=c(0,1),
     xlab="Chromosome (HG19: chr1~chr22)",ylab="540K methylation level",
     main="T47D (540K)")
lines(seq(1,44000,K)+(K/2),A[1,],col='black',lwd=2)
#axis(1,at=c(0,25,50),labels=c('','TSS','5 kbp'))  # col.axis="blue"
axis(1,at=seq(1,44000,2000)+1000,labels=seq(1,22,1))
axis(2);box();grid()
plot(1:44000,NZ[2,],pch='.',col='violet',axes=FALSE,
     xlim=c(1,44000),ylim=c(0,1),
     xlab="Chromosome (HG19: chr1~chr22)",ylab="RRBS methylation level",
     main="T47D (RRBS)")
lines(seq(1,44000,K)+(K/2),A[2,],col='black',lwd=2)
axis(1,at=seq(1,44000,2000)+1000,labels=seq(1,22,1))
axis(2);box();grid()

q(save="no")
