library(edgeR)
library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)
library(grid)
library(data.table)
library(limma)
library(locfdr)

options(stringsAsFactors=F)
setwd("~/projects/022017_deniz/")

dat=read.csv("./051017_stanford.tsv",sep="\t",header=T)
dat=dat[dat[,3]>=2,]
dat.mat=as.matrix(dat[,c(-1,-2,-3)])
rownames(dat.mat)=unlist(strsplit(dat[,1],"\\|"))[seq(2,nrow(dat)*3,3)]
symbols=dat[,1:3]
rownames(symbols)=rownames(dat.mat)

#Setup design
samples=rep(c("Control","Puro","RNase"),each=3)
reps=rep(c(1,2,3),3)
pairs=c(1,2,3,4,4,4,1,2,3) #pairing info
design=model.matrix(~0+pairs+samples)
colnames(design)=c("block","Ctrl","Puro","RNase")

#Normalization
dat.voom=voom(dat.mat,design,normalize.method="quantile")
plotMDS(dat.voom) #Using quantile normalization based on MDS plot clustering of replicates
dat.voom.corMat=cor(dat.voom$E)
colnames(dat.voom.corMat)=paste(samples,reps,sep="_")
rownames(dat.voom.corMat)=paste(samples,reps,sep="_")
write.table(dat.voom.corMat,file="correlation_matrix.csv",sep="\t",row.names=T,col.names=T,quote=F)
  
makeCont=function(contStrings, contColumnNames) {
  cmdstr=paste("makeContrasts(",paste(contStrings,collapse=","),",levels=design)",sep="")
  cmdcont=eval(parse(text=cmdstr))
  colnames(cmdcont)=contColumnNames
  return(cmdcont)
}

#Pairwise contrast
cont=apply(t(combn(colnames(design)[-1],2)),1,paste,collapse="-")
contNames=apply(t(combn(colnames(design)[-1],2)),1,paste,collapse="")
cont=paste("(",cont,")",sep="")
cont=makeCont(cont,contNames)

#LM fit and eBayes
fit=lmFit(dat.voom,design)
cont=makeContrasts(RNase-Ctrl,Puro-Ctrl,levels=design)
colnames(cont)=c("RNaseCtrl","PuroCtrl")
fit2=eBayes(contrasts.fit(fit,cont),proportion=0.15)

#Unmod T-stat/pval
fit2$unmodT=fit2$coefficients/fit2$stdev.unscaled/fit2$sigma
fit2$unmodPval=2*pt(-abs(fit2$unmodT), fit2$df.residual)

#Calculate 95% confidence intervals for all coefs
fit2$margin.error=sqrt(fit2$s2.post)*fit2$stdev.unscaled*qt(0.975,df=fit2$df.total) #Probably would want to do this again based on eNull
fit2$CI.L=fit2$coef-fit2$margin.error
fit2$CI.R=fit2$coef+fit2$margin.error

fit2.CI.L=fit2$CI.L
fit2.CI.R=fit2$CI.R
colnames(fit2.CI.L)=paste("CI.L.",colnames(fit2.CI.L),sep="")
colnames(fit2.CI.R)=paste("CI.R.",colnames(fit2.CI.R),sep="")

#Local FDR / empirical null modeling
locfdr.RNase=locfdr(fit2$t[,1],type=1,df=9,nulltype=1,bre=75,mlests=c(1,7)) #Using MLE, based off of starting points found in CME
locfdr.Puro=locfdr(fit2$t[,2],type=1,df=9,nulltype=1,bre=75)
fit2$locfdr=fit2$p.value
fit2$locfdr[,1]=locfdr.RNase$fdr
fit2$locfdr[,2]=locfdr.Puro$fdr

#Get p-value, two sided (go same distance both direction from null t-val center)
getLeftPRev=function(locfdrObj,tval) {
  center=locfdrObj$fp0["mlest","delta"]
  tval.abs=abs(tval-center)
  x.1=approx(x=locfdrObj$mat[,"x"],y=cumsum(locfdrObj$mat[,"f0"]/sum(locfdrObj$mat[,"f0"])),center-tval.abs,rule=2)$y
  x.2=approx(x=locfdrObj$mat[,"x"],y=1-cumsum(locfdrObj$mat[,"f0"]/sum(locfdrObj$mat[,"f0"])),center+tval.abs,rule=2)$y
  return(x.1+x.2)
}

#FDR and NPV calculations only counting for lost interactions as "positive hits"

getLeftFdr=function(locfdrObj,cutoff) {
  y=approx(y=locfdrObj$mat[,"x"],x=cumsum((locfdrObj$mat[,"f"]-locfdrObj$mat[,"p1f1"])/cumsum(locfdrObj$mat[,"f"])),cutoff,rule=2)$y
  return(y)
}

getLeftFdrRev=function(locfdrObj,tval) {
  x=approx(x=locfdrObj$mat[,"x"],y=cumsum((locfdrObj$mat[,"f"]-locfdrObj$mat[,"p1f1"])/cumsum(locfdrObj$mat[,"f"])),tval,rule=2)$y
  if(x>=1) {
    x=1
  }
  return(x)
}

getRightNPV=function(locfdrObj,cutoff) {
  ycenter=order(locfdrObj$mat[,"f0"],decreasing=T)[1]
  yp1f1=c((locfdrObj$mat[,"p1f1"])[1:ycenter],rep(0,nrow(locfdrObj$mat)-ycenter))
  y=approx(rev(locfdrObj$mat[,"x"]),x=cumsum(rev((locfdrObj$mat[,"f"]-yp1f1)))/cumsum(rev(locfdrObj$mat[,"f"])),cutoff,rule=2)$y
  
  return(y)
}

getRightNPVRev=function(locfdrObj,tval) {
  ycenter=order(locfdrObj$mat[,"f0"],decreasing=T)[1]
  yp1f1=c((locfdrObj$mat[,"p1f1"])[1:ycenter],rep(0,nrow(locfdrObj$mat)-ycenter))
  x=approx(x=rev(locfdrObj$mat[,"x"]),y=cumsum(rev((locfdrObj$mat[,"f"]-yp1f1)))/cumsum(rev(locfdrObj$mat[,"f"])),tval,rule=2)$y
  
  return(x)
}

getP1f1=function(locfdrObj, centerProp) {
  locfdrObj.f0choose.center=order(locfdrObj$mat[,"f0"],decreasing=T)[1]
  locfdrObj.f0choose.center.props=data.frame()
  for (x in 1:min(nrow(locfdrObj$mat)-locfdrObj.f0choose.center,locfdrObj.f0choose.center)) {
    y=data.frame(x,sum(locfdrObj$mat[,"p1f1"][(locfdrObj.f0choose.center-x):(locfdrObj.f0choose.center+x)]/sum(locfdrObj$mat[,"f"])))
    colnames(y)=c("x","p")
    locfdrObj.f0choose.center.props=rbind(locfdrObj.f0choose.center.props,y)
  }
  locfdrObj.f0choose.center.props=locfdrObj.f0choose.center.props[order(locfdrObj.f0choose.center.props$p,decreasing=T),]
  locfdrObj.f0choose=locfdrObj.f0choose.center.props[locfdrObj.f0choose.center.props$p<=centerProp,1][1]
  locfdrObj.f0choose.tmin=min(locfdrObj$mat[(locfdrObj.f0choose.center-locfdrObj.f0choose):(locfdrObj.f0choose.center+locfdrObj.f0choose),"x"])
  locfdrObj.f0choose.tmax=max(locfdrObj$mat[(locfdrObj.f0choose.center-locfdrObj.f0choose):(locfdrObj.f0choose.center+locfdrObj.f0choose),"x"])
  return(c(locfdrObj.f0choose.tmin,locfdrObj.f0choose.tmax))
}

#Benjamini-Hochberg FDR estimates for theoretical; will not use
fit2$p.value.adj=apply(fit2$p.value,2,p.adjust,method="BH")

#Empirical null p-value
fit2.pval.enull=data.frame(sapply(fit2$t[,1],getLeftPRev,locfdrObj=locfdr.RNase),sapply(fit2$t[,2],getLeftPRev,locfdrObj=locfdr.Puro))
fit2.leftFdr.enull=data.frame(sapply(fit2$t[,1],getLeftFdrRev,locfdrObj=locfdr.RNase),sapply(fit2$t[,2],getLeftFdrRev,locfdrObj=locfdr.Puro))
fit2.rightNpv.enull=data.frame(sapply(fit2$t[,1],getRightNPVRev,locfdrObj=locfdr.RNase),sapply(fit2$t[,2],getRightNPVRev,locfdrObj=locfdr.Puro))
colnames(fit2.pval.enull)=paste("loc.leftPval.eNull.",colnames(fit2$p.value),sep="")
colnames(fit2.leftFdr.enull)=paste("loc.leftFdr.eNull.",colnames(fit2$p.value),sep="")
colnames(fit2.rightNpv.enull)=paste("loc.RightNpv.eNull.",colnames(fit2$p.value),sep="")

#Some formatting and making a "master table"
fit2.pval=fit2$p.value
fit2.pval.adj=fit2$p.value.adj
fit2.pval.locfdr=fit2$locfdr
fit2.coef=fit2$coef
fit2.t=fit2$t
colnames(fit2.pval)=paste("P.Value.",colnames(fit2.pval),sep="")
colnames(fit2.pval.adj)=paste("adj.P.Val.",colnames(fit2.pval.adj),sep="")
colnames(fit2.pval.locfdr)=paste("loc.fdr.",colnames(fit2.pval.locfdr),sep="")
colnames(fit2.coef)=paste("beta.",colnames(fit2.coef),sep="")
colnames(fit2.t)=paste("t.",colnames(fit2.t),sep="")

sigTable=data.frame(cbind(rownames(fit2),fit2.coef,fit2.t,fit2.pval.enull,fit2.leftFdr.enull,fit2.rightNpv.enull)) #The master table
sigTable[,-1]=sapply(sigTable[,-1],as.numeric)
rownames(sigTable)=rownames(fit2$coef)
colnames(sigTable)[1]="UniProt"

write.table(sigTable,file="masterTable.tsv",row.names=F,col.names=T,quote=F,sep="\t")
write.table(cbind(dat.mat[sigTable[,1],],symbols[sigTable[,1],1:3],sigTable),file="masterTable2.tsv",row.names=F,col.names=T,quote=F,sep="\t")

#F0 mean and SD
rnase.null2sd=sigTable$beta.RNaseCtrl[(sigTable$t.RNaseCtrl<=getP1f1(locfdr.RNase,0.01)[2])&(sigTable$t.RNaseCtrl>=getP1f1(locfdr.RNase,0.01)[1])]
puro.null2sd=sigTable$beta.PuroCtrl[(sigTable$t.PuroCtrl<=getP1f1(locfdr.Puro,0.01)[2])&(sigTable$t.PuroCtrl>=getP1f1(locfdr.Puro,0.01)[1])]


nullMean.df=data.frame(c("RNase.f0.mean","RNase.f0.sd","Puro.f0.mean","Puro.f0.sd"),c(mean(rnase.null2sd),sd(rnase.null2sd),mean(puro.null2sd),sd(puro.null2sd)))
colnames(nullMean.df)=c("name","value")
write.table(nullMean.df,file="f0_mean_sd.tsv",sep="\t",quote=F,row.names=F,col.names=F)

#Make locfdr plots
locfdr.RNase.mat=locfdr.RNase$mat
locfdr.RNase.mat[,"f0"]=locfdr.RNase.mat[,"f"]-locfdr.RNase.mat[,"p1f1"]
locfdr.Puro.mat=locfdr.Puro$mat
locfdr.Puro.mat[,"f0"]=locfdr.Puro.mat[,"f"]-locfdr.Puro.mat[,"p1f1"]

locfdr.RNase.melt=melt(data.frame(locfdr.RNase.mat),id.vars="x")
locfdr.Puro.melt=melt(data.frame(locfdr.Puro.mat),id.vars="x")

locfdr.RNase.NPV.cutoff=getRightNPV(locfdr.RNase,0.99)
locfdr.Puro.NPV.cutoff=getRightNPV(locfdr.Puro,0.99)

locfdr.RNase.FDR.cutoff=getLeftFdr(locfdr.RNase,0.15)
locfdr.Puro.FDR.cutoff=getLeftFdr(locfdr.Puro,0.15)

sigTable$NPV_cutoff.RNase=(sigTable$t.RNaseCtrl>=locfdr.RNase.NPV.cutoff)
sigTable$NPV_cutoff.Puro=(sigTable$t.PuroCtrl>=locfdr.Puro.NPV.cutoff)

sigTable$FDR_cutoff.RNase=(sigTable$t.RNaseCtrl<=locfdr.RNase.FDR.cutoff)
sigTable$FDR_cutoff.Puro=(sigTable$t.PuroCtrl<=locfdr.Puro.FDR.cutoff)

ggplot(sigTable,aes(x=-beta.RNaseCtrl,y=-log(loc.leftPval.eNull.RNaseCtrl,10),colour=paste(NPV_cutoff.RNase,FDR_cutoff.RNase)))+theme_classic()+geom_point()

ggplot(sigTable,aes(x=-beta.PuroCtrl,y=-log(loc.leftPval.eNull.PuroCtrl,10),colour=paste(NPV_cutoff.Puro,FDR_cutoff.Puro)))+theme_classic()+geom_point()

ggplot(sigTable,aes(x=-beta.RNaseCtrl,y=abs(t.RNaseCtrl),colour=paste(NPV_cutoff.RNase,FDR_cutoff.RNase)))+theme_classic()+geom_point()

pdf("rnase.pdf",width=6,height=4)
ggplot()+theme_classic()+
  geom_histogram(data=sigTable,aes(x=t.RNaseCtrl),binwidth=mean(diff(locfdr.RNase$mat[,1])),colour="black",fill="grey",alpha=0.5,size=0.3)+
  geom_density(data=subset(locfdr.RNase.melt,variable%in%c("p1f1")),stat="identity",aes(x=x,y=value,fill=variable),size=0.5,alpha=0.5)+
  geom_line(data=subset(locfdr.RNase.melt,variable%in%c("f","f0")),aes(x=x,y=value,group=variable,colour=variable,linetype=variable),size=1,alpha=0.75)+
  geom_vline(xintercept=locfdr.RNase.FDR.cutoff,colour="red")+
  geom_vline(xintercept=locfdr.RNase.NPV.cutoff,colour="blue")+
  xlab("T-statistic")+ylab("Count")+
  scale_y_continuous(expand=c(0,0))+
  scale_colour_manual(values=c("darkgreen","blue"),name=NULL,labels=c("Mixture","Null"))+
  scale_linetype_discrete(labels=c("Mixture","Null"),name=NULL)+
  scale_fill_manual(values=c("red"),name=NULL,labels=c("Non-null"))+
  annotate("text",x=diff(range(sigTable$t.RNaseCtrl))*0.1+min(sigTable$t.RNaseCtrl),y=max(locfdr.RNase$mat[,"counts"])*0.9,label=paste("1-p0=",round(1-locfdr.RNase$fp0["mlest","p0"],digits=2),sep=""))+
  theme(legend.position=c(0.9,0.8))
dev.off()

pdf("puro.pdf",width=6,height=4)
ggplot()+theme_classic()+
  geom_histogram(data=sigTable,aes(x=t.PuroCtrl),binwidth=mean(diff(locfdr.Puro$mat[,1])),colour="black",fill="grey",alpha=0.5,size=0.3)+
  geom_density(data=subset(locfdr.Puro.melt,variable%in%c("p1f1")),stat="identity",aes(x=x,y=value,fill=variable),size=0.5,alpha=0.5)+
  geom_line(data=subset(locfdr.Puro.melt,variable%in%c("f","f0")),aes(x=x,y=value,group=variable,colour=variable,linetype=variable),size=1,alpha=0.75)+
  geom_vline(xintercept=locfdr.Puro.FDR.cutoff,colour="red")+
  geom_vline(xintercept=locfdr.Puro.NPV.cutoff,colour="blue")+
  xlab("T-statistic")+ylab("Count")+
  scale_y_continuous(expand=c(0,0))+
  scale_colour_manual(values=c("darkgreen","blue"),name=NULL,labels=c("Mixture","Null"))+
  scale_linetype_discrete(labels=c("Mixture","Null"),name=NULL)+
  scale_fill_manual(values=c("red"),name=NULL,labels=c("Non-null"))+
  annotate("text",x=diff(range(sigTable$t.PuroCtrl))*0.1+min(sigTable$t.PuroCtrl),y=max(locfdr.Puro$mat[,"counts"])*0.9,label=paste("1-p0=",round(1-locfdr.Puro$fp0["mlest","p0"],digits=2),sep=""))+
  theme(legend.position=c(0.9,0.8))
dev.off()


ggplot(sigTable,aes(x=-beta.Puggplot(sigTable,aes(x=-beta.RNaseCtrl,y=-log(loc.leftPval.eNull.RNaseCtrl,10),colour=paste(NPV_cutoff.RNase,FDR_cutoff.RNase)))+theme_classic()+geom_point()
roCtrl,y=-log10(P.Value.PuroCtrl),colour=paste(NPV_cutoff.Puro,FDR_cutoff.Puro)))+theme_classic()+geom_point()


#### Extra stuff ####

#Run revious experiment
silac=read.csv("silac_rnase.csv",sep="\t",header=F)
rownames(silac)=silac[,1]
samples=rep(c("Control","RNase"),each=2)
reps=rep(c(1,2),2)
design=model.matrix(~0+samples)
colnames(design)=c("Ctrl","RNase")

#Voom+quantile normalize
dat.voom=voom(silac[,-1],design,normalize.method="quantile")

#LM fit and eBayes
fit=lmFit(dat.voom,design)
cont=makeContrasts(RNase-Ctrl,levels=design)
colnames(cont)=c("RNaseCtrlSILAC")
fit2=eBayes(contrasts.fit(fit,cont))

#Calculate 95% confidence intervals for all coefs
fit2$margin.error=sqrt(fit2$s2.post)*fit2$stdev.unscaled*qt(0.975,df=fit2$df.total)
fit2$CI.L=fit2$coef-fit2$margin.error
fit2$CI.R=fit2$coef+fit2$margin.error
fit2.CI.L=fit2$CI.L
fit2.CI.R=fit2$CI.R
colnames(fit2.CI.L)=paste("CI.L.",colnames(fit2.CI.L),sep="")
colnames(fit2.CI.R)=paste("CI.R.",colnames(fit2.CI.R),sep="")

#Local FDR
locfdr.RNase=locfdr(fit2$t[,1],type=1,df=10,nulltype=3)
fit2$locfdr=fit2$p.value
fit2$locfdr[,1]=locfdr.RNase$fdr

#Add Benjamini-Hochberg FDR estimates
fit2$p.value.adj=apply(fit2$p.value,2,p.adjust,method="BH")

#Some formatting and making a "master table"
fit2.pval=fit2$p.value
fit2.pval.adj=fit2$p.value.adj
fit2.pval.locfdr=fit2$locfdr
fit2.coef=fit2$coef
fit2.t=fit2$t
colnames(fit2.pval)=paste("P.Value.",colnames(fit2.pval),sep="")
colnames(fit2.pval.adj)=paste("adj.P.Val.",colnames(fit2.pval.adj),sep="")
colnames(fit2.pval.locfdr)=paste("loc.fdr.",colnames(fit2.pval.locfdr),sep="")
colnames(fit2.coef)=paste("beta.",colnames(fit2.coef),sep="")
colnames(fit2.t)=paste("t.",colnames(fit2.t),sep="")

sigTable2=data.frame(cbind(rownames(fit2),fit2.coef,fit2.CI.L,fit2.CI.R,fit2.t,fit2.pval,fit2.pval.adj,fit2.pval.locfdr)) #The master table
sigTable2[,-1]=sapply(sigTable2[,-1],as.numeric)
rownames(sigTable2)=rownames(fit2$coef)
colnames(sigTable2)[1]="UniProt"

sigTable.join=merge(sigTable,sigTable2,by="UniProt",all=T)
RNaseCtrlMin=min(sigTable.join$t.RNaseCtrl,na.rm=T)-5
RNaseCtrlSILACMin=min(sigTable.join$t.RNaseCtrlSILAC,na.rm=T)-3
sigTable.join$t.RNaseCtrl[is.na(sigTable.join$t.RNaseCtrl)]=RNaseCtrlMin
sigTable.join$t.RNaseCtrlSILAC[is.na(sigTable.join$t.RNaseCtrlSILAC)]=
RNaseCtrlSILACMin
t.tmt.fdr.low=max(sigTable.join$t.RNaseCtrl[(sigTable.join$loc.fdr.RNaseCtrl<=0.1)&(sigTable.join$beta.RNaseCtrl<=0)],na.rm=T)
t.tmt.fdr.high=min(sigTable.join$t.RNaseCtrl[(sigTable.join$loc.fdr.RNaseCtrl<=0.1)&(sigTable.join$beta.RNaseCtrl>=0)],na.rm=T)
t.silac.fdr.low=max(sigTable.join$t.RNaseCtrlSILAC[(sigTable.join$loc.fdr.RNaseCtrlSILAC<=0.1)&(sigTable.join$beta.RNaseCtrlSILAC<=0)],na.rm=T)
t.silac.fdr.high=min(sigTable.join$t.RNaseCtrlSILAC[(sigTable.join$loc.fdr.RNaseCtrlSILAC<=0.1)&(sigTable.join$beta.RNaseCtrlSILAC>=0)],na.rm=T)

sigTable.join$group="Neither"
sigTable.join[rownames(subset(sigTable.join,(loc.fdr.RNaseCtrl<=0.1)&(loc.fdr.RNaseCtrlSILAC<=0.1))),"group"]="Both"
sigTable.join[rownames(subset(sigTable.join,(loc.fdr.RNaseCtrl<=0.1)&(loc.fdr.RNaseCtrlSILAC>0.1))),"group"]="Only new RNase"
sigTable.join[rownames(subset(sigTable.join,(loc.fdr.RNaseCtrl>0.1)&(loc.fdr.RNaseCtrlSILAC<=0.1))),"group"]="Only prev RNase"
sigTable.join[rownames(subset(sigTable.join,t.RNaseCtrl==RNaseCtrlMin)),"group"]="Not detected in new RNase"
sigTable.join[rownames(subset(sigTable.join,t.RNaseCtrlSILAC==RNaseCtrlSILACMin)),"group"]="Not detected in prev RNase"
ggplot(sigTable.join,aes(x=t.RNaseCtrl,y=t.RNaseCtrlSILAC,fill=group))+theme_classic()+geom_point(aes(shape=group),alpha=0.7)+geom_vline(xintercept=c(t.tmt.fdr.low,t.tmt.fdr.high))+geom_hline(yintercept=c(t.silac.fdr.low,t.silac.fdr.high))+scale_fill_manual(values=c("red","orange","black","black","darkgreen","blue"))+scale_shape_manual(values=c(21,21,3,4,21,21))
ggplot(sigTable.join,aes(x=beta.RNaseCtrl,y=beta.RNaseCtrlSILAC,fill=group))+theme_classic()+geom_point(aes(shape=group),alpha=0.7)+scale_fill_manual(values=c("red","orange","black","black","darkgreen","blue"))+scale_shape_manual(values=c(21,21,3,4,21,21))


#### PEPTIDE DIST ####
orflen=data.frame(fread("./orflen.txt",sep="\t",header=F))
rownames(orflen)=orflen[,1]

orflen2=data.frame(fread("./uniprot_trembl_use.len",sep="\t",header=F))
rownames(orflen2)=orflen2[,1]

noRp=data.frame(fread("./all_peptides_norp.bed.csv",sep="\t",header=F))
allPept=data.frame(fread("./all_peptides.bed",sep="\t",header=F))
rp=data.frame(fread("./all_peptides_rp.bed.csv",sep="\t",header=F))

plotDens=function(allPept,sizeRange,binning,center,uniqueProt,geneNorm) {
  allPept=allPept[!allPept[,1]=="",]
  allPept$len=rbind(orflen,orflen2)[allPept[,1],2]
  sizeMin=sizeRange[1]
  sizeMax=sizeRange[2]
  allPept=subset(allPept,(len>=sizeMin)&(len<=sizeMax))
  allPept[,2]=allPept[,2]/allPept$len*10
  allPept[,3]=allPept[,3]/allPept$len*10
  allPept$center=(allPept[,2]+allPept[,3])/2
  
  allPept.unique=allPept
  if(uniqueProt) {
    allPept.unique=unique(allPept)
  }
  
  allPept.unique=allPept.unique[!apply(is.na(allPept.unique),1,any),]
  avgsize=mean(allPept.unique[,3]-allPept.unique[,2])
  
  if(binning) {
    binMatrix=matrix(nrow=length(unique(allPept.unique[,1])),ncol=13,0)
  } else {
    allPept.unique[,2]=allPept.unique[,2]*allPept.unique$len/10
    allPept.unique[,3]=allPept.unique[,3]*allPept.unique$len/10
    binMatrix=matrix(nrow=length(unique(allPept.unique[,1])),ncol=max(allPept.unique$len)+2,0)
  }
  
  rownames(binMatrix)=unique(allPept.unique[,1])
  if(center) {
    for(i in 1:nrow(allPept.unique)) {
      if(binning) {
        pos=.bincode(allPept.unique[i,2]:allPept.unique[i,3],0:ncol(binMatrix))
      } else {
        pos=allPept.unique[i,2]:allPept.unique[i,3]
      }
      uniprotid=allPept.unique[i,1]
      binMatrix[uniprotid,pos]=binMatrix[uniprotid,pos]+1
    }
  } else {
    for(i in 1:nrow(allPept.unique)) {
      if(binning) {
        pos=.bincode(max(c(0,(allPept.unique[i,5]-avgsize/2))):min(c((allPept.unique[i,5]+avgsize/2)),ncol(binMatrix)-2),-1:ncol(binMatrix))
      } else {
      pos=allPept.unique[i,2]:allPept.unique[i,3]
    }
      uniprotid=allPept.unique[i,1]
      binMatrix[uniprotid,pos]=binMatrix[uniprotid,pos]+1
    }
  }

  binMatrix=binMatrix[,1:(ncol(binMatrix))]
  if(geneNorm) {
    binMatrix=binMatrix/matrix(rep(rowSums(binMatrix),ncol(binMatrix)),nrow=nrow(binMatrix),ncol=ncol(binMatrix),byrow=F)
  }
  binMatrix.df=data.frame(x=1:ncol(binMatrix),colSums(binMatrix))
  colnames(binMatrix.df)=c("bin","count")
  binMatrix.df[,1]=binMatrix.df[,1]-1
  
  #fig=ggplot(binMatrix.df,aes(x=bin,y=count))+theme_classic()+geom_line()+scale_y_continuous(expand=c(0,0))+xlab("Bin")+ylab("Count")
#  return(fig)
  return(binMatrix.df)
}

pdf("density_nounique_nocenter_norownorm.pdf",width=6,height=4)
plotDens(F,F,F)
dev.off()
pdf("density_unique_nocenter_norownorm.pdf",width=6,height=4)
plotDens(T,F,F)
dev.off()
pdf("density_unique_center_norownorm.pdf",width=6,height=4)
plotDens(T,T,F)
dev.off()
pdf("density_unique_center_rownorm.pdf",width=6,height=4)
plotDens(T,T,T)
dev.off()
pdf("density_nounique_center_norownorm.pdf",width=6,height=4)
plotDens(F,T,F)
dev.off()
pdf("density_nounique_center_rownorm.pdf",width=6,height=4)
plotDens(F,T,T)
dev.off()
pdf("density_nounique_nocenter_rownorm.pdf",width=6,height=4)
plotDens(F,F,T)
dev.off()
pdf("density_unique_nocenter_rownorm.pdf",width=6,height=4)
plotDens(T,F,T)
dev.off()

plotDens(allPept,600,700,T,F,F,T)
bins.df=merge(plotDens(noRp,c(0,Inf),T,F,T,T),plotDens(rp,c(0,Inf),T,F,T,T),all=T,by="bin")
bins.df[2,2:3]=bins.df[2,2:3]+bins.df[1,2:3]
bins.df=bins.df[2:(nrow(bins.df)-2),]
colnames(bins.df)=c("bin","noRp","rp")
bins.df=melt(bins.df,id.vars="bin")
bins.df$bin=factor(bins.df$bin,levels=unique(bins.df$bin))
colnames(bins.df)=c("bin","type","weighted_count")


pdf("norp.pdf",width=6,heigh=4)
ggplot(subset(bins.df,type=="noRp"),aes(x=bin,y=weighted_count))+theme_classic()+geom_bar(stat="identity",fill="darkred")+xlab("Bin")+ylab("Weighted count")+scale_y_continuous(expand=c(0,0))
dev.off()

pdf("rp.pdf",width=6,heigh=4)
ggplot(subset(bins.df,type=="rp"),aes(x=bin,y=weighted_count))+theme_classic()+geom_bar(stat="identity",fill="royalblue")+xlab("Bin")+ylab("Weighted count")+scale_y_continuous(expand=c(0,0))
dev.off()
