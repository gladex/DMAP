
fun_getGenesFC <- function(DMG,LOCI,HPRO) {
  load("08_19CL_RNAseq_LogFC.rda")
  #head(RNA,3)
  #  ENTREZID SYMBOL    logFC
  #1     6387 CXCL12 3.142078
  #2    79570 NKAIN1 2.646481
  
  #colnames(DMG)=c("SYMBOL","Methy")
  #head(DMG)
  #   SYMBOL    Methy
  #1 MIR200A 67.10608
  #2 MIR200B 67.10608
  
  load("Human_TSGs_1207.rda")
  #head(TSG[,1:2],3)
  # GeneID GeneSymbol
  #1    43       ACHE
  #2    95       ACY1
  
  # INEFF: rnag=RNA[RNA[,2] %in% DMG[,1], ]
  rnag=merge(RNA,DMG,by=c("SYMBOL"))
  #head(rnag) # 279x4
  #  SYMBOL ENTREZID       logFC    Methy
  #1   ABL1       25 -0.03474763 27.13594
  #2  ACAP3   116983 -0.25227989 22.86098
  #3    ACD    65057  0.18210375 25.15140
  #4  ACOT4   122970  0.57052322 20.18026
  #5   ACP1       52 -0.04805761 42.59468
  #6  ACPL2    92370 -0.09587496 25.15694

  ## LOCI
  # Cds=1;  Gene=2; Intron=3; Prom=4; Utr3=5; Utr5=6
  ## HPRO (Hyper-/Hypo-Methylation)
  # hyper=1; hypo=0
  ## TSG
  # Tsgi=1; Tsgi=0
  if(LOCI=='Cds'){
    idx=1
  } else if (LOCI=='Gene') {
    idx=2
  } else if (LOCI=='Intron'){
    idx=3
  } else if (LOCI=='Prom'){
    idx=4
  } else if (LOCI=='Utr3'){
    idx=5
  } else if (LOCI=='Utr5'){
    idx=6 
  }
  Loci=c(rep(idx,nrow(rnag)))
  
  if(HPRO=='hyper'){
    idh=1
  } else if (HPRO=='hypo') {
    idh=0
  }
  HPRO=c(rep(idh,nrow(rnag)))
  
  TSGi=as.factor(rnag[,2] %in% TSG[,1])  # ID
  ts=cbind(rnag,TSGi,Loci)
  
  return(ts)
}