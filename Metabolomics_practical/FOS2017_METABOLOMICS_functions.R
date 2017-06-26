
### Defining functions:

## Defining a function for 
print.flag.report <- function(phen){
  IND <- grep("abnormal_macromolecule_a|low_glucose|low_glutamine_high_glutamate|low_protein_content|high_citrate|high_ethanol|high_lactate|high_pyruvate|serum_sample|unidentified_small_molecule_a|unidentified_small_molecule_b|unknown_acetylated_compound|isopropyl_alcohol|polysaccharides|aminocaproic_acid|fast",colnames(phen))
  Dummi <- colSums(phen[,IND])
  ANYFLAG <- length(which(rowSums(phen[,IND])>0))
  NOFLAG <- length(which(rowSums(phen[,IND])==0))
  pander:::pander(data.frame(QCflag=c(names(Dummi),"ANYFLAG","NOFLAG"),counts=c(as.vector(Dummi),ANYFLAG,NOFLAG),stringsAsFactors=FALSE))
}

### Defining a function for plotting a heatmap indcating missing & zero values:
plot.na.heatmap  <- function(dat){
  Dummi <- matrix(NA,ncol=ncol(dat),nrow=nrow(dat))
  Dummi[which(dat==0)] <- 0
  Dummi[which(dat!=0)] <- 1
  layout(mat=matrix(c(1,2,3,4),ncol=2,byrow=TRUE),widths=c(4,1),heights=c(1,4))
  par(xaxs = "i")  # THIS IS EVIL MAGIC!
  par(yaxs = "i")
  par(mar=c(0,3,2,0))
  XCOUNT <- nrow(Dummi)-colSums(Dummi,na.rm=TRUE)
  XPERC  <- XCOUNT/nrow(Dummi)*100
  YCOUNT <- ncol(Dummi)-rowSums(Dummi,na.rm=TRUE)
  YPERC  <- YCOUNT/ncol(Dummi)*100
  if(max(XPERC)<10){
    ylim <- c(0,10)
  } else {
    ylim <- c(0,max(pretty(XPERC)))
  }
  par(mar=c(0,3,4,0))
  barplot(XPERC,axes=TRUE,main="Missingness & Zeros Plot",col="lightblue",border="lightblue",las=2,cex.main=4,ylim=ylim,cex.axis=1.5,font.axis=2)
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1))
  legend("topright",fill=c("grey90","grey30","red"),legend=c("missing","value","zero"))
  par(mar=c(3,3,0,0))
  image(t(Dummi[nrow(Dummi):1,]),axes=FALSE,col=c(rgb(1,0,0,1),"grey30"),xlim=c(0,1),ylim=c(0,1))
  mtext(text="measurements",side=2,font=2,cex=1.5)
  mtext(text="samples",side=1,font=2,cex=1.5)
  par(mar=c(3,0,0,0))
    if(max(YPERC)<10){
    xlim <- c(0,10)
  } else {
    xlim <- c(0,max(pretty(YPERC)))
  }
  barplot(rev(YPERC),horiz=TRUE,axes=TRUE,col="lightblue",border="lightblue",cex.axis=1.5,font.axis=2)
}

## Defining a function for counting missing values:
count.miss <- function(dat,on_sample=TRUE){
  #IND <- grep("abnormal_macromolecule_a|low_glucose|low_glutamine_high_glutamate|low_protein_content|high_citrate|high_ethanol|high_lactate|high_pyruvate|serum_sample|unidentified_small_molecule_a|unidentified_small_molecule_b|unknown_acetylated_compound|isopropyl_alcohol|polysaccharides|aminocaproic_acid|fast",rownames(dat))
  #dat <- dat[-IND,]
  Dummi <- matrix(NA,ncol=ncol(dat),nrow=nrow(dat))
  Dummi[which(dat==0)] <- 0
  Dummi[which(dat!=0)] <- 1
  if(on_sample){
    COUNTS <- list(miss=colSums(is.na(Dummi)),zero=colSums(Dummi==0,na.rm=TRUE),
                   miss_or_zero=nrow(Dummi)-colSums(Dummi==1,na.rm=TRUE),Nmeas=nrow(Dummi),Nsamp=ncol(Dummi))
  } else {
    COUNTS <- list(miss=rowSums(is.na(Dummi)),zero=rowSums(Dummi==0,na.rm=TRUE),
                   miss_or_zero=ncol(Dummi)-rowSums(Dummi==1,na.rm=TRUE),Nmeas=nrow(Dummi),Nsamp=ncol(Dummi))
  }
  return(COUNTS)
}

## Defining afunction for
print.miss.report <- function(dat,on_sample=TRUE,type="missingOrZero"){
  ## Defining a function to report on missingness given counts:
  rep.miss <- function(counts,Nmeas,Nsamp,on_sample){
    A <- table(counts)
    B <- data.frame(as.numeric(names(A)),
        as.numeric(names(A))/ifelse(on_sample,Nmeas,Nsamp)*100,
        as.vector(A),
        as.vector(A)/ifelse(on_sample,Nsamp,Nmeas)*100,stringsAsFactors=FALSE)
    if(on_sample){
      colnames(B) <- c("# missMeas","[%]","# samp","[%]")
    } else {
      colnames(B) <- c("# missSamp","[%]","# meas","[%]")
    }
    return(B)
  }
  COUNTS <- count.miss(dat,on_sample=on_sample)
  if(type=="missingOrZero"){
    res <- rep.miss(COUNTS[["miss_or_zero"]],Nmeas=COUNTS$Nmeas,Nsamp=COUNTS$Nsamp,on_sample=on_sample)
  }
  if(type=="missing"){
    res <- rep.miss(COUNTS[["miss"]],Nmeas=COUNTS$Nmeas,Nsamp=COUNTS$Nsamp,on_sample=on_sample)
  }
  if(type=="zero"){
    res <- rep.miss(COUNTS[["zero"]],Nmeas=COUNTS$Nmeas,Nsamp=COUNTS$Nsamp,on_sample=on_sample)
  }
  pander:::pander(res,justify = c('left',rep('right',3)),caption=paste0(type," in ",ifelse(on_sample,"samples","measurements")),digits=3)
}

## Defining a function for Rank Inverse Normal transformation:
RIN <- function(x){
  ## Defining a function for Rank Inverse Normal transformation of a single measurement:
  single.RIN <- function(x){
    x <- rank(x, "keep")
    x  <- (x - 0.5) / length(x)
    return(qnorm(x))
  }
  if(is.matrix(x)){
    x <- t(apply(x,1,single.RIN))
    return(x)
  } else {
    if(is.vector(x)){
      return(single.RIN(x))
    }
  }
}

## Defining a function for performing a metabo-was:
do.metabowas <- function(phen,dat,test_variable="age",covariates=c("sex")){
  vars <- phen[, c(test_variable, covariates)]
  vars <- na.omit(vars)
  dat <- dat[, match(rownames(vars), colnames(dat))]
  design <- model.matrix(~ ., vars)
  fit <- lmFit(dat, design)
  fit <- eBayes(fit)
  result <- limma::topTable(fit, coef=2, number=Inf)
  colnames(result) <- c("estimate","AveEpr","tstat","pval","pval.adj")
  result <- result[,c("estimate","tstat","pval","pval.adj")]
  return(result)
}

## Defining a function for plotting beta vs beta's:
plot.beta_beta <- function(res1,res2,pcut=0.05,ptype="bonf",xlim=NULL,ylim=xlim,main=NULL,xlab=NULL,ylab=NULL,
      col=c("darkorchid","brown3","cornflowerblue","black"),pch=c(22,21,21,20),cex=c(1.5,1,1,1)){
  # Merge results:
  res1 <- res1[,c("estimate","pval")]
  colnames(res1) <- c("estimate1","pval1")
  res2 <- res2[,c("estimate","pval")]
  colnames(res2) <- c("estimate2","pval2")
  res <- merge(res1,res2,by.x=0,by.y=0)
  if(is.null(xlim)){
    xlim <- c(-1,1) * rep(max(abs(pretty(c(res$estimate1,res$estimate2)))),2)
  }
  if(is.null(ylim)){
    ylim <- xlim
  }
  plot.new()
  plot.window(xlim=xlim,ylim=ylim)
  # Fit:
  LM <- lm(estimate2~estimate1,dat=res)
  # Predict:
  newx <- seq(xlim[1],xlim[2],length.out=1000)
  prd <- predict(LM,newdata=data.frame(estimate1=newx),interval=c("confidence"),level=0.95,type="response")
  # 95%
  polygon(c(rev(newx), newx), c(rev(prd[ ,3]), prd[ ,2]), col = 'grey80', border = NA)
  lines(newx,prd[,2],col="grey60",lty=3,lwd=2)
  lines(newx,prd[,3],col="grey60",lty=3,lwd=2)
  # Trend:
  abline(LM,col="black",lwd=3,lty=2)
  points(newx,newx,lwd=1,lty=3,type="l")
  points(rep(0,2),ylim,lwd=1,lty=3,type="l")
  points(xlim,rep(0,2),lwd=1,lty=3,type="l")
  if(ptype!="norm"){
    if(ptype=="fdr"){
      res$pval1 <- p.adjust(res$pval1,"fdr")
      res$pval2 <- p.adjust(res$pval2,"fdr")
    }
    if(ptype=="bonf"){
      res$pval1 <- p.adjust(res$pval1,"bonferroni")
      res$pval2 <- p.adjust(res$pval2,"bonferroni")
    }
  }
  # Points:
  IND <- list(both=which(res$pval1<=pcut & res$pval2<=pcut),first=which(res$pval1<=pcut & res$pval2>pcut),
              second=which(res$pval1>pcut & res$pval2<=pcut),none=which(res$pval1>pcut & res$pval2>pcut))
  for(i in 1:length(IND)){
      points(res$estimate1[IND[[i]]],res$estimate2[IND[[i]]],pch=pch[i],bg=col[i],cex=cex[i],lwd=2)
  }
  # Axes:
  axis(1)
  axis(2,las=1)
  # Title:
  if(!is.null(main)){
    title(main=main)
  }
  if(!is.null(xlab)){
    title(xlab=xlab)
  }
  if(!is.null(ylab)){
    title(ylab=ylab)
  }
  # Beta:
  text(xlim[1]+0.05*(xlim[2]-xlim[1]),ylim[2]-0.1*(ylim[2]-ylim[1]),
       substitute(paste(beta , "=" , X),list(X=formatC(coefficients(summary(LM))["estimate1","Estimate"],digits=3))),font=2,adj=c(0,1))
  # R2:
  text(xlim[1]+0.05*(xlim[2]-xlim[1]),ylim[2]-0.15*(ylim[2]-ylim[1]),
       substitute(paste(R^2 , "=" , X),list(X=formatC((cor.test(res$estimate1,res$estimate2)$estimate)^2,digits=3))),font=2,adj=c(0,1))
}

