# Process glmnet objects (calculate CCR, tune lambda and alpha)
# introduced multi-core parallel option (foreach) - call doParallel/Snow outside of this script
# removed jack-knifing and hnet
# 
# 01-Jun-14 Yen Low
###############################################################################
#Call script by:
#source("C:/Users/yenlow/Documents/Dropbox/MML_PIT/SOFTWARE/R/logit.R")

#RScriptPath=Sys.getenv("RScriptPath")
#RScriptPath="~/scripts"
source(paste(RScriptPath,"/R/utils.R",sep="")) #installnewpackage
source(paste(RScriptPath,"/R/sampling.R",sep=""))
installnewpackage(c("glmnet","foreach","reshape"))
require(glmnet)
require(foreach)
require(reshape)
#library(hierNet)

#calculate CCR from linear model object
getCCRfrlm<-function(mod,dsact){
  xtab=table(round(mod$fitted.values),dsact)
  cm=confusionMatrix(xtab,positive="1")
  CCR=mean(cm$byClass[c("Sensitivity","Specificity")])
  print(CCR)
  return(cm)
}

#hierNet
#tunes lambda in hierNet.logistic
optlambda<-function(x,y,nfold=5,nlam=10,strong=F,stand.main=stand.main){
  #tune lambda (fit to all possible lambda values; default 20 possible values
  hnettune=hierNet.logistic.path(x,y,nlam=nlam,strong=strong,diagonal=F,
                                 stand.main=stand.main,stand.int=F)
  #tune by n-fold CV
  hnetcv=hierNet.cv(hnettune,x,y,nfold=nfold,trace=0)
  return(hnetcv$lamhat.1se)
}

#get coefficients from hierNet object (across folds)
coeffAcrossFolds_hnet<-function(hnet){
  coeffarray=array(NA,dim=c(dim(hnet[[1]]$th),5))
  for(d in 1:length(hnet)){
    coeffarray[,,d]=hnet[[d]]$th
    #fill diagonal with main effects
    diag(coeffarray[,,d])=hnetmod[[d]]$bp-hnetmod[[d]]$bn
  }
  varnames=names(hnet[[d]]$mx)
  dimnames(coeffarray)=list(varnames,varnames,paste("fold",1:length(hnet),sep=""))
  
  #add upp tri to lower tri
  #note equation is set up to split interaction between 2 halve so add them back
  temp=coeffarray
  for(k in 1:dim(coeffarray)[3]){
    for(i in 1:nrow(coeffarray)){
      for(j in 1:(i-1)){
        temp[i,j,k]=sum(coeffarray[i,j,k],coeffarray[j,i,k])
        temp[j,i,k]=temp[i,j,k] #make symmetrical
      }
    }
  }
  
  return(temp)
}


######################################
#glmnet
#select optimal lambda and alpha in glmnet
plot.optimallambda<-function(cvfit,alphas,lambda="lambda.1se"){
  mincoeff=lapply(cvfit,coef,s=lambda)
  mincoeffmat=mincoeff[[1]][,"1"]  #fill up first row of mincoeffmat
  for(b in 2:length(mincoeff)) mincoeffmat=cbind(mincoeffmat,mincoeff[[b]][,"1"])
  colnames(mincoeffmat)=alphas
  plot(alphas,c(range(mincoeffmat,na.rm=T),0,0,0),pch="",ylab="coeffient value")
  abline(h=0,col="red")
  for(c in 1:nrow(mincoeffmat)) lines(alphas,mincoeffmat[c,],col=c)
  cat("\n","fold",i,"Number of zero coefficients\n")
  apply(mincoeffmat,2,function(x) print(sum(x==0)))
}

#get coefficients from lm object (across folds)
coeffAcrossFolds<-function(glmnet){
  coeff=lapply(glmnet,function(x) as.matrix(coef(x)))
  coeffmat=coeff[[1]]
  for(d in 2:length(coeff)) coeffmat=cbind(coeffmat,coeff[[d]])
  colnames(coeffmat)=paste("fold",1:length(glmnet),sep="")
  return(coeffmat)
}

#get optimal beta coefficients cv.glmnet at s=lambda.1se
coeffAtlambda<-function(cvfit,lambda="lambda.1se"){
  coeff_interpretation=coef(cvfit,s=lambda)
  coeff_interpretmat=as.matrix(coeff_interpretation)
  coeff_interpretmat[-(coeff_interpretation@i+1),]=NA #replace . with NA
  coeff=as.vector(coeff_interpretmat)
  names(coeff)=rownames(coeff_interpretmat)
  return(coeff)
}


#form sq matrix from vector of cross pair terms
pairvec2sqmat<-function(pairvec,varorder,upptri=T,outputfile="sqmat.txt",backquote=T) {
  temp=as.data.frame(pairvec)
  colnames(temp)="value"
  if(!backquote){ #add backquotes around every variable
    newnames=paste("`",names(pairvec),"`",sep="")
    names(pairvec)=gsub(":","`:`",newnames) #replace : with `:`
  }	
  temp$varname=sub("`:`","`split`",names(pairvec))
  vars=strsplit(temp$varname,"split")
  temp$rowvar=unlist(lapply(vars,`[`,1))
  temp$colvar=unlist(lapply(vars,`[`,2))
  temp$colvar[is.na(temp$colvar)]=temp$rowvar[is.na(temp$colvar)]
  temp2=as.matrix(cast(temp,colvar~rowvar,value="value"))
  #reduced matrix(low.tri) with p-1 dimensions (no intercept)
  #replace `for formula with ""(nothing)
  rownames(temp2)=gsub("`","",rownames(temp2))
  colnames(temp2)=gsub("`","",colnames(temp2))
  #reorder by correlated variables together
  rowID=match(varorder,rownames(temp2))
  colID=match(varorder,colnames(temp2))
  sqmat=temp2[rowID,colID]
  #replace blank column and row with appropriate variable name
  rownames(sqmat)[is.na(rownames(sqmat))]=varorder[is.na(rowID)]
  colnames(sqmat)[is.na(colnames(sqmat))]=varorder[is.na(colID)]
  #fill lower tri with some values in upp triangle
  #check that upp triangle is not all NA, else transpose
  if(sum(is.na(sqmat[upper.tri(sqmat)]))==sum(upper.tri(sqmat))) sqmat=t(sqmat)
  for(p in 1:(nrow(sqmat)-1)){
    for(q in (p+1):ncol(sqmat)){
      #make symmetrical by replacing NA with reflected values
      if(!is.na(sqmat[p,q]) & is.na(sqmat[q,p])){
        sqmat[q,p]=sqmat[p,q]
      }else if(is.na(sqmat[p,q]) & !is.na(sqmat[q,p])){
        sqmat[p,q]=sqmat[q,p]
      }
    }
  }
  if(upptri==FALSE) sqmat[upper.tri(sqmat)]=NA #fill upp tri?
  if(!is.null(outputfile)) write.table(sqmat,file=outputfile,sep="\t",quote=F,col.names=NA)
  return(sqmat)
}


#####################
#generate CI by bootstrapping
#use median alpha and lambda values (specific global best lambda or use cv.glmnet to choose from 20 values)
genCI<-function(xmat,yvec,ntrials=100,family="binomial",lambda,alpha=1,penalty.factor=1){
  #loop through 1 to ntrials
  modnet=foreach(j=1:ntrials,.packages="glmnet") %dopar% {
    if(is.vector(yvec)){  #for most models
      #print("yvec is a vector for most regression except Cox")
      mod_ID=sample(1:length(yvec),replace=T)
      yvec_subset=yvec[mod_ID]
    }else if(ncol(yvec)==2){  #for cox model where yvec is a matrix or Surv obj (2 columns)
      #print("yvec has 2 columns (time, status) for Cox regression")
      mod_ID=sample(1:nrow(yvec),replace=T)
      yvec_subset=yvec[mod_ID,]
    }else stop("yvec must be a vector OR matrix/Surv obj of 2 columns (time, status) for Cox regression")
    
    if(!is.null(lambda)){
      glmnet(xmat[mod_ID,],yvec_subset,thresh=1e-5,
             alpha=alpha,penalty.factor=penalty.factor,lambda=lambda,
             family=family,standardize=F)
    }else{
      cv.glmnet(xmat[mod_ID,],yvec_subset,thresh=1e-5,nfold=3,
                alpha=alpha,penalty.factor=penalty.factor,nlambda=20,
                family=family,standardize=F,parallel=T)
    }
  } #end of outer j foreach loop
  
  if(!is.null(lambda)) {
    temp=sapply(modnet,coef)
  }else{
    temp=sapply(modnet,function(x) coef(x,s=x$lambda.1se))      
  }
  coeffarray=temp[[1]]
  for(p in 2:length(temp)) coeffarray=cBind(coeffarray,temp[[p]])
  
  list(coeffarray=coeffarray,type="glmnet")
}


############################
#extract CI limits
setCL<-function(CIobj,CL=c(0.025,0.5,0.975),maxNA=0.5,verbose=FALSE){
  betaCIlim=t(apply(CIobj$coeffarray,1,quantile,probs=CL,na.rm=T,names=F))
  colnames(betaCIlim)=c("lowlim","median","upplim")
  beta_nonZero=!apply(betaCIlim,1,function(x) x[1]==0 & x[3]==0)
  fracNA=(rowSums(is.na(CIobj$coeffarray))/ncol(CIobj$coeffarray))<maxNA
  beta_nonZero=((beta_nonZero+fracNA)==2)
  #set main effects indices to TRUE
  #lowtridiag=nrow(CIobj$coeffarray)
  #nvar=floor(uniroot(function(x) x^2 + x -2*(lowtridiag-1),c(1,sqrt(lowtridiag*2)))$root)
  #beta_nonZero[1:(nvar+1)]=TRUE
  #if(verbose==TRUE) print(betaCIlim[beta_nonZero,])
  #   if(CIobj$type=="hnet"){
  #     betaCIlim=apply(CIobj$coeffarray,c(2,1),quantile,probs=CL,na.rm=T,names=F)
  # 		dimnames(betaCIlim)[[1]]=c("lowlim","median","upplim")
  # 		beta_nonZero=!apply(betaCIlim,c(2,3),function(x) x[1]==0 & x[3]==0)
  # 		#print(betaCIlim["lowlim",,]*beta_nonZero)
  # 	}
  
  return(list(betaCIlim=betaCIlim,beta_nonZero=beta_nonZero))
}


##############################
#generate beta table with 3 columns: 95%CI lower limit, beta estimate, 95%CI upper limit
genBetaCI<-function(mod,outfile_beta="beta",outfile_OR="OR"){
  #get beta coeff and SE
  #if glm model
  if("glm" %in% attributes(mod)$class){
    results=summary(mod)$coefficients
    betas=results[,1]
    se=results[,2]
  }else if(attributes(mod)$class=="ridgeLogistic"){
    #if ridge model (ridge package)
    results=summary(mod)$summaries$summary1$coefficients
    #unscale beta coeff and SE
    betas=results[-1,2]/mod$scales
    se=results[-1,3]/mod$scales
  }else{
    stop("mod must be a model generated by lm, glm, glmnet, or ridge (from ridge package)")
  }
  
  betatable=as.data.frame(cbind(betas-1.96*se, betas, betas+1.96*se))
  colnames(betatable)=c("lowlim","median","upplim")
  
  #get significant coeff
  if("glm" %in% attributes(mod)$class){
    sigEff=betatable[results[,4]<=0.05,]
  }else if(attributes(mod)$class=="ridgeLogistic"){
    sigEff=betatable[results[-1,5]<=0.05,]
  }
  
  ORtable=exp(sigEff)
  
  if(is.character(outfile_beta)){
    write.table(betatable,file=paste(outfile_beta,".txt",sep=""),sep="\t",quote=F,col.names=NA)
    write.table(betatable[order(betatable[,"median"],decreasing=T),],file=paste(outfile_beta,"_sorted.txt",sep=""),sep="\t",quote=F,col.names=NA)
    if(is.character(outfile_OR)){
      write.table(ORtable,file=paste(outfile_OR,"_sig.txt",sep=""),sep="\t",quote=F,col.names=NA)
    }else{
      print("Set outfile_OR to filenames to save significant OR estimates and CIs as .txt file. Set to NULL for no output.")
    }
  }else{
    print("Set outfile_beta to filenames to save beta estimates and CIs as .txt file. Set to NULL for no output.")
  }
  list(betatable=betatable,ORsig=ORtable)
}



forestplot<-function(d,orderoflevels=NULL,xlab="Odds Ratio",ylab="",xmin=0.25, xmax=3, xgap=0.25, sortByLabels=FALSE){
  require(ggplot2)
  if(!is.null(orderoflevels)){
    d$x=factor(rownames(d),levels=orderoflevels)
  }else{
    if(sortByLabels==FALSE) d$x=factor(rownames(d),levels=rownames(d)) else d$x=factor(rownames(d),levels=sort(rownames(d)))
  }
  p=ggplot(d, aes(x=x, y=median, ymin=lowlim, ymax=upplim)) + 
    geom_pointrange() + 
    geom_hline(aes(yintercept=1), lty=2) +
    coord_flip() +
    ylab(xlab) +
    xlab(ylab) + #switch because of the coord_flip() above
    theme_bw() +
    scale_y_log10(breaks=seq(xmin,xmax,by=xgap),limits=c(xmin,xmax))
  #+ scale_y_continuous(limits=c(0, 100))
  #+ scale_colour_manual(values = c("black","red"))
  #+ opts(  panel.border=theme_rect(colour="black",size=1),
  #         panel.grid.major = theme_blank(),
  # 	      legend.position="none")
  plot(p)
  return(p)
}
