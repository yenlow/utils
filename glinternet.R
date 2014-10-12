#######functions for glinternet
# 21-Jul-14 Yen Low

source(paste(RScriptPath,"/R/logit.R",sep="")) #installnewpackage

#to eliminate rows and/or columns coded 0
stripDummy<-function(x){
  #cat(nrow(x),ncol(x),"\n")  
  if(nrow(x)==2 & ncol(x)==2){
    #eliminate row and column
    ans=x[-1,-1]
  }else if(nrow(x)==2 & ncol(x)!=2){
    #eliminate dummy row
    ans=x[-1,]
  }else if(nrow(x)!=2 & ncol(x)==2){
    #eliminate dummy column
    ans=x[,-1]
  }else{
    ans=x
  }
  return(ans)
}

#function for generating suffixes (categorical codes)
getSuffix<-function(x,y){
  #if matrix
  if(length(x)==2){
    suffix=apply(expand.grid(0:(x[1]-1),0:(x[2]-1)),1,paste,collapse="_",sep="")
  }else if(length(x)==1 & x>1){
  #if vector
    if(is.null(y)){
      suffix=0:(x-1)
    }else if(sum(!grepl("x[0-9]",y))==0){
      suffix=as.numeric(gsub("x","",y))-1
    }else if(sum(!grepl("cat[0-9]+_",y))==0){
      suffix=as.numeric(gsub("cat[0-9]+_","",y))
    }else{
      suffix=0:(x-1)
    }
  }else{
    suffix=NULL
  }
  
  return(suffix)
}



###########################
#correct suffixes of coefficients
correctCoeffSuffix<-function(coeff){
  #get dimensions of coeffef
  nameCoeff=lapply(coeff,names)
  lengthCoeff=lapply(coeff,length)
  dimCoeff=lapply(coeff,dim)
  dimCoeff_final=mapply(function(x,y) if(!is.null(y)) x=y else x=x,lengthCoeff,dimCoeff)
  
  #name elements of list with appropriate suffix
  suffix=mapply(getSuffix,dimCoeff_final,nameCoeff)
  if(!is.list(suffix)) suffix=list(suffix)
  for(i in 1:length(coeff)) names(coeff[[i]])=suffix[[i]]
  return(coeff)
}


###########################
#get main and interaction effects from coeffList from glinternet
getMainIntEff<-function(glinmod,xnames,numCont){
  
  #extract coefficients from glinmod
  coeffList=coef(glinmod,lambdaIndex=2)
  
  #extract main effects (ensure continuous before categorical)
  maineff_cont=unlist(coeffList[[2]]$mainEffectsCoef$cont)  #main effects (continuous)
  #eliminate first level (coded 0) if binary
  maineff_cat=lapply(coeffList[[2]]$mainEffectsCoef$cat,function(x) if(length(x)==2) x[-1] else x) 
  maineff=c(maineff_cont,maineff_cat)
  #offset categorical main effects by numCont
  names(maineff)=xnames[c(coeffList[[2]]$mainEffects$cont,coeffList[[2]]$mainEffects$cat+numCont)]
  
  #extract interaction effects
  inteff_contcont=unlist(coeffList[[2]]$interactionsCoef$contcont)
  #eliminate first level (coded 0) if binary
  inteff_contcat=lapply(coeffList[[2]]$interactionsCoef$catcont,function(x) if(length(x)==2) x[-1] else x)
  #eliminate first row and first col (coded 0)
  inteff_catcat=lapply(coeffList[[2]]$interactionsCoef$catcat,stripDummy)
  inteff=c(inteff_contcont,inteff_contcat,inteff_catcat)
  #offset binary indices by 2
  contcontInd=coeffList[[2]]$interactions$contcont
  catcontInd=cbind(coeffList[[2]]$interactions$catcont[,1]+numCont,
                   coeffList[[2]]$interactions$catcont[,2])
  catcatInd=coeffList[[2]]$interactions$catcat+numCont
  names(inteff)=apply(matrix(xnames[rbind(contcontInd,catcontInd,catcatInd)],ncol=2,byrow=F),1,paste,collapse=":")

  maineffSuffixed=correctCoeffSuffix(maineff)
  inteffSuffixed=correctCoeffSuffix(inteff)

  list(maineff=maineffSuffixed,inteff=inteffSuffixed)
}


###################
#generate CI for glinternet coefficients
genCI_glin<-function(xmat,yvec,ntrials=100,family="binomial",lambda,numLevels,numToFind){
  #loop through 1 to ntrials
  glinmodList=foreach(j=1:ntrials,.packages="glinternet") %dopar% {
    mod_ID=sample(1:length(yvec),replace=T)
    yvec_subset=yvec[mod_ID]
    glinternet(xmat[mod_ID,],yvec_subset,numLevels=numLevels,family=family,
               lambda=lambda,tol=1e-3,maxIter=1000,numCores=1,numToFind=numToFind)

  } #end of outer j foreach loop
  return(glinmodList)
}


getBetaCIArray<-function(glinmodList,xnames,numCont){
  #keep post processin gof models outside of foreach loop for modular troubleshooting
  ##reprocess coeffList (list of results  per trial) into matrix
  #get first trial results to initiate coeffarray
  coeffarray=as.matrix(unlist(getMainIntEff(glinmodList[[1]],xnames,numCont=numCont)))
  #ensure maineff before inteff
  rownames(coeffarray)=gsub("maineff.","amxx0.",rownames(coeffarray)) 
  rownames(coeffarray)=gsub("inteff.","bixx1.",rownames(coeffarray))
  
  #fill coeffarray with other trial results
  for(p in 2:length(glinmodList)){
    coeff=unlist(getMainIntEff(glinmodList[[p]],xnames,numCont=numCont))
    #remove NA values
    vec=coeff[!grepl("NA$",names(coeff))]
    #ensure maineff before inteff
    temp=gsub("maineff.","amxx0.",names(vec)) 
    temp=gsub("inteff.","bixx1.",temp)
    
    coeffarray=merge(coeffarray,cbind(temp,vec),by.x=0,by.y=1,all=T,sort=T)
    rownames(coeffarray)=coeffarray[,"Row.names"]
    coeffarray=coeffarray[,-grep("Row.names",colnames(coeffarray))]
    colnames(coeffarray)=1:p
  }
  coeffarray=as.matrix(coeffarray)
  mode(coeffarray)="numeric"
  rownames(coeffarray)=gsub("amxx0.|bixx1.","",rownames(coeffarray))
  return(coeffarray)
}



###replace medians with beta estimates from glinmod
replaceMedianWithModelBeta<-function(betaTable,glinmod,xnames,numCont){
  betas=unlist(getMainIntEff(glinmod,xnames,numCont))
  names(betas)=gsub("^maineff.|^inteff.","",names(betas))
  temp=merge(betaTable$betaCIlim[betaTable$beta_nonZero,],as.matrix(betas),by=0,all.x=T,sort=F)
  temp$na=temp$median
  temp$na[!is.na(temp$V1)]=temp$V1[!is.na(temp$V1)]  #replace only if not NA
  beta_replaced=temp[,c("lowlim","na","upplim")]  #offset by 1 cos first col is rownames
  colnames(beta_replaced)=c("lowlim","median","upplim")
  rownames(beta_replaced)=temp[,"Row.names"]
  return(beta_replaced)
}



pairvec2sqmat_glin<-function(pairvec,varorder=NULL,flipvar=NULL,upptri=T,outputfile="sqmat.txt",backquote=F,returnLongForm=FALSE){
  
#   #convert x2 suffixes to .1 suffixes
#   names(pairvec)=gsub(".x6$",".5",names(pairvec))
#   names(pairvec)=gsub(".x5$",".4",names(pairvec))
#   names(pairvec)=gsub(".x4$",".3",names(pairvec))
#   names(pairvec)=gsub(".x3$",".2",names(pairvec))
#   if(!is.null(specialstring)){
#     specialID=grep(specialstring,names(pairvec))
#     names(pairvec)[specialID]=gsub(".x2$",".1",names(pairvec)[specialID])
#   }
#   names(pairvec)=gsub("\\.x2$","",names(pairvec))
  
  #flip special cases
  flippedvar=sapply(strsplit(flipvar,":"), function(x) paste(x[2:1],collapse=":",sep=""))
  for(i in 1:length(flipvar)) names(pairvec)=gsub(flipvar[i],flippedvar[i],names(pairvec))
  
  #get suffixes
  suffstring=gregexpr("\\.[0-9]*_*[0-9]$",names(pairvec))
  suffpos=sapply(suffstring,function(x) x[[1]])
  suffixes=mapply(function(x,y) if(y!=-1) strsplit(substring(x,y+1),"_") else NA,names(pairvec),suffpos)
  
  #split pairwise interaction into rowvar and col var
  vars=strsplit(gsub("\\.[0-9]*_*[0-9]$","",names(pairvec)),":")
  rowvar=mapply(function(x,y) if(!is.na(y[1])) paste(x[1],y[1],sep=".") else x[1],vars,suffixes)
  colvar=mapply(function(x,y) if(!is.na(x[2])){ if(!is.na(y[2])) paste(x[2],y[2],sep=".") else  x[2]} else NA,vars,suffixes)
  
  newintnames=mapply(function(x,y) if(!is.na(y)) paste(x,y,sep=":") else x, rowvar,colvar)
  names(pairvec)=newintnames
  
  #form beta square matrix
  if(is.null(varorder)) varorder=sort(unique(na.omit(c(rowvar,colvar))))
  betamat=pairvec2sqmat(pairvec=pairvec,varorder=varorder,upptri=upptri,outputfile=outputfile,backquote=backquote)
  if(returnLongForm==TRUE){
    longform=as.data.frame(cbind(rowvar,colvar,pairvec))
    return(longform)
  }else return(betamat)
}


##create model matrix for cluster robust CI adjustment
createModelMatrixOfVarWithLevels<-function(varWithLevels,data){
  temp=data[,varWithLevels[1:2]]
  for(i in 1:length(varWithLevels)){
    tab=table(data[,varWithLevels[i]])
    minmaxlevels=range(as.integer(names(tab)),na.rm=T)
    data[,varWithLevels[i]]=factor(data[,varWithLevels[i]],levels=minmaxlevels[1]:minmaxlevels[2])
    temp2=model.matrix(as.formula(paste("~",varWithLevels[i],"-1",sep="")),data)
    colnames(temp2)=paste(varWithLevels[i],minmaxlevels[1]:minmaxlevels[2],sep=".")
    temp=cbind(temp,temp2)
  }
  temp=temp[,-(1:2)]
  return(temp)
}