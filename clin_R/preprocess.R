# Functions for preprocessing data matrix:
# 	1. rmCol: 
# 	2. rmRow:
#	3. recode:
#	4. reorderTable
#
# 21-Jan-14 Yen Low
#
# To use:
# source("/mnt/hgfs/Dropbox/scripts/R/preprocess.R")
###############################################################################
#RScriptPath=Sys.getenv("RScriptPath")
#source(paste(RScriptPath,"/scripts/R/utils.R",sep=""))
#installnewpackage("car")
require("car")

#TODO
#recode
#http://rprogramming.net/recode-data-in-r/
#SchoolData$Grade<-recode(SchoolData$Grade,"5=6")# Recode grade 5 to grade 6
#SchoolData$Grade<-recode(SchoolData$Grade,"`Grade Five`=5") #text must be in ` ` or ' '
#SchoolData$Grade<-recode(SchoolData$Grade,"c(1,2,3,4,5)='Five or Less'")
#SchoolData$Grade<-recode(SchoolData$Grade,"5=6;6=7")

reorderTable<-function(data,colOrder,rowOrder=NULL,outfile="table.txt"){
	if(is.null(rowOrder)) rowOrder=1:nrow(data)
	reorderedData=data[rowOrder,colOrder]
	if(!is.null(outfile)){
		cat(colOrder,"\n",sep="\t",file=outfile)
		write.table(reorderedData,file=outfile,row.names=F,col.names=F,
					sep="\t",quote=F,append=T)
	}
	return(reorderedData)
}


#Preprocess XA files
#rmCol removes near-constant and highly correlated columns (i.e. variables, descriptors)
#rmRow removes rows (i.e. items, samples, compounds) identical or near-identical to each other
#
#code modified from Rajarshi Guha http://rguha.net/code/R/#idtest
#Inputs:	1. x=descriptor matrix (or data frame)
#			2. maxFracConstant=0.9 (optional; Default to 0.9)
#			3. minSD=0.001 (optional; Default to 0.001)
#			4. maxR2=0.9 (optional; Default to 0.9)
#			5. dist=0 (optional; Default to 0)
#29-Jul-12 Yen Low
##################################################################


rmCol<-function(x,maxFracConstant=0.9,minSD=0.001,maxR2=0.9,maxNA=0.3,preferNonCfrag=FALSE){
	if (maxFracConstant > 1 || maxFracConstant <= 0) 	stop(" 0 <= maxFracConstant < 1")
	if (minSD < 0) 	stop(" minSD >= 0")
	if (maxR2 > 1 || maxR2 <= 0) 	stop(" 0 <= maxR2 < 1")
	if (maxNA > 1 || maxNA <= 0) 	stop(" 0 <= maxNA < 1")
	if (!is.matrix(x) && !is.data.frame(x)) stop("Must supply a data.frame or matrix")
	
	nrow_ori=nrow(x)
	#remove columns with > maxNA portion of NA
	dropidx=apply(x, 2, function(x) { sum(is.na(x))/nrow_ori > maxNA })
	x=x[, !dropidx]
	cat("Drop variables with > maxNA ratio:",sum(dropidx),"\n")
	
	#remove constant and near-constant descriptors (FractionConstant > maxFracConstant)
	dropidx=apply(x, 2, function(x) { sum(x==0,na.rm=T)/sum(!is.na(x)) > maxFracConstant })
	x=x[, !dropidx]
	cat("Drop constant/near-constant variables:",sum(dropidx),"\n")
	
	#remove descriptors with low SD
	colSD=apply(x,2,sd,na.rm=T)
	keepColSD=(colSD>=minSD)
	x=x[,keepColSD]
	cat("Drop low-variance variables:",sum(!keepColSD),"\n")
	
	#remove highly correlated descriptors (R2 > maxR2)
	r2cut = sqrt(maxR2)
	keepCol=1:ncol(x)
	drop.idx=c()
	if(!is.null(x)){
		cormat=cor(x)
		cor.idx=which(abs(cormat)>r2cut,arr.ind=T)
		if(!is.null(dim(cor.idx))){
			#remove duplicates and identity diagonal
			coridx=matrix( cor.idx[cor.idx[,1] > cor.idx[,2]], ncol=2)
			freqmat=matrix(NA,ncol=ncol(coridx),nrow=nrow(coridx))
			if(nrow(coridx)>0){
				#set preferNonCfrag=TRUE if fragments with more non carbon atoms are preferred
				#more interpretable
				if(preferNonCfrag){
					NumofC=NumofnonC=fraglength=matrix(NA,ncol=2,nrow=nrow(coridx))
					for(i in 1:nrow(coridx)){
						print(colnames(x)[coridx[i,]])
						fraglength[i,]=floor((nchar(colnames(x)[coridx[i,]])+1)/2)
						NumofC[i,]=unlist(lapply(gregexpr("C",colnames(x)[coridx[i,]]),length))
					}
					NumofnonC=fraglength-NumofC
					#return the shorter fragment to be dropped
					#if equal length, return the fragment with fewer non-C atoms to be dropped
					ID1or2=apply(fraglength,1,	function(x) ifelse(x[1]!=x[2],which.min(x),apply(NumofnonC,1,which.min))) 
				}else{
					#remove the descriptor with the most correlations (helps to min descriptors)
					freqofcor=table(coridx)
					freqmat[,1]=apply(coridx,1,function(x) freqofcor[names(freqofcor)==x[1]])
					freqmat[,2]=apply(coridx,1,function(x) freqofcor[names(freqofcor)==x[2]])
					ID1or2=apply(freqmat,1,function(x) ifelse(x[1]!=x[2],which.min(x),1+round(runif(1))))
					#randomly remove either descriptor of the correlation pair
					#drop.idx=ifelse(runif(nrow(coridx)) > .5, coridx[,1], coridx [,2])
				}
				
				for(i in 1:nrow(coridx)) drop.idx[i]=coridx[i,ID1or2[i]]
				if(length(drop.idx)>0) keepCol=(1:ncol(x))[-unique(drop.idx)]
			}
		}
	}
	
	cat("Drop highly correlated variables:",length(unique(drop.idx)),"\n")
	return(x=x[,keepCol])
}


rmRow<-function(x,dist=0,maxNA=0.5,dup=T){
	if (maxNA > 1 || maxNA <= 0) 	stop(" 0 <= maxNA < 1")
	if (!is.matrix(x) && !is.data.frame(x)) stop("Must supply a data.frame or matrix")
	
	ncol_ori=ncol(x)
	#remove columns with > maxNA portion of NA
	dropidx=apply(x, 1, function(x) { sum(is.na(x))/ncol_ori >= maxNA })
	x=x[!dropidx,]
	cat("Drop rows with > maxNA ratio:",sum(dropidx),"\n")
	
	if(dup){
		#remove identical rows
		if(dist==0){
			dupID=duplicated(x)
			cat("Drop duplicated rows:",sum(dupID),"\n")
			x=x[!dupID,]
		}else if(dist>0){
			#remove near-identical rows (euclidean distance < dist)
			keepRow=1:nrow(x)
			drop.idx=c()
			if(!is.null(x)){
				distmat=as.matrix(dist(x,diag=T,upp=T))
				dist.idx=which(distmat<dist,arr.ind=T)
				if(!is.null(dist.idx)){
					#remove duplicates and identity diagonal
					distidx=matrix( dist.idx[dist.idx[,1] > dist.idx[,2]], ncol=2)
					freqmat=matrix(NA,ncol=ncol(distidx),nrow=nrow(distidx))
					if(nrow(distidx)>0){
						#remove the row with the most similar descriptor profile
						#remove row that is least frequently similar (helps to min number of rows)
						#for rows with equal frequencies, randomly remove either row 
						freqofdist=table(distidx)
						freqmat[,1]=apply(distidx,1,function(x) freqofdist[names(freqofdist)==x[1]])
						freqmat[,2]=apply(distidx,1,function(x) freqofdist[names(freqofdist)==x[2]])
						ID1or2=apply(freqmat,1,function(x) ifelse(x[1]!=x[2],which.min(x),1+round(runif(1))))
					}
					for(i in 1:nrow(distidx)) drop.idx[i]=distidx[i,ID1or2[i]]
					if(length(drop.idx)>0) keepRow=(1:nrow(x))[-unique(drop.idx)]
				}
			}
			cat("Drop duplicated rows:",length(unique(drop.idx)),"\n")
			x=x[keepRow,]
		}else stop("dist>=0") #invalid distance threshold
	}
	return(x)
}

