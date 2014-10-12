# Functions for self-controlled study design
# Compare labs before and after drug use (self-control design; within the same patient)
# set parameters for lab component, washout, blackoutDays, minExposureTime, minLabs
# input data: user_yenlow.mirtazapine_labGlucoseBeforeAfterMir
# 
# 18-Nov-13 Yen Low
#
# To use:
# source("/mnt/hgfs/Dropbox/scripts/R/compareLabsBeforeAfter.R")
###############################################################################

#function to compare labs by blackout period and component type
compareLabsByPatient<-function(	labdata,comp,baselineDays=0,
								washout=7,blackoutDays=3,minExposureTime=30,labTime=30,labWindow=7,
								alt,var.eq=F,minLabs=30,savediff=F,pre="mean"){
	lab1Component=subset(labdata,component==comp & indexTime>=washout & PairedLabBeforeAndAfter==1 & exposureTime>=minExposureTime,
						select=c("pid","indexTime","lastExposureTime","exposureTime","timeoffset","component","ord_num","afterDrug"))
	n=nrow(lab1Component)
	print(n)
	if(nrow(lab1Component)>=minLabs){
		lab1Component$afterDrug=NA
		lab1Component$afterDrug[lab1Component$timeoffset>lab1Component$indexTime+labTime-labWindow
							&	lab1Component$timeoffset<lab1Component$indexTime+labTime+labWindow]=1
		lab1Component$afterDrug[lab1Component$timeoffset<lab1Component$indexTime+blackoutDays]=NA
		lab1Component$afterDrug[lab1Component$timeoffset<lab1Component$indexTime
							& 	lab1Component$timeoffset>lab1Component$indexTime-baselineDays]=0
		nlid=sum(lab1Component$afterDrug==1,na.rm=T)
		cat("nlids post-drug: ",nlid, "\n")
		if(nlid>=1){
			if(pre=="mean"){
			#if mean_pre - mean_post
				labs=aggregate(ord_num~afterDrug+pid,lab1Component,mean,na.rm=T) #will be paired before and after
			}else if(pre=="last"){
			#get mean_post
				labs=aggregate(ord_num~pid,lab1Component,mean,subset=lab1Component$afterDrug==1,na.rm=T)				
				labs$afterDrug=1
				#get last_pre
				temp=lab1Component[!is.na(lab1Component$afterDrug) & lab1Component$afterDrug==0,c("pid","ord_num","timeoffset")]
				temp=temp[order(temp$pid,-temp$timeoffset),] #sort by pid and then latest timeoffset
				labs_pre=temp[!duplicated(temp[,"pid"]),c("pid","ord_num")]
				labs_pre$afterDrug=0
			
				labs=rbind(labs_pre,labs)
			}else{
				stop("Specify aggregation function for labs before drug index date. Pre='mean' or 'last'.")
			}
			
	
			#lab1Component=subset(labBeforeAfterMir,component==3009 & indexTime>=7 & PairedLabBeforeAndAfter==1 & exposureTime>=0,
			#	select=c("pid","indexTime","lastExposureTime","exposureTime","timeoffset","component","ord_num","afterDrug"))
			labsxtab=reshape(labs,v.names="ord_num",idvar="pid",timevar="afterDrug",dir="wide") 
			#ensure labs exist in desired time point
			if(ncol(labsxtab)==3 & sum(!is.na(rowSums(labsxtab)))>=minLabs){
				diff=labsxtab[,"ord_num.1"]-labsxtab[,"ord_num.0"]
				labsxtab=cbind(labsxtab,diff)
				ttest=try(t.test(labsxtab[,"ord_num.1"],labsxtab[,"ord_num.0"],
								alternative=alt,paired=T,var.equal=var.eq,na.action=na.exclude))
				wilcoxtest=try(wilcox.test(labsxtab[,"ord_num.1"],labsxtab[,"ord_num.0"],
								alternative=alt,paired=T))
				if(class(ttest)!="try-error"& class(wilcoxtest)!="try-error"){
					print(ttest)
					print(wilcoxtest)
				}else{ #try error
					print("Too many pids with NA labs. Insufficient values for ttest!")
					ttest=NULL
					wilcoxtest=NULL
				}
			}else{   #not enough columns or rows for tests
				cat("No lab readings before or after indexTime in desired labTimePoint or <",minLabs,"patients for ttest!\n")
				ttest=NULL
				wilcoxtest=NULL
			}
		}else{  #not enough rows
			cat(n,"patients with labs after subsetting. Insufficient for comparison!\n")
			ttest=NULL
			wilcoxtest=NULL
			labsxtab=NULL
		}
		if(savediff) list(n=n,ttest=ttest,wilcoxtest=wilcoxtest,lab=labsxtab) else list(n=n,ttest=ttest,wilcoxtest=wilcoxtest)
	}else{
		cat("Zero rows after subsetting!\n")
		list(n=0,ttest=NULL,wilcoxtest=NULL)
	}
}

##########
#compare labs by blackout period and component type and exposure time; plots diff by blackout time
#calls compareLabsByPatient
#set baselineDays to some large number to go far back in time
compareLabsByLabTimes<-function(labdata,components,baselineDays=9999,
											washout=7,blackoutDays=3,minExposureTime=NULL,
											labTimePoints,labWindow=7,
											alt="less",var.eq=T,minLabs=30,savediff=F,pre="mean"){
	results=resultsByDay=list()
	for(i in 1:length(components)){
		for(j in 1:length(labTimePoints)){
			if(is.null(minExposureTime)){ #use labTimes as minExposureTime#use labTimes as minExposureTime#use labTimes as minExposureTime#use labTimes as minExposureTime#use labTimes as minExposureTime#use labTimes as minExposureTime
				resultsByDay[[j]]=compareLabsByPatient(labdata,components[i],baselineDays=baselineDays,
						washout=washout,blackoutDays,minExposureTime=labTimePoints[j]-labWindow,
						labTime=labTimePoints[j],labWindow=labWindow,
						alt=alt,var.eq=var.eq,minLabs=minLabs,savediff=savediff,pre=pre)
			}else if(minExposureTime>=0) { 
				resultsByDay[[j]]=compareLabsByPatient(labdata,components[i],baselineDays=baselineDays,
						washout=washout,blackoutDays,minExposureTime=minExposureTime,
						labTime=labTimePoints[j],labWindow=labWindow,
						alt=alt,var.eq=var.eq,minLabs=minLabs,savediff=savediff,pre=pre)
			}else{
				stop("minExposureTime must be a non-negative number or NULL (for exposure times up to labTime)")
			}
		}
		results[[i]]=resultsByDay
	}
	
#put relevant results in a matrix
	innerlength=length(results[[1]])
	mat=matrix(NA,nrow=length(results)*innerlength,ncol=14)
	for(p in 1:length(results)){
		for(q in 1:innerlength){
			if(!is.null(results[[p]][[q]][["lab"]])){
				mean_pre=mean(results[[p]][[q]][["lab"]]$ord_num.0,na.rm=T)
				mean_post=mean(results[[p]][[q]][["lab"]]$ord_num.1,na.rm=T)
				diffInMeans=mean_post - mean_pre				
				sd_pre=sd(results[[p]][[q]][["lab"]]$ord_num.0,na.rm=T)
				sd_post=sd(results[[p]][[q]][["lab"]]$ord_num.1,na.rm=T)
				if("diff" %in% colnames(results[[p]][[q]][["lab"]])){
					sd_meandiff=sd(results[[p]][[q]][["lab"]]$diff,na.rm=T)
					mediandiff=median(results[[p]][[q]][["lab"]]$diff,na.rm=T)
				}else{
					mediandiff=NA
					sd_meandiff=NA
				}
			}else{
				mean_pre=NA
				mean_post=NA
				diffInMeans=NA
				mediandiff=NA
				sd_pre=NA
				sd_post=NA
				sd_meandiff=NA
			}
			if (!is.null(results[[p]][[q]][["ttest"]])){
				df=unlist(results[[p]][[q]][["ttest"]]["parameter"])+1
				meandiff=unlist(results[[p]][[q]][["ttest"]]["estimate"])
				pval_t=unlist(results[[p]][[q]][["ttest"]]["p.value"])
				pval_wilcox=unlist(results[[p]][[q]][["wilcoxtest"]]["p.value"])
			}else{
				df=NA
				meandiff=NA
				pval_t=NA
				pval_wilcox=NA
			}
			mat[(p-1)*innerlength+q,1]=components[p]
			mat[(p-1)*innerlength+q,2]=labTimePoints[q]
			mat[(p-1)*innerlength+q,3]=results[[p]][[q]][["n"]]
			mat[(p-1)*innerlength+q,4]=df
			mat[(p-1)*innerlength+q,5]=mean_pre
			mat[(p-1)*innerlength+q,6]=mean_post
			mat[(p-1)*innerlength+q,7]=diffInMeans			
			mat[(p-1)*innerlength+q,8]=meandiff
			mat[(p-1)*innerlength+q,9]=mediandiff
			mat[(p-1)*innerlength+q,10]=sd_pre
			mat[(p-1)*innerlength+q,11]=sd_post
			mat[(p-1)*innerlength+q,12]=sd_meandiff
			mat[(p-1)*innerlength+q,13]=pval_t
			mat[(p-1)*innerlength+q,14]=pval_wilcox
		}
	}
	colnames(mat)=c("component","labTime","nlid","npid",
					"mean_pre","mean_post","diffInMeans","meanPairedDiff","medianPairedDiff",
					"sd_pre","sd_post","sd_meandiff",
					"pval_t","pval_wilcox")
	
	list(results=results,summarymat=mat)
}


###########################
#get patient Info
getPid=function(x){
	pid=list()
	counter=0
	for(p in 1:length(x$results)) for(q in 1:length(x$results[[p]])){
			counter=counter+1
			temp=x$results[[p]][[q]]$lab$pid
			if(is.null(temp)) temp=NA
			pid[[counter]]=temp
		}
	return(pid)
}


#get patient characteristics summary
getPatientSummary<-function(pid,patientInfo,saveData=F){
	piddem=patientInfo[patientInfo$pid %in% pid,]
	piddem$exposureTime=piddem$lastExposureTime-piddem$indexTime
	contVar=c(	colMeans(piddem[,c(	"indexTime","lastExposureTime",
									"exposureTime","ageIndexTime")],na.rm=T),
				apply(piddem[,c(	"indexTime","lastExposureTime",
									"exposureTime","ageIndexTime")],2,sd,na.rm=T))
	names(contVar)=c(	"mean_indexTime","mean_lastExposureTime",
						"mean_exposureTime","mean_ageIndexTime",
						"sd_indexTime","sd_lastExposureTime",
						"sd_exposureTime","sd_ageIndexTime")
	f_male=sum(piddem$gender=="MALE")/nrow(piddem)
	#f_race=prop.table(table(piddem$race))
	names(f_male)="male"
	if(!saveData) piddem=NULL
	temp=c(contVar,f_male)
	list(summary=list(temp,pidData=piddem))
}

