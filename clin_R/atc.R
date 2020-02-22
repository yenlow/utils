# Functions for:
# 	1. mapping drugs (RxNORM) to ATC
# 	2. aggregating drugs by ATC class (2nd or 3rd level recommended)
# input data: user_yenlow.mirtazapine_labGlucoseBeforeAfterMir
# 
# 07-Jan-14 Yen Low
#
# To use:
# source("/mnt/hgfs/Dropbox/scripts/R/atc.R")
###############################################################################

#RScriptPath=Sys.getenv("RScriptPath")
#source(paste(RScriptPath,"/scripts/R/utils.R",sep="")) #installnewpackage

#installnewpackage("RMySQL")
#library(RMySQL)

#######connection info to database (to get RxNORM, ATC from RxNORM2013 or terminology3)
#username="yenlow"
#passwd="yenlow81"
#host="ncbolabs-db1.stanford.edu"
#
#drv = dbDriver("MySQL")
##get RxNORM table
#con=dbConnect(drv,dbname="user_yenlow",user=username,host=host,password=passwd)
#rxnorm=dbGetQuery(con,"SELECT code,str FROM umls2011ab.MRCONSO where SAB='RXNORM' and tty in ('IN','PIN');")
#dbDisconnect(con)
#
##get mapping key between RxNORM and ATC
#con=dbConnect(drv,dbname="user_yenlow",user=username,host=host,password=passwd)
#atc=dbGetQuery(con,"SELECT rxcui, code FROM rxnorm.RXNCONSO where SAB='ATC' and tty in ('IN','PT');")
#dbDisconnect(con)
#
##get ATC names
#con=dbConnect(drv,dbname="user_yenlow",user=username,host=host,password=passwd)
#atcNames=dbGetQuery(con,"SELECT * FROM terminology3.atc2str;")
#dbDisconnect(con)
#colnames(atcNames)=c("atc","drugclass")
#atcNames$drugclass=tolower(atcNames$drugclass)
#
#rxnormAtcTable=merge(atc,rxnorm,by=1,sort=F,all=T)
#rxnormAtcTable$name=tolower(rxnormAtcTable$name)
#colnames(rxnormAtcTable)=c("rxcui","atc","name")
#save(rxnormAtcTable,atcNames,file="/mnt/hgfs/Dropbox/scripts/R/rxnorm2atc.RData")
load(paste(RScriptPath,"/R/rxnorm2atc.RData",sep=""))


##########function converting drug code 1 (e.g. rxcui) to drug code 2 (e.g. atc or name)
#examples
#drugs=c(44,7052,18,161,2670) #in rxcui
#x2y(drugs,x="rxcui",y="name")
#x2y(drugs,x="rxcui",y="atc")
#x2y(drugs,x="rxcui",y=c("atc","name"))

x2y<-function(drugs, x="rxcui",y="name"){
	return(unique(rxnormAtcTable[rxnormAtcTable[,x] %in% drugs,c(x,y)]))
}

##########function converting drug name to ATC
#or reverse: ATC to drug name
#(uses ATC table instead of generic rxnormAtcTable)
drug2atc<-function(drugs){
	return(atcNames[atcNames$drugclass %in% tolower(drugs),])
}

#ATC to drug or drug class
atc2drugclass<-function(atcs){
	return(atcNames[atcNames$atc %in% atcs,])
}


#############function for retrieving ATC level of drugs
#Example:
#getATClevel(drugs,x="rxcui",atclevel=2,names=T)

getATClevel<-function(drugs,x="rxcui",atclevel=2,names=T){
	if(!is.logical(names)) stop("Error: Names must be TRUE or FALSE") else if(names) atc=x2y(drugs,x=x,y=c("atc","name")) else atc=x2y(drugs,x=x,y="atc")
	if		(atclevel==1){
		atc$atcL=substr(atc$atc,1,1)
	}else if(atclevel==2){
		atc$atcL=substr(atc$atc,1,3)
	}else if(atclevel==3){
		atc$atcL=substr(atc$atc,1,4)
	}else if(atclevel==4){
		atc$atcL=substr(atc$atc,1,5)
	}else{
		stop("Error: atclevel must be between 1 to 4")
	}
	results=merge(atc,atcNames,by.x="atcL",by.y="atc")
	return(results)
}



#############function to aggregate indiv drugs into ATC class
#Example:
#drugs=unique(rxnormAtcTable$rxcui[1:1000])[1:100]
#drugmat=matrix(round(runif(1000)),10,100,dimnames=list(1:10,drugs))
#temp=aggATC(drugmat,x="rxcui",atclevel=2,atclabel=NULL)

aggATC<-function(drugmat,x="rxcui",atclevel=2,atclabel=NULL,pattern=NULL){
	#get ATC level
	atctable=getATClevel(colnames(drugmat),x=x,atclevel=atclevel,names=F)
	#get unique ATC level with >1 drug (no point mapping 1 drug to ATC)
	atcL=names(which(table(atctable$atcL)>1))
	if(!is.null(pattern)) atcL=atcL[grep(pattern,atcL)]
	#create blank output table
	atcmat=matrix(NA,nrow=nrow(drugmat),ncol=length(atcL),dimnames=list(rownames(drugmat),atcL))
	for(i in 1:length(atcL)){
		colID=colnames(drugmat)%in% atctable[atctable$atcL==atcL[i],x]
		atcmat[,i]=as.numeric(rowSums(drugmat[,colID],na.rm=T)>0)
	}
	if(!is.null(atclabel)){
		if(atclabel==TRUE){
			colnames(atcmat)=atcNames[atcNames$atc %in% atcL,"drugclass"]
		}else if(length(atclabel)==ncol(atcmat)){
			colnames(atcmat)=atclabel
		}else{
			stop("Error: atclabel must be NULL (for ATC level codes, default), TRUE (for drugclass), a vector of atc levels (check length)")
		}
	}
	return(atcmat)
}


