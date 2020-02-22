# Functions for:
# 	1. calculating the Charlson Index based on 17 diseases (Deyo's adaptation)
# INPUTS: 2 columns (patient id and icd 9 code, numeric with decimal places)
# OUTPUTS: comorbmat=comorbidity matrix (rows as patients, 17 columns for 17 corbidities)
# and Charlson index of 17 comorbidities
#
# To use:
# source("https://www.dropbox.com/s/jl4by4wg5fc0wdd/comorbidities.R")
# OR source(paste(RScriptPath,"/scripts/R/comorbidities.R",sep=""))
#where RScriptPath=Sys.getenv("RScriptPath")
# 04-Feb-14 Yen Low
###############################################################################

###Calculate the Deyo index calculated from presence of 17 comorbidities
#"MI","CHF","PVD","CVD","DEMENTIA","COPD","RHEUM","PUD","MILD.LIVER","DM","DM.COMP","PLEGIA","RENAL","MALIGNANCY","SEVERE.LIVER", "METS", "HIV"
# INPUTS: 2 columns (patient id and icd 9 code, numeric with decimal places)
# OUTPUTS: comorbmat=comorbidity matrix (rows as patients, 17 columns for 17 corbidities)
# and Charlson index of 17 comorbidities
# wts=1 except for:
#2 each: Hemiplegia, moderate or severe kidney disease, diabetes with end organ damage, tumor, leukemia, lymphoma.
#3 each: Moderate or severe liver disease.
#6 each: Malignant tumor, metastasis, AIDS
#http://en.wikipedia.org/wiki/Comorbidity
#adapted from comorbidities R package (Don't use! Very slow due to bad indexing)
#TODO: adapt the other comorbidity indices (Elixhauser, AHRQ which modified the Elixhauser)
comorbidities.deyo <- function(code, verbose=TRUE) {
	
	#get unique patient IDs
	pid=unique(code[,1])
	
	#special handling for V codes (E codes are not used)
	pid_pvdtxt=as.character(unique(code[grep("^V43.4",code[,2]),1]))
	
	#create vectors of 17 comorbidities
	mi <- c(seq(from=410, to=410.9, by=0.01), 412)
	chf <- c(seq(from=428, to=428.9, by=0.01))
	pvd <- c(443.9, 441, 441.9, 785.4) #v code v43.4 not included in this list
	cvd <- c(seq(from=430, to=438, by=0.01))
	dementia <- c(seq(from=290, to=290.9, by=0.01))
	copd <- c(seq(from=490, to=496, by=0.01), seq(from=500, to=505, by=0.01), 506.4)
	rheum <- c(710, 710.1, 710.4, seq(from =714, to=714.2, by=0.01), 714.81, 725)
	pud <- c(seq(from=531, to=534.9, by=0.01))
	mild.liver <- c(571.2, 571.5, 571.6, seq(from=571.4, to=571.49, by=0.01))
	dm <- c(seq(from=250,to=250.3,by=0.01), 250.7)
	dm.comp <- c(seq(from=250.4, to=250.6, by=0.01)) #2 point items start here
	plegia <- c(344.1, seq(from=342, to=342.9, by=0.01))
	renal <- c(seq(from=582, to=582.9, by=0.01), seq(from=583, to=583.7, by=0.01), 585, 586,seq(from=588, to=588.9, by=0.01))
	malignancy <- c(seq(from=140, to=172.9, by=0.01), seq(from=174, to=195.8, by=0.01), seq(from=200, to=208.9, by=0.01))
	severe.liver <- c(seq(from=572.2, to=572.8, by=0.01),seq(from=456, to=456.21, by=0.01)) # 3 point item
	mets <- c(seq(from=196, to=199.1, by=0.01)) # 6 point items
	hiv <- c(seq(from=42, to=44.93, by=0.01))
	
	#create list of 17 deyo comorbidities
	deyo.list <- list(mi,chf,pvd,cvd,dementia,copd,rheum,pud,mild.liver,dm,dm.comp,plegia,renal,malignancy,severe.liver,mets,hiv)
	
	#create comorbidity matrix (rows as patients, 17 columns for 17 corbidities)
	outmat=matrix(0, nrow=length(pid), ncol=17,
			dimnames=list(pid,c("MI","CHF","PVD","CVD","DEMENTIA","COPD","RHEUM","PUD","MILD.LIVER","DM","DM.COMP","PLEGIA","RENAL","MALIGNANCY","SEVERE.LIVER", "METS", "HIV")))
	for (k in 1:length(deyo.list)) outmat[as.character(code[code[,2] %in% deyo.list[[k]],1]),k]=1
	if(length(pid_pvdtxt)>0) outmat[pid_pvdtxt,"PVD"]=1 #special handling if V43.4 was present
	
	#weighted sum of comorb
	wts=c(rep(1,10),rep(2,4),3,rep(6,2))
	deyoIndex=outmat %*% wts
	names(deyoIndex)=rownames(outmat)
	
	if(verbose==TRUE) print(colSums(outmat),na.rm=T)
	
	list(deyoIndex=deyoIndex,comorbmat=as.data.frame(outmat))
}


