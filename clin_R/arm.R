# Functions for:
# 	1. pruning itemsets (get k supersets meeting minimum lift and confidence
# 	2. aggregating drugs by ATC class (2nd or 3rd level recommended)
# input data: user_yenlow.mirtazapine_labGlucoseBeforeAfterMir
# 
# 21-Jan-14 Yen Low
#
# To use:
# source("/mnt/hgfs/Dropbox/scripts/R/arm.R")
###############################################################################

#RScriptPath=Sys.getenv("RScriptPath")
source(paste(RScriptPath,"/R/utils.R",sep=""))

installnewpackage(c("arules","arulesViz"))
require(arules)
require(arulesViz)

#Prune itemsets (get k supersets meeting minimum lift and confidence)
topKRules<-function(rules,metric="lift",k,minlift=1,minconf=0.5,pruned=TRUE,verbose=FALSE){
	if(is.numeric(k)) if(k>nrow(rules@quality)) k=nrow(rules@quality)
	if(nrow(rules@quality)!=0){
		if(is.numeric(k)){ 
			if(k>0){
				topID=order(rules@quality[,metric],decreasing=T)[1:k]
				rules_top=subset(rules[topID],subset=lift > minlift & confidence > minconf)
			}
		}else if(minlift>=1){
			rules_top=subset(rules,subset=lift > minlift & confidence > minconf)
		}
	}
	if(exists("rules_top")){
		if(pruned==TRUE & length(rules_top)>0){
			#find which pair of rules are subsets
			subset.matrix=is.subset(rules_top)
#keep only upp triangle
			subset.matrix[lower.tri(subset.matrix,diag=T)]=NA
#find and remove redundant rules
			redundant=colSums(subset.matrix, na.rm=T) >= 1
			which(redundant)
			rules_pruned=rules_top[!redundant]
			rules_top=rules_pruned
		}
		if(verbose==TRUE) print(inspect(rules_top))
		return(rules_top)
	}else{
		stop("Error: No rules found")
	}
}


#############
list2adjmat<-function(rules,weighted=TRUE){
	listOfItems=LIST(rules@lhs,decode=F)
	itemNames=sort(unique(unlist(listOfItems)))
	lens=unlist(lapply(listOfItems,length))
	g=list()
	gtotal=data.frame(from=NA,to=NA)
	for(i in 1:length(lens)){
		g[[i]]=graph.full(lens[i],loops=F)
		V(g[[i]])$name=listOfItems[[i]]
		gtotal=rbind(gtotal,get.data.frame(g[[i]]))
	}
	#gtotal=graph.union(unlist(g,recursive=F))
	gtotal=gtotal[-1,]
	adjmat=table(factor(gtotal$from, levels=itemNames),factor(gtotal$to, levels=itemNames))
	if(weighted==FALSE) adjmat[adjmat>=1]=1
	adjmatNames=dimnames(rules@lhs)[[2]][as.numeric(dimnames(adjmat)[[1]])]
	adjmatNames[adjmatNames=="hormoneTNBC"]="TNBC"
	adjmatNames[adjmatNames=="AFFILIATIONonly SU"]="Stanford"
	adjmatNames[adjmatNames=="histotypemalignant"]="malignant"
	dimnames(adjmat)=list(adjmatNames,adjmatNames)
	return(adjmat)
}

#########legacy code
#lift=interestMeasure(basket_rules,c("lift"),trans)
#conf=interestMeasure(basket_rules,c("allConfidence"),trans)
#measures=cbind(quality(basket_rules),lift,conf)
#topID=order(lift,decreasing=T)[1:10]
#cbind(inspect(basket_rules[topID]),measures[topID,-1])