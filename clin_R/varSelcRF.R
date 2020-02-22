#adapted from varSelRF for randomForest to cForest
#02-Jan-14 Yen Low
#
#Usage: RScriptPath=Sys.getenv("RScriptPath")
#source(paste(RScriptPath,"/scripts/R/varSelcRF.R",sep=""))
#######################################

RScriptPath=Sys.getenv("RScriptPath") #Rscriptpath is an environment variable, path="/mnt/hgfs/Dropbox"
source(paste(RScriptPath,"/scripts/R/utils.R",sep=""))

installnewpackage(c("caret","party"))
require(caret)
require(party)

varSelcRF<-function (xdata, Class, c.sd = 1, mtryFactor = 1, ntree = 500, 
					ntreeIterat = 500, vars.drop.num = NULL, vars.drop.frac = 0.2, 
					whole.range = FALSE, recompute.var.imp = FALSE, verbose = FALSE, 
					returnFirstForest = TRUE,
					subsetID=NULL,samplewt=NULL,conditionalimpt=FALSE) {
    if (!is.factor(Class)) 
        stop("Class should be a factor")
    if ((is.null(vars.drop.num) & is.null(vars.drop.frac)) | 
        (!is.null(vars.drop.num) & !is.null(vars.drop.frac))) 
        stop("One (and only one) of vars.drop.frac and vars.drop.num must be NULL and the other set")
    max.num.steps <- dim(xdata)[2]
    num.subjects <- dim(xdata)[1]
	if (is.null(subsetID)) subsetID=1:length(Class)
	if (is.null(colnames(xdata))) colnames(xdata)=paste("V",1:dim(xdata)[2],sep="")
	
    n.vars <- vars <- OOB.rf <- OOB.sd <- rep(NA, max.num.steps)
	df=as.data.frame(xdata)
	df$Class=Class
	
    oobError <- function(rf) {
		ypred=factor(predict(rf,OOB=T),levels=0:1)
		yobs=factor(Class[subsetID],levels=0:1)
        stats <- confusionMatrix(ypred,yobs,"1")
		#Calc error=1-balanced accuracy (better to use balanced acc for model selection)
		1-mean(stats$byClass[1:2])
		#Calc error=1-accuracy (danger of selecting trivial models)		
		#1-stats$overall[["Accuracy"]]
    }
	
    mtry <- floor(sqrt(ncol(xdata)) * mtryFactor)
    rf <- cforest(Class~.,data=df,subset=subsetID,control=cforest_unbiased(ntree=ntree,mtry=mtry))
    if (returnFirstForest) 
        FirstForest <- rf
    else FirstForest <- NULL
    m.iterated.ob.error <- m.initial.ob.error <- oobError(rf)
    sd.iterated.ob.error <- sd.initial.ob.error <- sqrt(m.iterated.ob.error * 
        (1 - m.iterated.ob.error) * (1/num.subjects))
		
    if (verbose) cat("\n","Initial OOB error: ",round(m.initial.ob.error,4)," +/- ",round(sd.initial.ob.error, 4))
	
    importances <- varimp(rf, conditional=conditionalimpt)
    selected.vars <- names(sort(importances, decreasing = TRUE))
    ordered.importances <- importances[selected.vars]
    initialImportances <- importances
    initialOrderedImportances <- ordered.importances
	
    j <- 1
    n.vars[j] <- dim(xdata)[2]
    vars[j] <- paste("`",colnames(xdata),"`",collapse = " + ",sep="")	
    OOB.rf[j] <- m.iterated.ob.error
    OOB.sd[j] <- sd.iterated.ob.error
    
	var.simplify <- TRUE
	while (var.simplify) {
        if (verbose) {
			cat("\n",length(selected.vars)," variables\n")
			#print(selected.vars)
			#print(vars[j])
			cat("OOB error = ",round(OOB.rf[j],4)," +/- ",round(OOB.sd[j],4),"\n")
            print("gc inside loop of varSelRF")
            print(gc())
        }
        else {
            gc()
        }
        last.rf <- rf
        last.vars <- selected.vars
        previous.m.error <- m.iterated.ob.error
        previous.sd.error <- sd.iterated.ob.error
        if (length(selected.vars) <= 2) {
            var.simplify <- FALSE
            break
        }
        if (recompute.var.imp & (j > 1)) {
            importances <- varimp(rf,conditional=conditionalimpt)
            tmp.order <- names(sort(importances, decreasing = TRUE))
            selected.vars <- tmp.order
            ordered.importances <- importances[tmp.order]
        }
        num.vars <- length(selected.vars)
        if (is.null(vars.drop.num)) vars.drop=round(num.vars * vars.drop.frac)
        else vars.drop <- vars.drop.num
        if (num.vars >= (vars.drop + 2)) {
            selected.vars <- selected.vars[1:(num.vars - vars.drop)]
            ordered.importances <- ordered.importances[1:(num.vars-vars.drop)]
        }
        else {
            selected.vars <- selected.vars[1:2]
            ordered.importances <- ordered.importances[1:2]
        }
        if (length(selected.vars) < 2){
            var.simplify <- FALSE
            break
        }
        mtry <- floor(sqrt(length(selected.vars)) * mtryFactor)
        if (mtry > length(selected.vars)) 
            mtry <- length(selected.vars)
			fml=as.formula(paste("Class ~ ", paste("`",selected.vars,"`",collapse=" + ",sep="")))
        rf <- cforest(fml,data=df,subset=subsetID,control=cforest_unbiased(ntree=ntreeIterat,mtry=mtry))
        m.iterated.ob.error <- oobError(rf)
        sd.iterated.ob.error <- sqrt(m.iterated.ob.error * (1 - 
            m.iterated.ob.error) * (1/num.subjects))

		j <- j + 1
        n.vars[j] <- length(selected.vars)
        vars[j] <- paste("`",selected.vars,"`",collapse = " + ",sep="")
        OOB.rf[j] <- m.iterated.ob.error
        OOB.sd[j] <- sd.iterated.ob.error
        if (!whole.range & 
			((m.iterated.ob.error > (m.initial.ob.error + c.sd * sd.initial.ob.error)) |
			(m.iterated.ob.error > (previous.m.error + c.sd * previous.sd.error)))) 
            var.simplify <- FALSE
    }

#If whole.range=TRUE continue dropping variables until a forest with only two variables is built, 
#and choose the best model from the complete series of models. 
#If FALSE, stop the iterations if current OOB error > initial OOB error + c.sd*OOB standard error
#or current OOB error > previous OOB error + c.sd*OOB standard error)
    if (!whole.range) {
        selected.vars <- last.vars
        out <- list(selec.history = data.frame(Number.Variables=n.vars,
												Vars.in.Forest=vars,
												OOB=OOB.rf,sd.OOB=OOB.sd)[1:j,],
			rf.model = last.rf, selected.vars = selected.vars, 
            selected.model = paste("`",selected.vars,"`",collapse = " + ",sep=""), 
            best.model.nvars = length(selected.vars), initialImportances = initialImportances, 
            initialOrderedImportances = initialOrderedImportances, 
            ntree = ntree, ntreeIterat = ntreeIterat, mtryFactor = mtryFactor, 
            firstForest = FirstForest)
        class(out) <- "varSelRF"
        return(out)
    }
    else {
        n.vars <- n.vars[1:j]
        vars <- vars[1:j]
        OOB.rf <- OOB.rf[1:j]
        OOB.sd <- OOB.sd[1:j]
        min.oob.ci <- min(OOB.rf) + c.sd * OOB.sd[which.min(OOB.rf)]
        best.pos <- which(OOB.rf <= min.oob.ci)[which.min(n.vars[which(OOB.rf <= min.oob.ci)])]
        selected.vars <- unlist(strsplit(vars[best.pos], " + ", fixed = TRUE))
        out <- list(selec.history = data.frame(Number.Variables = n.vars, 
            Vars.in.Forest = vars, OOB = OOB.rf, sd.OOB = OOB.sd), 
            rf.model = NA, selected.vars = selected.vars, 
			selected.model = paste("`",selected.vars,"`",collapse = " + ",sep=""),best.model.nvars = n.vars[best.pos], 
            initialImportances = initialImportances, initialOrderedImportances = initialOrderedImportances, 
            ntree = ntree, ntreeIterat = ntreeIterat, mtryFactor = mtryFactor, 
            firstForest = FirstForest)
        class(out) <- "varSelRF"
        return(out)
    }
}

