#Modified from varSelRF from varSelRF package (to pass more parameters into randomForest)
#02-Jan-14 Yen Low
#
#Usage: RScriptPath=Sys.getenv("RScriptPath")
#source(paste(RScriptPath,"/scripts/R/varSelRF.R",sep=""))

varSelRF<-function (xdata, Class, c.sd = 1, mtryFactor = 1, ntree = 5000, 
    ntreeIterat = 2000, vars.drop.num = NULL, vars.drop.frac = 0.2, 
    whole.range = TRUE, recompute.var.imp = FALSE, verbose = FALSE, 
    returnFirstForest = TRUE, fitted.rf = NULL, keep.forest = FALSE,
	sampsize=sampsize, strata=Class) 
{
    if (!is.factor(Class)) 
        stop("Class should be a factor")
    if ((is.null(vars.drop.num) & is.null(vars.drop.frac)) | 
        (!is.null(vars.drop.num) & !is.null(vars.drop.frac))) 
        stop("One (and only one) of vars.drop.frac and vars.drop.num must be NULL and the other set")
    max.num.steps <- dim(xdata)[2]
    num.subjects <- dim(xdata)[1]
    if (is.null(colnames(xdata))) 
        colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep = "")
    n.vars <- vars <- OOB.rf <- OOB.sd <- rep(NA, max.num.steps)
    oobError <- function(rf) {
        ooo <- rf$confusion[, -dim(rf$confusion)[2]]
        s.ooo <- sum(ooo)
        diag(ooo) <- 0
        sum(ooo)/s.ooo
    }
    if (!is.null(fitted.rf)) {
        if (ncol(fitted.rf$importance) < 2) 
            stop("The fitted rf was not fitted with importance = TRUE")
        n.ntree <- fitted.rf$ntree
        mtry <- fitted.rf$mtry
        n.mtryFactor <- mtry/sqrt(ncol(xdata))
        if ((n.ntree != ntree) | (n.mtryFactor != mtryFactor)) 
            warning("Using as ntree and mtry the parameters obtained from fitted.rf", 
                immediate. = TRUE)
        ntree <- n.ntree
        mtryFactor <- n.mtryFactor
        rm(n.ntree, n.mtryFactor)
        rf <- fitted.rf
    }
    else {
        mtry <- floor(sqrt(ncol(xdata)) * mtryFactor)
        rf <- randomForest(x = xdata, y = Class, ntree = ntree, sampsize=sampsize,
            mtry = mtry, importance = TRUE, keep.forest = keep.forest)
    }
    if (returnFirstForest) 
        FirstForest <- rf
    else FirstForest <- NULL
    m.iterated.ob.error <- m.initial.ob.error <- oobError(rf)
    sd.iterated.ob.error <- sd.initial.ob.error <- sqrt(m.iterated.ob.error * 
        (1 - m.iterated.ob.error) * (1/num.subjects))
    if (verbose) {
        print(paste("Initial OOB error: mean = ", round(m.initial.ob.error, 
            4), "; sd = ", round(sd.initial.ob.error, 4), sep = ""))
    }
    importances <- importance(rf, type = 1, scale = FALSE)
    selected.vars <- order(importances, decreasing = TRUE)
    ordered.importances <- importances[selected.vars]
    initialImportances <- importances
    initialOrderedImportances <- ordered.importances
    j <- 1
    n.vars[j] <- dim(xdata)[2]
    vars[j] <- paste(colnames(xdata), collapse = " + ")
    OOB.rf[j] <- m.iterated.ob.error
    OOB.sd[j] <- sd.iterated.ob.error
    var.simplify <- TRUE
    while (var.simplify) {
        if (verbose) {
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
            importances <- importance(rf, type = 1, scale = FALSE)
            tmp.order <- order(importances, decreasing = TRUE)
            selected.vars <- selected.vars[tmp.order]
            ordered.importances <- importances[tmp.order]
        }
        num.vars <- length(selected.vars)
        if (is.null(vars.drop.num)) 
            vars.drop <- round(num.vars * vars.drop.frac)
        else vars.drop <- vars.drop.num
        if (num.vars >= (vars.drop + 2)) {
            selected.vars <- selected.vars[1:(num.vars - vars.drop)]
            ordered.importances <- ordered.importances[1:(num.vars - 
                vars.drop)]
        }
        else {
            selected.vars <- selected.vars[1:2]
            ordered.importances <- ordered.importances[1:2]
        }
        if ((length(selected.vars) < 2) | (any(selected.vars < 
            1))) {
            var.simplify <- FALSE
            break
        }
        mtry <- floor(sqrt(length(selected.vars)) * mtryFactor)
        if (mtry > length(selected.vars)) 
            mtry <- length(selected.vars)
        if (recompute.var.imp) 
            rf <- randomForest(x = xdata[, selected.vars], y = Class, sampsize=sampsize,
                importance = TRUE, ntree = ntree, mtry = mtry, 
                keep.forest = keep.forest)
        else rf <- randomForest(x = xdata[, selected.vars], y = Class, sampsize=sampsize,
            importance = FALSE, ntree = ntreeIterat, mtry = mtry, 
            keep.forest = keep.forest)
        m.iterated.ob.error <- oobError(rf)
        sd.iterated.ob.error <- sqrt(m.iterated.ob.error * (1 - 
            m.iterated.ob.error) * (1/num.subjects))
        if (verbose) {
            print(paste("..... iteration ", j, "; OOB error: mean = ", 
                round(m.iterated.ob.error, 4), "; sd = ", round(sd.iterated.ob.error, 
                  4), "; num. vars = ", length(selected.vars), 
                sep = ""))
        }
        j <- j + 1
        n.vars[j] <- length(selected.vars)
        vars[j] <- paste(colnames(xdata)[selected.vars], collapse = " + ")
        OOB.rf[j] <- m.iterated.ob.error
        OOB.sd[j] <- sd.iterated.ob.error
        if (!whole.range & ((m.iterated.ob.error > (m.initial.ob.error + 
            c.sd * sd.initial.ob.error)) | (m.iterated.ob.error > 
            (previous.m.error + c.sd * previous.sd.error)))) 
            var.simplify <- FALSE
    }
    if (!whole.range) {
        if (!is.null(colnames(xdata))) 
            selected.vars <- sort(colnames(xdata)[last.vars])
        else selected.vars <- last.vars
        out <- list(selec.history = data.frame(Number.Variables = n.vars, 
            Vars.in.Forest = vars, OOB = OOB.rf, sd.OOB = OOB.sd)[1:j, 
            ], rf.model = last.rf, selected.vars = selected.vars, 
            selected.model = paste(selected.vars, collapse = " + "), 
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
        best.pos <- which(OOB.rf <= min.oob.ci)[which.min(n.vars[which(OOB.rf <= 
            min.oob.ci)])]
        selected.vars <- sort(unlist(strsplit(vars[best.pos], 
            " + ", fixed = TRUE)))
        out <- list(selec.history = data.frame(Number.Variables = n.vars, 
            Vars.in.Forest = vars, OOB = OOB.rf, sd.OOB = OOB.sd), 
            rf.model = NA, selected.vars = selected.vars, selected.model = paste(selected.vars, 
                collapse = " + "), best.model.nvars = n.vars[best.pos], 
            initialImportances = initialImportances, initialOrderedImportances = initialOrderedImportances, 
            ntree = ntree, ntreeIterat = ntreeIterat, mtryFactor = mtryFactor, 
            firstForest = FirstForest)
        class(out) <- "varSelRF"
        return(out)
    }
}
