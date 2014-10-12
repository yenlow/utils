#http://www.ddiez.com/teac/surv/time-dep-coxph.R
# To use:
# source("/mnt/hgfs/Dropbox/scripts/R/timeDepCoxph.R")
#######################################################3
d.f=ds_matched
col.time="days"
col.delta="dead"
col.cov=c("AGE_AT_DX","daysDM","hormone","white","LUMPECTOMY","MASTECTOMY","EVIDENCE_OF_CHEMO","EVIDENCE_OF_HORMONE","EVIDENCE_OF_RT")
td.cov="Metformin365"
temp=time.dep.coxph(d.f, col.time, col.delta, col.cov, td.cov, transform=log,
		method='breslow', output.model=TRUE, output.data.frame=TRUE, verbose=TRUE)

time.dep.coxph <- function(d.f, col.time, col.delta, col.cov, td.cov, transform=log,
		method='breslow', output.model=TRUE, output.data.frame=FALSE, verbose=TRUE){
	if( length(col.cov) > 1){
		cov.names <- colnames(d.f[,col.cov])
	} else {
		if(is.numeric(col.cov)){
			cov.names <- colnames(d.f)[col.cov]
		} else {
			cov.names <- col.cov
		}
	}
	
	surv <- d.f[ ,col.time]
	obs <- d.f[ ,col.delta]
	n.old.rows <- dim(d.f)[1]
	max.event <- as.integer( max(surv)+1 )
	n.entries <- sum(as.integer(surv)) + n.old.rows
	t.cov <- matrix(0, nrow=n.entries, ncol=length(col.cov)+4)
	colnames(t.cov) <- c("start", "stop", "observ", cov.names, "time.dep.cov")
	row <- 0
	hold <- 0
	for(i in 1:n.old.rows){
		for(j in 1:max.event){
			if(surv[i] > (j-1)){
				row <- row+1
				t.cov[row,"start"] <- j-1
				t.cov[row,"stop"] <- j
				if(surv[i] <= t.cov[row,"stop"])
					t.cov[row,"observ"] <- obs[i]
				t.cov[row,cov.names] <- as.numeric(d.f[i,col.cov])
				t.cov[row,"time.dep.cov"] <- d.f[i,td.cov]*transform(j)
			}
		}
		if(verbose){
			temp <- as.integer(100*row/n.entries)
			if(temp > hold+4 & temp < 100){
				cat("data frame is ",temp,"% complete.\n", sep="")
				hold <- temp
			}
		}
	}
	t.cov <- as.data.frame(t.cov)
	if(verbose) cat("Data frame completed.\n")
	if(output.model){
		cat("The covariates in the output will be of the same order as specified in 'col.cov'.\n")
		if(verbose) cat("\nFitting model.\n\n")
		. <- as.matrix(t.cov[,c(cov.names,"time.dep.cov")])
		coxph.fit <- coxph(Surv(t.cov$start, t.cov$stop, t.cov$observ) ~ ., data=t.cov, method=method)
		if(!output.data.frame){
			coxph.fit
		} else {
			list(coxph.fit=coxph.fit, time.dep.df=t.cov)
		}
	} else if(output.data.frame) t.cov
}
