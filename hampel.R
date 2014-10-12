# Hampel Identifier for Outlier Detection
# 
# 02-Dec-13
# http://global.oup.com/us/companion.websites/fdscontent/uscompanion/us/static/companion.websites/9780195089653/TextFiles/hampelproc.txt
# Usage: 
#RScriptPath=Sys.getenv("RScriptPath")
#source(paste(RScriptPath,"/scripts/R/hampel.R",sep=""))
###############################################################################


hampel <- function(x, t = 3, RemoveNAs = FALSE,two.sided=TRUE){
	#
	#  This procedure returns an index of x values declared
	#  outliers according to the Hampel detection rule, if any
	#
	mu <- median(x, na.rm = RemoveNAs)
	sig <- mad(x, na.rm = RemoveNAs)
	indx <- which( abs(x - mu) > t*sig)
	if(!two.sided)  indx <- which( x - mu > t*sig)
	#
	indx
}

