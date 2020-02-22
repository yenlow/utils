# Functions:
# 	1. installnewpackage: installnewpackage ONLY if they have not been installed
# 	2. lsos: lists objects in memory, sorted by size
#   3. gci: gc() until memory is clean
# 
# 21-Jan-14 Yen Low
#
# To use:
# RScriptPath=Sys.getenv("RScriptPath") #Rscriptpath is an environment variable, path="/mnt/hgfs/Dropbox"
# source(paste(RScriptPath,"/scripts/R/utils.R",sep=""))
###############################################################################

############### installnewpackage ################# 
#Check if required packages have been installed.
# If not, install the missing packages

#EXAMPLES:
##required packages are "caret","ROCR","boot"
#reqpackages=c("caret","ROCR","boot")
#installnewpackage(reqpackages)

installnewpackage<-function(reqpackages){
	newpackages=reqpackages[!(reqpackages %in% installed.packages()[,"Package"])]
	if(length(newpackages)>0) install.packages(newpackages)
}


#################### lsos #######################
# Function for listing objects in memory, sorted by size
# From http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
		decreasing=FALSE, head=FALSE, n=5) {
	napply <- function(names, fn) sapply(names, function(x)
					fn(get(x, pos = pos)))
	names <- ls(pos = pos, pattern = pattern)
	obj.class <- napply(names, function(x) as.character(class(x))[1])
	obj.mode <- napply(names, mode)
	obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
	obj.size <- napply(names, object.size)
	obj.dim <- t(napply(names, function(x)
						as.numeric(dim(x))[1:2]))
	vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
	obj.dim[vec, 1] <- napply(names, length)[vec]
	out <- data.frame(obj.type, obj.size, obj.dim)
	names(out) <- c("Type", "Size", "Rows", "Columns")
	if (!missing(order.by))
		out <- out[order(out[[order.by]], decreasing=decreasing), ]
	if (head)
		out <- head(out, n)
	out
}
# shorthand
lsos <- function(..., n=10) {
	.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}


################ collectGarbage ############################
## The following function collects garbage until the memory is clean.
## Usage: 1. immediately call this function after you call a function or
##        2. rm()
gci <- function(){
	while (gc()[2,4] != gc()[2,4]){}
}
