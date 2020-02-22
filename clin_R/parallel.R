# Functions for parallel processing
# 
# 01-Aug-14 Yen Low
###############################################################################
#Call script by:
#source("<RScriptPath>/R/parallel.R")

#Remember to module load openmpi in cluster (Sherlock)
#module loadopenmpi/1.6.5

#RScriptPath=Sys.getenv("RScriptPath")
source(paste(RScriptPath,"/R/utils.R",sep="")) #installnewpackage
installnewpackage(c("Rmpi","doMPI","doParallel"))
require(Rmpi)   #use Rmpi and doMPI for computing cluster (e.g. Sherlock)
require(doMPI)
#require(doParallel)    #use doParallel for multicore server (e.g dev2)

#getNodeInfo (for )
getNodeInfo <- function() {
  info <- Sys.info()[c("nodename", "machine")]
  print(paste("Node: ", info[1], " | CPU type: ", info[2]))
}