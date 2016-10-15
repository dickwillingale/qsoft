# qR start up - the environment variable QSOFT must be set
# This file to be put in your home directory as file .Rprofile
# If will be executed when R or Rscript starts
# Dick Willingale 2016-Sep-28
.First <- function() {
# Pick up top directory
  qsoft<- Sys.getenv("QSOFT")
  qversion<- Sys.getenv("QVERSION")
# Load the shareable libraries
  dyn.load(paste(qsoft,"/R_libraries/qfitsR.so",sep=""))
  dyn.load(paste(qsoft,"/R_libraries/imagesR.so",sep=""))
  dyn.load(paste(qsoft,"/R_libraries/astroR.so",sep=""))
  dyn.load(paste(qsoft,"/R_libraries/xsrtR.so",sep=""))
  dyn.load(paste(qsoft,"/R_libraries/xscatR.so",sep=""))
# Source the R definitions
  source(paste(qsoft,"/R_libraries/qfits.R",sep=""))
  source(paste(qsoft,"/R_libraries/images.R",sep=""))
  source(paste(qsoft,"/R_libraries/astro.R",sep=""))
  source(paste(qsoft,"/R_libraries/xsrt.R",sep=""))
  source(paste(qsoft,"/R_libraries/xscat.R",sep=""))
#
# Quick get out of jail
  qq<-function() {
    q(save="no")
    invisible()
  }
  cat("Q version",qversion,": qfits images astro xsrt xscat loaded\n")
}
