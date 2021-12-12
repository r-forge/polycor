# last modified 2021-12-11 by J. Fox

hetcor <- function(data, ..., ML=FALSE, std.err=TRUE, 
         use=c("complete.obs", "pairwise.complete.obs"), bins=4,
         pd=TRUE, parallel=FALSE, ncores=detectCores(logical=FALSE),
         thresholds=FALSE){
  UseMethod("hetcor")
  }
