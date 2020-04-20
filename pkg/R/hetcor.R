# last modified 2020-04-20 by J. Fox

hetcor <- function(data, ..., ML=FALSE, std.err=TRUE, 
         use=c("complete.obs", "pairwise.complete.obs"), bins=4,
         pd=TRUE, parallel=FALSE, ncores=detectCores(logical=FALSE)){
  UseMethod("hetcor")
  }
