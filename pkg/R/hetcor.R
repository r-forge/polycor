# last modified 2020-04-19 by J. Fox

hetcor <- function(data, ..., ML=FALSE, std.err=TRUE, 
         use=c("complete.obs", "pairwise.complete.obs"), bins=4, pd=TRUE, ncores=1){
  UseMethod("hetcor")
  }
