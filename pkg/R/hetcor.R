# last modified 2016-10-05 by J. Fox

"hetcor" <-
function(data, ..., ML=FALSE, std.err=TRUE, 
         use=c("complete.obs", "pairwise.complete.obs"), bins=4, pd=TRUE){
  UseMethod("hetcor")
  }
