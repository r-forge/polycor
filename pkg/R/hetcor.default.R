# last modified 5 Dec 07 by J. Fox

"hetcor.default" <-
function(data, ..., ML=FALSE, std.err=TRUE, bins=4, pd=TRUE){
  dframe <- data.frame(data, ...)
  names(dframe)[1] <- deparse(substitute(data))
  hetcor(dframe, ML=ML, std.err=std.err, bins=bins, pd=pd)
  }
