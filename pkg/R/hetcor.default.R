# last modified 2020-04-19 by J. Fox

hetcor.default <- function (data, ..., ML = FALSE, std.err = TRUE, 
	          use=c("complete.obs", "pairwise.complete.obs"), 
	          bins = 4, pd = TRUE, parallel=FALSE, ncores=detectCores(logical=FALSE)) 
{
	use <- match.arg(use)
	dframe <- data.frame(data, ...)
	if (!missing(...)) names(dframe)[1] <- deparse(substitute(data))
	hetcor(dframe, ML = ML, std.err = std.err, use=use, bins = bins, pd = pd, ncores=ncores)
}

