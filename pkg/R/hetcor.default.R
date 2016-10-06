# last modified 2016-10-05 by J. Fox

hetcor.default <-
	function (data, ..., ML = FALSE, std.err = TRUE, 
	          use=c("complete.obs", "pairwise.complete.obs"), bins = 4, pd = TRUE) 
{
	use <- match.arg(use)
	dframe <- data.frame(data, ...)
	if (!missing(...)) names(dframe)[1] <- deparse(substitute(data))
	hetcor(dframe, ML = ML, std.err = std.err, use=use, bins = bins, pd = pd)
}

