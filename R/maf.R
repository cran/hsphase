.maf <- function(snp)
{
	result <- .Call("MAFC", snp, PACKAGE = "hsphase")
	result
} 