pogc <- function(oh, genotypeError)
{
    diag(oh) <- NA
    halfsib <- apply(oh, 1, function(x) names(which(x < genotypeError)))
    halfsib <- halfsib[unlist(lapply(halfsib, length)) > 2]
    if (length(halfsib) > 0)
    {
        repeated <- lapply(halfsib, length)
        pedigree <- data.frame(unlist(halfsib), rep(names(halfsib), repeated))
        colnames(pedigree) <- c("offspring", "parent")
        return(pedigree)
    }
    print("Could not find any groups!")
} 
