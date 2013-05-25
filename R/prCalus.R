.prCalus <- function(oh, genotype)
{
	maf_geno <- function(x)
	{
		z <- length(which(x == 0))
		o <- length(which(x == 1))
		maf <- (z * 2 + o)/(sum(!is.na(x)) * 2)
		if (!is.na(maf)) 
			if (maf > 0.5) 
				maf <- 1 - maf
		maf
	}
	
	p <- apply(genotype, 2, .maf)
	p <- p[!is.na(p)]
	q <- 1 - p
	
	maxsnpnooh <- (sum(p^2 * q^2) + sum(2 * (p^2 * q^2)))/2
	maxsnpnooh <- maxsnpnooh - (.1 * maxsnpnooh)
	
	cat("id group \n", file = "temp.txt")
	rhsr_rc <- function(oh, maxsnpnooh = maxsnpnooh)
	{
		print("----")
		## d <- dist(oh, method = 'manhattan')
		d <- as.dist(.fastdist(oh))
		if (length(d) > 2)
		{
			fit <- hclust(d, method = "ward")
			groups <- cutree(fit, k = 2)
			a <- which(groups == 1)
			b <- which(groups == 2)
			
			
			if (length(a) > 2)
			{
				subohA <- oh[a, a]
				maxSubohA <- max(subohA[lower.tri(subohA)])
			} else
			{
				maxSubohA <- 0
			}
			
			if (length(b) > 2)
			{
				subohB <- oh[b, b]
				maxSubohB <- max(subohB[lower.tri(subohB)])
			} else
			{
				maxSubohB <- 0
			}
			
			if (maxSubohA > maxsnpnooh && length(a) > 2)
			{
				
				rhsr_rc(oh[a, a], maxsnpnooh)
			} else
			{
				
				write.table(data.frame(names(a), round(abs(rnorm(1) * 10^5))), "temp.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
			}
			if (maxSubohB > maxsnpnooh && length(b) > 2)
			{
				rhsr_rc(oh[b, b], maxsnpnooh)
			} else
			{
				write.table(data.frame(names(b), round(abs(rnorm(1) * 10^6))), "temp.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
			}
		} else
		{
			if (!is.integer(oh)) 
				write.table(data.frame(rownames(oh), round(abs(rnorm(1) * 10^6))), "temp.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
		}
	}
	result <- rhsr_rc(oh, maxsnpnooh)
	result <- read.table("temp.txt", header = TRUE)
	file.remove("temp.txt")
	result
} 
