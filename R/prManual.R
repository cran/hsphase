# Copyright (C) 2014 Mohammad H. Ferdosi
#
# HSPhase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# HSPhase program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http:#www.gnu.org/licenses/>.
#' Manual recursive clustering using an OH matrix and a fixed threshold
#'
#' Performs a recursive hierarchical clustering on an opposing-homozygotes (OH)
#' matrix using Ward clustering and splits clusters until the maximum within-
#' cluster OH value is below a user-supplied threshold (\code{maxsnpnooh}).
#'
#' The function returns a two-column data frame with individual IDs and a group
#' code. Group codes are generated randomly (via \code{rnorm()}) and therefore
#' are not stable across runs.
#'
#' @param oh A numeric matrix representing opposing-homozygotes (OH) counts
#'   between individuals. Row and column names should be individual IDs. The
#'   matrix is expected to be square and symmetric.
#' @param maxsnpnooh Numeric scalar giving the maximum allowed within-group OH
#'   value. Groups with a larger maximum OH value (and size > 2) are split
#'   recursively.
#'
#' @details
#' The recursion proceeds as follows:
#' \enumerate{
#'   \item Compute a pairwise distance object from \code{oh} using \code{.fastdist}
#'   and convert it to a \code{dist} object.
#'   \item Apply hierarchical clustering using \code{\link[stats]{hclust}} with
#'   \code{method = "ward.D"}.
#'   \item Cut the dendrogram into two groups using \code{\link[stats]{cutree}}.
#'   \item For each group, compute the maximum within-group OH value; if it
#'   exceeds \code{maxsnpnooh} and the group has more than two individuals,
#'   recurse into that subgroup. Otherwise, write group assignments and stop.
#' }
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{id}{Individual ID (character).}
#'   \item{group}{An integer-like group code (generated randomly; not reproducible).}
#' }
#'
#' @section Side effects:
#' This function writes to and reads from a file named \code{"temp.txt"} in the
#' current working directory, and then deletes it. It also prints
#' \code{maxsnpnooh} to the console.
#'
#' @seealso \code{\link[stats]{hclust}}, \code{\link[stats]{cutree}},
#'   \code{\link[stats]{as.dist}}
#'
#' @keywords internal
.prManual <- function(oh, maxsnpnooh)
{
	cat("id group \n", file = "temp.txt")
	rhsr_rc <- function(oh, maxsnpnooh = maxsnpnooh)
	{
		print("----")
		## d <- dist(oh, method = 'manhattan')
		d <- as.dist(.fastdist(oh))
		if (length(d) > 2)
		{
			fit <- hclust(d, method = "ward.D")
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
				
				utils::write.table(data.frame(names(a), round(abs(rnorm(1) * 10^5))), "temp.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
			}
			if (maxSubohB > maxsnpnooh && length(b) > 2)
			{
				rhsr_rc(oh[b, b], maxsnpnooh)
			} else
			{
				utils::write.table(data.frame(names(b), round(abs(rnorm(1) * 10^6))), "temp.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
			}
		} else
		{
			if (!is.integer(oh)) 
				utils::write.table(data.frame(rownames(oh), round(abs(rnorm(1) * 10^6))), "temp.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
		}
	}
	result <- rhsr_rc(oh, maxsnpnooh)
	result <- utils::read.table("temp.txt", header = TRUE)
	file.remove("temp.txt")
	print(maxsnpnooh)
	
	result
} 
