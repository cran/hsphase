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
#' Reconstruct half-sib groups by recursive clustering with a recombination stop rule
#'
#' Internal helper used by \code{\link{rpoh}} (reconstruction pedigree of half-sib
#' families). This function recursively splits individuals into two clusters using
#' hierarchical clustering on a distance derived from the provided opposing
#' homozygote (OH) matrix, and then decides whether each cluster should be split
#' further by checking the maximum number of recombination events inferred within
#' that cluster.
#'
#' The recursive splitting stops for a cluster when the maximum recombination
#' count in that cluster is \code{<= maxRec}. Final group assignments are written
#' to a temporary file and then read back as a two-column data frame.
#'
#' @param genotypeMatrix Numeric genotype matrix (individuals in rows, SNPs in
#' columns) coded as `0`, `1`, `2` (and typically `9` for missing), as used by
#' hsphase. This matrix is subset recursively when splitting clusters.
#'
#' @param oh A square opposing-homozygote matrix for the same individuals as
#' \code{genotypeMatrix} (rownames/colnames are individual IDs). Typically
#' produced by \code{\link{ohg}}. This matrix is subset recursively along with
#' \code{genotypeMatrix}.
#'
#' @param forwardVectorSize Integer. Passed to \code{\link{bmh}} when computing
#' recombination blocks inside each candidate cluster.
#'
#' @param excludeFP Logical. Passed to \code{\link{bmh}}.
#'
#' @param nsap Integer. Passed to \code{\link{bmh}}.
#'
#' @param maxRec Integer. Maximum allowed recombination count (within a cluster)
#' before the cluster is recursively split again.
#'
#' @return A \code{data.frame} with two columns:
#' \itemize{
#'   \item \code{id}: individual IDs
#'   \item \code{group}: an integer-like group label assigned by the recursive
#'   procedure
#' }
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Converts \code{oh} to a distance object via \code{as.dist(.fastdist(oh))}
#'   and performs hierarchical clustering (\code{hclust}, Ward method).
#'   \item Splits into \code{k = 2} clusters via \code{cutree}.
#'   \item For each cluster with at least 4 individuals, computes recombination
#'   counts as \code{recombinations(bmh(subGenotype, ...))} and uses the maximum
#'   recombination count as a stop/split criterion.
#'   \item If \code{max(recombinations) > maxRec}, the cluster is split again
#'   recursively; otherwise, individuals in that cluster are assigned a new group
#'   label and written to a temporary file.
#' }
#'
#' @section Implementation notes:
#' \itemize{
#'   \item This function uses a fixed temporary filename \code{"temp.txt"} in the
#'   current working directory and deletes it at the end. This is not safe under
#'   parallel execution or if the working directory is not writable.
#'   \item Group labels are generated using \code{rnorm()}, so results are not
#'   deterministic unless a seed is set and the recursion order remains identical.
#' }
#'
#' @seealso \code{\link{rpoh}}, \code{\link{ohg}}, \code{\link{bmh}},
#' \code{\link{recombinations}}, \code{\link{.fastdist}}
#'
#' @keywords internal
.rpohhsphase <- function(genotypeMatrix, oh, forwardVectorSize = 30, excludeFP = TRUE, nsap = 3, maxRec = 15)
{
	
	cat("id group \n",file="temp.txt")
	rhsr_rc <- function(genotypeMatrix ,oh)
	{
		d <- as.dist(.fastdist(oh)) 
		if(length(d)>2)
		{
			fit <- hclust(d, method = "ward")
			groups <- cutree(fit, k = 2)
			a <- which(groups==1)
			b <- which(groups==2)
			# print(length(a))
			if(length(a)>3)
			{
				rec_1 <- recombinations(bmh(genotypeMatrix[a,],excludeFP = excludeFP,nsap = nsap,forwardVectorSize = forwardVectorSize))
				rec_1 <- max(rec_1)
				#print(paste("rec_1",rec_1))
			}
			else
			{
				rec_1 <- 9
			}
			print(length(b))
			if(length(b)>3)
			{
				rec_2 <- recombinations(bmh(genotypeMatrix[b,],excludeFP = excludeFP,nsap = nsap,forwardVectorSize = forwardVectorSize))
				rec_2 <- max(rec_2)
				#print(paste("rec_1",rec_2))
			}
			else
			{
				rec_2 <- 9
			}
			if(rec_1>maxRec)
			{
				rhsr_rc(genotypeMatrix[a,],oh[a,a])
			}
			else
			{		
				utils::write.table(data.frame(names(a),round(abs(rnorm(1)*10^5))),"temp.txt",append = TRUE,col.names = FALSE,row.names = FALSE)
			}
			if(rec_2>maxRec)
			{
				rhsr_rc(genotypeMatrix[b,],oh[b,b])
			}
			else
			{				
				utils::write.table(data.frame(names(b),round(abs(rnorm(1)*10^6))),"temp.txt",append = TRUE,col.names = FALSE,row.names = FALSE)
			}
		}
		else
		{
			if(!is.integer(oh))
			{
				utils::write.table(data.frame(rownames(oh),round(abs(rnorm(1)*10^6))),"temp.txt",append = TRUE,col.names = FALSE,row.names = FALSE)			
			}
			print(paste(class(oh),nrow(oh)))
		
		}
	}
	
	
	result <- rhsr_rc(genotypeMatrix, oh)
	
	result <- utils::read.table("temp.txt",header = TRUE)
	file.remove("temp.txt")
	result
}
