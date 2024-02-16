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
#' Calculate Genotypic Distances
#'
#' Calculates a symmetric matrix of distances between genotypes, based on a given genotype matrix.
#' Each row in the `GenotypeMatrix` represents a genotype, and each column represents a marker.
#' The genotype is coded as 0 for AA, 1 for AB, and 2 for BB. Use 9 to represent missing data.
#'
#' @param GenotypeMatrix A matrix where each row represents a genotype and each column 
#' represents a marker. Genotypes should be coded as 0 for AA, 1 for AB, and 2 for BB, 
#' with 9 representing missing data.
#'
#' @return Returns a symmetric matrix of distances between the genotypes specified in the 
#' `GenotypeMatrix`. Row and column names of the returned matrix correspond to the row names 
#' of the `GenotypeMatrix`.
#'
#' @examples
#' # Simulate genotype data for 40 individuals across 1000 SNPs
#' # genotypes <- simulateHalfsib(numInd = 40, numSNP = 1000, recbound = 0:6, type = "genotype")
#' # Calculate the distance matrix
#' # dist_matrix <- fastdist(genotypes)
#' # Display the distance matrix
#' # print(dist_matrix)
#'
#' @export
.fastdist <- function(GenotypeMatrix)
{
	n <- nrow(GenotypeMatrix) * nrow(GenotypeMatrix)
	fMat <- matrix(as.integer(rep(0, n)), nrow = nrow(GenotypeMatrix))
	expandMat <- as.numeric(t(GenotypeMatrix))
	result <- .C("fastDist", expandMat = as.integer(expandMat), nrow = as.integer(nrow(GenotypeMatrix)), 
			ncol = as.integer(ncol(GenotypeMatrix)), result = fMat)$result
	result[upper.tri(result)] <- t(result)[upper.tri(result)]
	colnames(result)=rownames(GenotypeMatrix)
	rownames(result)=rownames(GenotypeMatrix)
	result
}
