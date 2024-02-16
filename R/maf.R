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

#' Calculate Minor Allele Frequency (MAF)
#'
#' This function calculates the minor allele frequency (MAF) for a given single nucleotide polymorphism (SNP) data. 
#' The SNP data should be coded numerically: 0 for homozygous for the first allele (AA), 
#' 1 for heterozygous (AB), and 2 for homozygous for the second allele (BB). Missing data should be coded as 9.
#'
#' @param snp A numeric vector representing the genotype of individuals for a single SNP. 
#' The genotype should be coded as 0 for AA, 1 for AB, and 2 for BB. Use 9 to represent missing data.
#' 
#' @return A numeric value representing the minor allele frequency (MAF) for the SNP data provided.
#'
#' @examples
#' snp_data <- c(0, 0, 1, 2, 2, 9)
#' maf_value <- .maf(snp_data)
#' print(maf_value)
#'
#' @export
#' @useDynLib hsphase, .registration=TRUE
#' @importFrom Rcpp .Call
.maf <- function(snp) {
  result <- .Call("MAFC", snp, PACKAGE = "hsphase")
  result
}
