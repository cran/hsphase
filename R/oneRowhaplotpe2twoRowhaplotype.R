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
#' Convert a one-row haplotype to a two-row haplotype
#'
#' Converts a haplotype matrix where each individual is represented by one row
#' and alleles are stored in alternating columns (1st allele, 2nd allele, ...)
#' into a two-row-per-individual representation.
#'
#' Internally, any allele code of `2` is converted to `0` before conversion.
#'
#' @param haplotype A haplotype object:
#' \itemize{
#'   \item a matrix with individuals in rows and allele columns in pairs
#'   (i.e. \code{ncol(haplotype)} must be even).
#' }
#'
#' @return An integer matrix in a two-row-per-individual format with
#' \code{2 * nrow(haplotype)} rows and \code{ncol(haplotype) / 2} columns.
#' Row names are interleaved using the original individual names.
#'
#' @export
.o2tH <- .one2tworowhap <- function(haplotype)
{
	if(is.integer(haplotype) && !is.matrix(haplotype))
		haplotype <- t(as.matrix(haplotype))
	haplotype[haplotype==2] <- 0
	result <- matrix(0,nrow = 2 * nrow(haplotype),ncol = ncol(haplotype) / 2)
	result[seq(1,nrow(result),2),] <- haplotype[,seq(1,ncol(haplotype),2),drop=FALSE]
	result[seq(2,nrow(result),2),] <- haplotype[,seq(2,ncol(haplotype),2),drop=FALSE]
	rownames(result) <- gdata::interleave(rownames(haplotype),paste(rownames(haplotype),1:nrow(haplotype),sep="_"))
  result
}
