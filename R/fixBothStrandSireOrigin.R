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
#' Fix conflicting allele assignments on both haplotype strands
#'
#' Internal helper that scans a haplotype/strand matrix in *pairs of rows*
#' (two rows per individual/segment) and detects loci where **both rows in a
#' pair are non-zero**. Such positions are treated as conflicting assignments
#' and are replaced with the sentinel value `100` on *both* rows of the pair.
#'
#' This function is a thin R wrapper around the native routine
#' \code{fixBothStrand} implemented in C++ and called via \code{.Call()}.
#'
#' @param groups A numeric matrix with an **even number of rows**. Rows are
#' interpreted in pairs: rows \code{1,2} form the first pair, \code{3,4} the
#' second, and so on. Columns represent markers/positions. Values are typically
#' strand/allele codes; `0` is treated as "unknown/empty", and any non-zero
#' value indicates an assigned state.
#'
#' @return A numeric matrix with the same dimensions as \code{groups}. For each
#' row-pair and column:
#' \itemize{
#'   \item if both entries are non-zero, both are set to `100`;
#'   \item otherwise, the original values are preserved.
#' }
#'
#' @details
#' The sentinel value `100` is used downstream as an indicator of an
#' unresolvable conflict between the two strands at that position.
#'
#' @seealso \code{\link{aio}}, \code{\link{ssp}}, \code{\link{bmh}} for the main
#' hsphase workflow where strand/haplotype matrices are produced and refined.
#'
#' @keywords internal
.fixBothStrand <- function(groups)
{
		result <- .Call( "fixBothStrand", groups, PACKAGE = "hsphase" )	
		result	
}

