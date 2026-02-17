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

#' Fix strand label rotation across consecutive block-structure columns
#'
#' Internal helper to enforce a consistent strand-label orientation across
#' adjacent columns of a block-structure matrix.
#'
#' The input typically encodes sire strand-of-origin labels per individual
#' (rows) and marker/SNP (columns), where `0` indicates unknown and non-zero
#' values indicate an assigned strand/state. The native algorithm compares each
#' column to the previous one and, when a "contrast" (swap of strand labels)
#' increases agreement, it relabels the next column to reduce apparent
#' strand-rotation between columns.
#'
#' This function is a thin R wrapper around the native routine
#' \code{fixRotation} implemented in C++ and called via \code{.Call()}.
#'
#' @param blockStructure A numeric matrix (typically individuals in rows, SNPs
#' in columns) representing a block/strand structure. Values are expected to be
#' small integers (commonly including `0`, `1`, `2` and possibly other internal
#' codes).
#'
#' @return A numeric matrix with the same dimensions as \code{blockStructure},
#' where some entries in column \code{i+1} may be relabeled to improve
#' consistency with column \code{i}. The transformation is applied iteratively
#' from left to right across columns.
#'
#' @details
#' At each step, the C++ code computes an agreement score between column
#' \code{i} and column \code{i+1} using only positions where both columns are
#' non-zero. It also computes the score after applying a contrast mapping to
#' column \code{i+1} (conceptually swapping strand labels `1` and `2`, leaving
#' `0` unchanged). If the contrasted version agrees more with column \code{i},
#' the function relabels column \code{i+1}.
#'
#' The relabeling performed by the native code is:
#' \itemize{
#'   \item `1 -> 3`
#'   \item `2 -> 1`
#'   \item `3 -> 2`
#' }
#' leaving other values unchanged. (These codes are part of hsphase's internal
#' block/strand encoding.)
#'
#' @seealso \code{\link{bmh}}, \code{\link{ssp}}, \code{\link{aio}} for creation
#' and downstream usage of block structures.
#'
#' @keywords internal

.fixRotation <- function(blockStructure)
{
		result <- .Call( "fixRotation", blockStructure, PACKAGE = "hsphase" )
		result
}
