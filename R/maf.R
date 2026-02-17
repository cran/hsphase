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

#' Calculate minor allele frequency (MAF)
#'
#' Calculates the minor allele frequency (MAF) for a single SNP coded as:
#' 0 = AA, 1 = AB, 2 = BB, and 9 = missing.
#'
#' @param snp A numeric vector of genotypes for one SNP. Values must be 0, 1, 2,
#'   or 9 (missing).
#'
#' @return A single numeric value: the minor allele frequency (MAF).
#'
#' @examples
#' snp_data <- c(0, 0, 1, 2, 2, 9)
#' .maf(snp_data)
.maf <- function(snp) {
  result <- .Call("MAFC", snp, PACKAGE = "hsphase")
  result
}
