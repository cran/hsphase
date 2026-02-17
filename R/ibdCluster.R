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
#' Cluster individuals by local IBD-like similarity in sliding windows
#'
#' Internal helper that applies a sliding-window clustering procedure across a
#' genotype matrix. For each window of markers, individuals are grouped using
#' an opposing-homozygote style criterion (controlled by \code{maxOH}) via a
#' native C++ routine.
#'
#' This function is a thin R wrapper around the native routine
#' \code{ibdCluster} implemented in C++ and called via \code{.Call()}.
#'
#' @param genotype A numeric matrix of genotypes with individuals in rows and
#' markers/SNPs in columns. This function treats genotype code `1` (heterozygote)
#' as missing by converting it to `9` prior to calling the native routine.
#'
#' @param cpus Integer scalar. Requested number of CPU threads. The underlying
#' C++ implementation uses OpenMP; the actual number of threads may depend on
#' how R and your compiler toolchain were built and configured.
#'
#' @param windowsSize Integer scalar giving the window size (in markers). Must
#' be between 1 and \code{ncol(genotype)}. Internally, the native routine uses
#' a 0-based offset window length of \code{windowsSize - 1}.
#'
#' @param maxOH Integer scalar. Maximum allowed opposing-homozygote
#' count (or a similar mismatch threshold) used by the native grouping
#' algorithm within each window.
#'
#' @return A list (R "list" object) of length \code{ncol(genotype) - (windowsSize - 1)}.
#' Each element corresponds to one sliding window position and contains an
#' integer vector of group assignments produced by the native implementation.
#' (Exact encoding of groups is defined by the underlying \code{mhMat} grouping
#' routine.)
#'
#' @details
#' The function scans windows of consecutive markers. For each window starting
#' at marker \code{i}, it clusters individuals based on the submatrix
#' \code{genotype[, i:(i+windowsSize-1)]} after recoding heterozygotes (`1`) to
#' missing (`9`). Computation is parallelized in C++ using OpenMP across window
#' start positions.
#'
#' @keywords internal
.ibdCluster <- function(genotype, cpus, windowsSize, maxOH)
{
    genotype[genotype == 1] <- 9
    if (windowsSize > ncol(genotype))
    {
        stop("The windows size must be equal or less than number of columns")
    }
    if (windowsSize < 1)
    {
        stop("The windows size must be greater than 1")
    }
	if(maxOH < 0)
	{
		stop("maxOH must be equal or greater than 0");
	}
    windowsSize <- windowsSize - 1
    result <- .Call("ibdCluster", genotype, cpus, windowsSize, maxOH, PACKAGE = "hsphase")
    result
} 
