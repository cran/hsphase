# Copyright (C) 2013 Mohammad H. Ferdosi
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
#' Fill gaps in a paternal strand vector (native routine wrapper)


#'
#' Internal wrapper around the native C routine \code{fillGap}. It fills
#' zero-valued gaps in a paternal-strand vector by propagating values from
#' neighbouring non-zero positions (see native implementation for exact rules).
#'
#' @param paternalStrandBMH An integer (or numeric) vector representing the
#'   paternal strand state (e.g., output from block-building routines). The
#'   vector is coerced to integer before calling native code.
#'
#' @return An integer vector of the same length as \code{paternalStrandBMH}
#'   containing the gap-filled strand values.
#'
#' @details
#' This function calls compiled code via \code{\link[base]{.C}}:
#' \code{.C("fillGap", ...)}.
#'
#' @seealso \code{\link[base]{.C}}
#'
#' @author mhf
#'
#' @keywords internal
.fillGap <- function(paternalStrandBMH)
{
	if(!is.vector(paternalStrandBMH))
		stop("paternalStrandBMH should be a VECTOR")
	fvec <- numeric(length(paternalStrandBMH))
	result <- .C("fillGap",as.integer(paternalStrandBMH), length(paternalStrandBMH), result = as.integer(fvec))$result
	result
}


#' Build haplotype blocks from a BMH result matrix (native routine wrapper)
#'
#' Internal wrapper around the native C routine \code{hblock}. It transforms a
#' BMH result matrix into a block representation, with an optional maximum block
#' size constraint.
#'
#' @param bmhResult A numeric/integer matrix containing BMH results (block
#'   matching/haplotype-block intermediate output). Must be a matrix.
#' @param MaxBlock Integer scalar. Maximum block size (default: 400).
#'
#' @return A matrix (same general shape as \code{bmhResult}) containing inferred
#'   block structure. Row and column names are propagated from \code{bmhResult}
#'   where available.
#'
#' @details
#' This function transposes and flattens \code{bmhResult} before passing it to
#' compiled code via \code{\link[base]{.C}}:
#' \code{.C("hblock", ...)}.
#'
#'
.hblock <- function(bmhResult, MaxBlock = 400)
{
    if (!is.matrix(bmhResult)) 
        stop("GenotypeMatrix should be a MATRIX")
    expandMat <- as.double(t(bmhResult))
    n <- ncol(bmhResult) * nrow(bmhResult)
    fMat <- matrix(as.integer(rep(0, n)), nrow = ncol(bmhResult))
    result <- .C("hblock", expandMat = as.integer(expandMat), nrow = as.integer(ncol(bmhResult)), ncol = as.integer(nrow(bmhResult)), 
        result = fMat, MB = as.integer(MaxBlock))$result
    colnames(result) <- rownames(bmhResult)

    if (!is.null(colnames(bmhResult))) 
        rownames(result) <- colnames(bmhResult)
    t(result)
}



hbp <- function(PhasedGenotypeMatrix, PhasedSireGenotype, strand = "auto")
{
	if (!is.matrix(PhasedGenotypeMatrix)) 
		stop("PhasedGenotypeMatrix should be a MATRIX")
	if (length(PhasedGenotypeMatrix[PhasedGenotypeMatrix != 0 & PhasedGenotypeMatrix != 1 & 
							PhasedGenotypeMatrix != 9]) > 0) 
		stop("PhasedGenotypeMatrix must contain only 0 and 1 or 9 for missing SNP")
	if (!is.matrix(PhasedSireGenotype)) 
		stop("PhasedSireGenotype should be a MATRIX")
	if (length(PhasedSireGenotype[PhasedSireGenotype != 0 & PhasedSireGenotype != 1 & PhasedSireGenotype != 
							9]) > 0) 
		stop("PhasedSireGenotype must contain only 0 and 1 or 9 for missing SNP")
	if (ncol(PhasedGenotypeMatrix) != ncol(PhasedSireGenotype)) 
		stop("Number of markers in sire and half-sib family must be the same")
	if (nrow(PhasedSireGenotype) != 2) 
		stop("PhasedSireGenotype must have 2 rows")
	
	METHODS <- c("auto", 1,2)
	method <- pmatch(strand, METHODS)
	if (is.na(method)) 
		stop("invalid strand")
	if (method == -1) 
		stop("ambiguous pedigree reconstruction method")
	if (method == 1) 
	{
		str = 0
	}
	else if(method == 2)
	{
		str = 1		
	}
	else if(method == 3)
	{
		str = 2
	}
	expandMat <- as.numeric(t(PhasedGenotypeMatrix))
	n <- ncol(PhasedGenotypeMatrix) * nrow(PhasedGenotypeMatrix)
	fMat <- matrix(as.integer(rep(0, n/2)), nrow = ncol(PhasedGenotypeMatrix))
	result <- .C("hbphased", expandMat = as.integer(expandMat), nrow = as.integer(nrow(PhasedGenotypeMatrix)), 
			ncol = as.integer(ncol(PhasedGenotypeMatrix)), result = fMat, siregenotype = as.integer(t(PhasedSireGenotype)), strand = as.integer(str))$result
	
	if (!is.null(rownames(PhasedGenotypeMatrix))) 
		colnames(result) <- rownames(PhasedGenotypeMatrix)[seq(from = 1, by = 2, to = nrow(PhasedGenotypeMatrix))]
	t(result)
}
phf <- function(GenotypeMatrix, blockMatrix, sirePhasedMatrix)
{
	if (is.null(GenotypeMatrix) | is.null(blockMatrix) | is.null(sirePhasedMatrix)) 
		stop("Invalid input!")
	if ((!is.matrix(GenotypeMatrix)) | (!is.matrix(blockMatrix)) | (!is.matrix(sirePhasedMatrix))) 
		stop("All inputs should be a MATRIX")
	if (length(GenotypeMatrix[GenotypeMatrix != 0 & GenotypeMatrix != 2 & GenotypeMatrix != 
							1 & GenotypeMatrix != 9]) > 0) 
		stop("GenotypeMatrix must contain only 0,1 and 2 or 9 for missing SNPs")
	if (length(blockMatrix[blockMatrix != 0 & blockMatrix != 2 & blockMatrix != 1]) > 0) 
		stop("blockMatrix must contain only 0,1 and 2")
	if (length(sirePhasedMatrix[sirePhasedMatrix != 0 & sirePhasedMatrix != 1 & sirePhasedMatrix != 
							9]) > 0) 
		stop("SireMatrix must contain only 0,1 and 9")
	expandMat <- as.numeric(GenotypeMatrix)
	n <- nrow(GenotypeMatrix) * ncol(GenotypeMatrix)
	fMat <- matrix(as.integer(rep(0, n)), ncol = ncol(GenotypeMatrix))
	result <- .C("phase", genotype = as.integer(GenotypeMatrix), nrow = as.integer(nrow(GenotypeMatrix)), 
			ncol = as.integer(ncol(GenotypeMatrix)), block = as.integer(blockMatrix), sirePhasedMatrix = as.integer(t(sirePhasedMatrix)), 
			result = fMat)$result
   colnames(result) <- colnames(GenotypeMatrix)
   result
}
#' Phase half-sib paternal haplotype using blocks and sire haplotypes (no offspring genotype needed)
#'
#' Internal helper that constructs a half-sib paternal haplotype matrix using:
#' \itemize{
#'   \item a block/strand-of-origin matrix (typically produced by \code{\link{bmh}})
#'   \item a 2-row phased sire haplotype matrix (typically produced by \code{\link{ssp}})
#' }
#'
#' For each marker (column) and individual (row), if the block code is:
#' \itemize{
#'   \item `1`: assign sire haplotype row 1 allele at that marker
#'   \item `2`: assign sire haplotype row 2 allele at that marker
#'   \item `0`: leave as missing (`9`)
#' }
#'
#' This function calls a native C routine (\code{phaseNogenotype}) via
#' \code{.C()}.
#'
#' @param blockMatrix An integer/numeric matrix of block assignments with
#' individuals in rows and markers in columns. Must contain only `0`, `1`, and `2`,
#' where `0` indicates unknown origin.
#'
#' @param sirePhasedMatrix An integer/numeric matrix with **two rows** (the sire
#' haplotypes) and the same number of columns as \code{blockMatrix}. Must contain
#' only `0`, `1`, and `9` (where `9` indicates missing).
#'
#' @return An integer matrix with the same dimensions as \code{blockMatrix},
#' containing the inferred paternal haplotype allele for each individual and
#' marker. Values are `0`/`1` for alleles and `9` for missing/unknown (e.g. where
#' \code{blockMatrix} is `0`).
#'
#' @details
#' The underlying C implementation initializes the entire result matrix to `9`
#' and then fills entries according to \code{blockMatrix}:
#' \itemize{
#'   \item if \code{blockMatrix[j,i] == 1}, then \code{result[j,i] = sirePhasedMatrix[1,i]}
#'   \item if \code{blockMatrix[j,i] == 2}, then \code{result[j,i] = sirePhasedMatrix[2,i]}
#' }
#'
#' @keywords internal
.phfnoGenotype <- function(blockMatrix, sirePhasedMatrix)
{
		
	if (length(blockMatrix[blockMatrix != 0 & blockMatrix != 2 & blockMatrix != 1]) > 0) 
		stop("blockMatrix must contain only 0,1 and 2")
	if (length(sirePhasedMatrix[sirePhasedMatrix != 0 & sirePhasedMatrix != 1 & sirePhasedMatrix != 
							9]) > 0) 
		stop("SireMatrix must contain only 0,1 and 9")
	expandMat <- as.numeric(blockMatrix)
	n <- nrow(blockMatrix) * ncol(blockMatrix)
	fMat <- matrix(as.integer(rep(0, n)), ncol = ncol(blockMatrix))
	.C("phaseNogenotype", nrow = as.integer(nrow(blockMatrix)), 
			ncol = as.integer(ncol(blockMatrix)), block = as.integer(blockMatrix), sirePhasedMatrix = as.integer(t(sirePhasedMatrix)), 
			result = fMat)$result
}
pm <- function(blockMatrix, method = "constant")
{
	## if(missing(method))
	##     stop("please set the method")
	METHODS <- c("constant", "relative")
	method <- pmatch(method, METHODS)
	
	if (is.na(method)) 
		stop("invalid method")
	if (method == -1) 
		stop("ambiguous  method")
	
	if (!is.matrix(blockMatrix)) 
		stop("blockMatrix should be a MATRIX")
	expandMat <- as.double(t(blockMatrix))
	n <- ncol(blockMatrix) * nrow(blockMatrix)
	fMat <- matrix(as.double(rep(0, n)), nrow = ncol(blockMatrix))
	res <- .C("pm", expandMat = as.integer(expandMat), nrow = as.integer(ncol(blockMatrix)), 
			ncol = as.integer(nrow(blockMatrix)), method = method,result = fMat)$result
	result <- t(res[-nrow(res), ])
	rownames(result) <- rownames(blockMatrix)
	colnames(result) <- colnames(blockMatrix)[-1]
	result
}


recombinations <- function(blockMatrix)
{
    if (!is.matrix(blockMatrix)) 
        stop("The inputs must be MATRIX")
    if (length(blockMatrix[blockMatrix != 0 & blockMatrix != 1 & blockMatrix != 2]) > 0) 
        stop("Inputs must contain only 0 and 1 or 2")
    mat <- as.numeric(t(blockMatrix))
    nSwitch <- integer(nrow(blockMatrix))
    .C("recombinations", mat = as.integer(mat), nrow = as.integer(nrow(blockMatrix)), ncol = as.integer(ncol(blockMatrix)), 
        result = nSwitch)$result
} 


