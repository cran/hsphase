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

ohg <- function(GenotypeMatrix)
{
	if (length(GenotypeMatrix[GenotypeMatrix != 9 & GenotypeMatrix != 0 & GenotypeMatrix != 2 & GenotypeMatrix != 
							1]) > 0) 
		stop("GenotypeMatrix must contain only 0, 1, 2 or 9")
	n <- nrow(GenotypeMatrix) * nrow(GenotypeMatrix)
	fMat <- matrix(as.integer(rep(0, n)), nrow = nrow(GenotypeMatrix))
	expandMat <- as.numeric(t(GenotypeMatrix))
	result <- .C("ohp", expandMat = as.integer(expandMat), nrow = as.integer(nrow(GenotypeMatrix)), 
			ncol = as.integer(ncol(GenotypeMatrix)), result = fMat)$result
	result[upper.tri(result)] <- t(result)[upper.tri(result)]
	rownames(result) <- colnames(result)  <- rownames(GenotypeMatrix)
	result
}
