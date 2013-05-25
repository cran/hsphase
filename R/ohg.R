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

ohg <- function(genotypeMatrix, cpus = 2)
{
  if (is.null(genotypeMatrix)) 
    stop("Invalid input!")
  if (!is.matrix(genotypeMatrix)) 
    stop("Genotype should be a MATRIX")
  if (length(genotypeMatrix[genotypeMatrix != 9 & genotypeMatrix != 0 & genotypeMatrix != 
                              2 & genotypeMatrix != 1]) > 0) 
    stop("Genotype must contain only 0, 1, 2 or 9")
  
  result <- .Call( "ohg", genotypeMatrix, cpus, PACKAGE = "hsphase" )
  rownames(result[[1]]) <- rownames(genotypeMatrix)
  colnames(result[[1]]) <- rownames(genotypeMatrix)
  result[[1]][lower.tri(result[[1]])] <- t(result[[1]])[lower.tri(result[[1]])]
  result[[1]]
}
