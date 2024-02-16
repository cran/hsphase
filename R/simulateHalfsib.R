# Copyright (C) 2014 Cedric Gondro and Mohammad H. Ferdosi  
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

#' Simulate Half-Sibling Genotypes
#'
#' This function simulates genotypes for a set of half-siblings based on specified parameters,
#' including the number of individuals, the number of SNPs, recombination boundaries, and the type of data to return.
#' It generates a sire genotype, maternal half-sib genotypes, and combines these to simulate offspring genotypes,
#' optionally returning phased genotypes based on recombination events.
#'
#' @param numInd Integer, the number of half-siblings to simulate.
#' @param numSNP Integer, the number of SNPs to simulate for each individual.
#' @param recbound Numeric vector, specifying the range of possible recombination events to simulate.
#' @param type Character string, specifying the type of data to return: "genotype" for genotypic data or any other string for phased genotypic data.
#'
#' @return Depending on the \code{type} parameter, this function returns a matrix of simulated genotypic data
#' for half-siblings. If \code{type} is "genotype", it returns unphased genotypic data; otherwise, it returns phased genotypic data.
#'
#' @examples
#' sim_genotypes <- .simulateHalfsib(numInd = 40, numSNP = 10000, recbound = 0:6, type = "genotype")
#' dim(sim_genotypes) # Should return 40 rows (individuals) and 100 columns (SNPs)
#'
#' @export

.simulateHalfsib <- function(numInd = 40, numSNP = 10000, recbound = 0:6, type = "genotype")
{
    sire <- matrix(sample(c(0, 1), numSNP * 2, replace = T), 2, numSNP)
    halfsibsM <- matrix(sample(c(0, 1), numInd * numSNP, replace = T), numInd, numSNP)
    init <- sample(c(1, 2), numInd, replace = T)
    halfsibsP <- matrix(NA, numInd, numSNP)
    recs <- numeric(numInd)
    for (i in 1:numInd)
    {
        strand <- init[i]
        temphap <- sire[strand, ]
        rec <- sample(recbound, 1)
        recs[i] <- rec
        if (rec > 0)
        {
            swap <- sort(sample(100:(numSNP - 100), rec))
            swap2 <- c(swap - 1, numSNP)
            swap <- c(1, swap)
            for (j in 1:length(swap))
            {
                temphap[swap[j]:swap2[j]] <- sire[strand, swap[j]:swap2[j]]
                if (strand == 1) 
                  strand <- 2
                else strand <- 1
            }
        }
        halfsibsP[i, ] <- temphap
    }
    halfsibs <- halfsibsP + halfsibsM
    phasedhalf <- matrix(NA, numInd * 2, numSNP)
    index1 <- seq(from = 1, to = numInd * 2, by = 2)
    index2 <- seq(from = 2, to = numInd * 2, by = 2)
    phasedhalf[index1, ] <- halfsibsP
    phasedhalf[index2, ] <- halfsibsM
    if (type == "genotype") 
        halfsibs
    else phasedhalf
} 
