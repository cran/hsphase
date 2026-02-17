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
##' Convert a two-row haplotype per individual to a one-row representation
##'
##' Converts a haplotype matrix where each individual is represented by two rows
##' (allele 1 and allele 2) into a single-row-per-individual representation.
##'
##' @param haplotype A matrix containing haplotypes with two rows per individual.
##' @return A matrix with one row per individual (allele 1 and allele 2 combined).
.ptr2por <- function(haplotype)
{
    t( gdata::interleave(t(haplotype[seq(from = 1, to = nrow(haplotype), by = 2), ]), t(haplotype[seq(from = 2, to = nrow(haplotype), by = 2), ])))
}
