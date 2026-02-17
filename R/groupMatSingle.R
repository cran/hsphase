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
groupMatSingle <- function(haplotype, windowsSize, cpus = 2, input = "haplotype", oh = 0)
{
    if (input == "haplotype")
    {
        haplotype[haplotype == 1] <- 2
        result <- hsphase::.ibdCluster(haplotype, cpus, windowsSize, oh)
        result <- do.call(cbind, result)
        hsphase::.fixRotation(hsphase::.fixBothStrand(result))
        
    } else if (input == "genotype")
    {
        result <- hsphase::.ibdCluster(haplotype, cpus, windowsSize, oh)
        result <- do.call(cbind, result)
        hsphase::.fixRotation(result)
    } else
    {
        print("Please select the input type")
    }
} 
