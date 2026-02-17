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
# fix the switch errors
addSwitch <- function(haplotypeMatrix, switchPoints, minLength)
{
 
  ID = rownames(haplotypeMatrix)
  if(nrow(haplotypeMatrix)/2!=length(switchPoints))    
  {
    stop("The number of elements in the switchPoints must be equal to the number of individuals")
  }
  
  for (i in 1:length(switchPoints))
    {
        switchPoints[[i]] <- c(0, switchPoints[[i]], ncol(haplotypeMatrix))
    }
    result <- .Call("switchAdd", haplotypeMatrix, switchPoints, minLength, PACKAGE = "hsphase")
    rownames(result) = ID
    result
}

