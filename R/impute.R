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

impute <- function(halfsib_parental_ld, halfsib_block_ld, sire_ld, sire_hd)
{
	snpld <- colnames(halfsib_parental_ld)
	snphd <- colnames(sire_hd)
	ldindex <- match(snpld,snphd)
	nald <- which(is.na(ldindex))
	if(length(nald)>0)
	{
		ldindex <- ldindex[-nald]
		halfsib_parental_ld <- halfsib_parental_ld[,-nald]
		halfsib_block_ld <- halfsib_block_ld[,-nald]
		sire_ld <- sire_ld[,-nald]
	}
	
	
	sire_hdformLD <- sire_hd[,ldindex]
	sire_ld[sire_ld==9] <- sire_hdformLD[sire_ld==9]
	nineIndex <- which(sire_hdformLD[1,]!=9)
	a <- summary(lm(sire_ld[1,nineIndex]~sire_hdformLD[1,nineIndex]))$r.squared
	b <- summary(lm(sire_ld[1,nineIndex]~sire_hdformLD[2,nineIndex]))$r.squared
	
	if(a<b)
	{
		sire_hd <- sire_hd[c(2,1),]
	}
	
	result <- matrix(0,ncol=ncol(sire_hd),nrow=nrow(halfsib_parental_ld))
	result[,ldindex] <- halfsib_block_ld
	result[result==1] <- 3
	result[result==2] <- 4
	result <- .hblock(result)
	.phfnoGenotype(result,sire_hd)
}



