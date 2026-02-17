// Copyright (C) 2014 Mohammad H. Ferdosi
//
// HSPhase is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// HSPhase program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
/*
 * SwitchDetector.cpp
 *
 *  Created on: 10/06/2014
 *      Author: mhf
 */

#include "SwitchDetector.h"
using namespace std;
using namespace Rcpp;

SEXP switchDetector(SEXP group_Matrix)
{
	NumericMatrix groupmatrix(group_Matrix);
	arma::mat groupMatrix(as<arma::mat>(groupmatrix).begin(), groupmatrix.nrow(), groupmatrix.ncol(), TRUE);
	groupMatrix.fill(8);

	HaplotypeStrand *myHaplotype = new HaplotypeStrand[groupmatrix.nrow() / 2];
	vector<vector<int> > switchResult;

	for (int i = 0; i < groupmatrix.nrow() / 2; i++)
	{
		myHaplotype[i].solveAll(100, i, groupmatrix);
		switchResult.push_back(myHaplotype[i].getSwitchLoci());
	}
	return wrap(switchResult);
}

HaplotypeStrand::HaplotypeStrand()
{
	Strand = 0;
	strand1 = 0;
	strand2 = 0;
	AddedStrand = 0;
	set3Value = 0;
	ncol = 0;
}

HaplotypeStrand::~HaplotypeStrand()
{
	delete[] strand1;
	delete[] strand2;
	delete[] AddedStrand;
}

int HaplotypeStrand::switchAddPoint(int point)
{
	int *temp1 = new int[ncol];
	int *temp2 = new int[ncol];

	for (int i = 0; i < ncol; i++)
	{
		temp1[i] = Strand[0][i];
		temp2[i] = Strand[1][i];
	}
	for (int i = point; i < ncol; i++)
	{
		temp1[i] = Strand[1][i];
		temp2[i] = Strand[0][i];
	}
	for (int i = point; i < ncol; i++)
	{
		Strand[0][i] = temp1[i];
		Strand[1][i] = temp2[i];
	}
	return 0;
}

int HaplotypeStrand::fillGap3()
{
	for (int i = 0; i < ncol; i++)
	{
		if (strand1[i] == set3Value)
		{
			int gapEnd = oppos3(i, strand1);
			if (gapEnd == strand1[i - 1])
			{
				fillGap3(i, strand1, gapEnd);
			}
		}
	}
	for (int i = 0; i < ncol; i++)
	{
		if (strand2[i] == set3Value)
		{
			int gapEnd = oppos3(i, strand2);
			if (gapEnd == strand2[i - 1])
			{
				fillGap3(i, strand2, gapEnd);
			}
		}
	}

	return 0;
}

int HaplotypeStrand::firstOccurance(int start, int value, int* strandA)
{
	for (int i = start; i < ncol; i++)
	{
		if (strandA[i] == value)
		{
			return i;
		}
	}
	return -1;
}

int HaplotypeStrand::min(int n1, int n2, int n3)
{
	if (n1 < n2 && n1 < n3)
	{
		return n1;
	}
	if (n2 < n1 && n2 < n3)
	{
		return n2;
	}
	if (n3 < n2 && n3 < n1)
	{
		return n3;
	}
	return 0;
}

int HaplotypeStrand::oppos3(int start, int* strand)
{
	for (int i = start; i < ncol; i++)
	{
		if (strand[i] != 3)
			return strand[i];
	}
	return 0;
}

int HaplotypeStrand::fillGap3(int start, int*strand, int fillValue)
{
	for (int i = start; i < ncol; i++)
	{
		if (strand[i] == set3Value)
		{
			strand[i] = fillValue;
		}
		else
		{
			break;
		}
	}
	return 0;
}

int HaplotypeStrand::replace(int oldValue, int newValue)
{
	for (int i = 0; i < ncol; i++)
	{
		if (strand1[i] == oldValue)
		{
			strand1[i] = newValue;
		}
		if (strand2[i] == oldValue)
		{
			strand2[i] = newValue;
		}
	}
	return 0;
}

void HaplotypeStrand::replaceSTR1(int oldValue, int newValue)
{
	for (int i = 0; i < ncol; i++)
	{
		if (strand1[i] == oldValue)
		{
			strand1[i] = newValue;
		}
	}
}

void HaplotypeStrand::replaceSTR2(int oldValue, int newValue)
{
	for (int i = 0; i < ncol; i++)
	{
		if (strand2[i] == oldValue)
		{
			strand2[i] = newValue;
		}
	}
}
void HaplotypeStrand::addStrands()
{
	AddedStrand = new int[ncol];
	for (int i = 0; i < ncol; i++)
	{
		AddedStrand[i] = strand1[i] + strand2[i];
	}
}

void HaplotypeStrand::switchIdentification()
{
	int CurrentStrand = 0;
	int previousIndex = 0;
	int initialI0 = 0;

	int oneOC = firstOccurance(0, 1, AddedStrand);
	int twoOC = firstOccurance(0, 2, AddedStrand);

	if (oneOC > twoOC && oneOC != -1 && twoOC != -1)
	{
		CurrentStrand = 2;
		initialI0 = twoOC;
	}
	else if (oneOC < twoOC && oneOC != -1 && twoOC != -1)
	{
		CurrentStrand = 1;
		initialI0 = oneOC;
	}
	else if (twoOC == -1)
	{
		CurrentStrand = 1;
		initialI0 = oneOC;
	}
	else if (oneOC == -1)
	{
		CurrentStrand = 2;
		initialI0 = twoOC;
	}

	else
	{
		Rprintf("Error %d %d\n", oneOC, twoOC);
	}

	for (int i = initialI0; i < ncol; i++)
	{
		if (AddedStrand[i] == 0 || AddedStrand[i] != CurrentStrand)
		{
			previousIndex = i;
		}
		if (AddedStrand[i] != CurrentStrand && AddedStrand[i] != 0)
		{
//			 condition ? value_if_true : value_if_false
			CurrentStrand == 2 ? CurrentStrand = 1 : CurrentStrand = 2;
			if (previousIndex != initialI0)
			{
//				Rprintf("--%d  %d\n", previousIndex, initialI0);
				switchLoci.push_back(previousIndex);
			}
		}
	}
	if (switchLoci.size() == 0)
	{
		switchLoci.push_back(0);
	}
}

/**
 * @param value3 common IBD in the haplotypes
 * @param i for loop i
 * @param groupmatrix input numeric matrix
 */
void HaplotypeStrand::solveAll(int value, int i, NumericMatrix & groupmatrix)
{
	setSet3Value(100);
	setStrand(groupmatrix, i);
	fillGap3();
	replace(100, 0);
	replace(2, 1);
	replaceSTR2(1, 2);
	addStrands();
	switchIdentification();
}
