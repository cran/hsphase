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
 * SwitchDetector.h
 *
 *  Created on: 10/06/2014
 *      Author: mhf
 */

#ifndef SWITCHDETECTOR_H_
#define SWITCHDETECTOR_H_

#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;
RcppExport SEXP switchDetector(SEXP groupMatrix);

class HaplotypeStrand
{
public:
	HaplotypeStrand();
	~HaplotypeStrand();

	int** getStrand() const
	{
		return Strand;
	}
	void setStrand(NumericMatrix matrix, int haploNum)
	{
		Strand = new int*[2];
		for (int i = 0; i < 2; i++)
			Strand[i] = new int[matrix.ncol()];

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < matrix.ncol(); j++)
			{
				Strand[i][j] = matrix(i + (haploNum * 2), j);
			}
		}
		this->ncol = matrix.ncol();
		strand1 = Strand[0];
		strand2 = Strand[1];
	}

	int switchAddPoint(int point);
	int fillGap3();
	int firstOccurance(int strat, int value, int* strand);
	int min(int n1, int n2, int n3);
	int oppos3(int start, int* strand);
	int fillGap3(int start, int*strand, int fillValue);
	int replace(int oldValue, int newValue);
	void replaceSTR1(int oldValue, int newValue);
	void replaceSTR2(int oldValue, int newValue);
	void addStrands();
	void switchIdentification();

	/**
	* @param value3 common IBD in the haplotypes
	* @param i for loop i
	* @param groupmatrix input numeric matrix
	*/
	void solveAll(int value3, int i, NumericMatrix & groupmatrix);

	const std::vector<int>& getSwitchLoci() const
	{
		return switchLoci;
	}

	void setSet3Value(int set3Value)
	{
		this->set3Value = set3Value;
	}

private:
	int** Strand;
	int *strand1;
	int *strand2;
	int *AddedStrand;
	int ncol;
	int set3Value;
	std::vector<int> switchLoci;

};

#endif /* SWITCHDETECTOR_H_ */
