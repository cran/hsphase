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
#include "ohg.h"

#define mat imat

using namespace Rcpp;
int make_nas(arma::mat & genotype)
{
	for (unsigned int i = 0; i < genotype.n_rows; i++)
	{
		for (unsigned int j = 0; j < genotype.n_cols; j++)
		{
			if (genotype.at(i, j) == 1)
			{
				genotype.at(i, j) = 9;
			}
		}

	}
	return 0;
}

arma::mat twoFreqs(const arma::mat & data)
{
	arma::mat result = arma::zeros<arma::mat>(1, data.n_rows);
	int frequency = 0;
	for (unsigned int i = 0; i < data.n_rows; i++)
	{
		frequency = 0;
		for (unsigned int j = 0; j < data.n_cols; j++)
		{
			if (data.at(i, j) == 2 || data.at(i, j) == -2)
			{
				frequency = frequency + 1;
			}
		}
		result.at(0, i) = frequency;
	}
	return result;

}

arma::mat vecMinusMats(const arma::mat &genotype, const int index)
{
	arma::mat result = genotype;
	arma::mat a = genotype.row(index);

	for (unsigned int i = index; i < genotype.n_rows; i++)
	{
		result.row(i) = a - genotype.row(i);
	}
	return result;
}
SEXP ohg(SEXP geno, SEXP threads)
{

//  int t = as<int> (threads);
	NumericMatrix genotype(geno);
	int n = genotype.nrow(), k = genotype.ncol();
	arma::mat Genotype(as<arma::mat>(genotype).begin(), n, k, TRUE);

	arma::mat result = arma::zeros<arma::mat>(genotype.nrow(), genotype.nrow());
	arma::mat res2 = arma::zeros<arma::mat>(genotype.nrow(), genotype.nrow());
	// make_nas(Genotype);

	#pragma omp parallel for private(res2) schedule(dynamic)  num_threads(as<int> (threads))
	for (int i = 0; i < genotype.nrow(); i++)
	{
	//	Rprintf("%d\n", i);
		res2 = vecMinusMats(Genotype, i);
		result.row(i) = twoFreqs(res2);
		res2.clear();
	}
	List res;
	res["mat+mat"] = result;
	return res;
}
