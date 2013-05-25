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

	int t = as<int> (threads);
	NumericMatrix genotype(geno);
	int n = genotype.nrow(), k = genotype.ncol();
	arma::mat Genotype(as<arma::mat>(genotype).begin(), n, k, TRUE);

	arma::mat result = arma::zeros<arma::mat>(genotype.nrow(), genotype.nrow());
	arma::mat res2 = arma::zeros<arma::mat>(genotype.nrow(), genotype.nrow());
	// make_nas(Genotype);

	#pragma omp parallel for private(res2) schedule(dynamic)  num_threads(t)
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
