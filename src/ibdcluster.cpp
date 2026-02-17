/*

 * ibdcluster.cpp
 *
 *  Created on: 30/12/2013
 *      Author: mhf
 */

#include "ibdcluster.h"
using namespace Rcpp;
using namespace std;

SEXP ibdCluster(SEXP geno, SEXP threads, SEXP windowsSize,SEXP oh)
{
	NumericMatrix genotype(geno);
	int mt = as<int>(threads);
  int maxoh = as<int>(oh);
	int winSize = as<int>(windowsSize);
	int n = genotype.nrow(), k = genotype.ncol();
	arma::mat Genotype(as<arma::mat>(genotype).begin(), n, k, TRUE);
	/*mhMat matGenotypeA;
	 vector<int> result = matGenotypeA.groupWrapCommand(Genotype.cols(0, winSize), maxoh);
	 */

	vector<vector<int> > ohResult;
	ohResult.resize(k - winSize);

//if( windowsSize*n > 5000)
  Rcout << "Number of threads: " << mt << endl;
  // num_threads(mt)
  #pragma omp parallel for
	for (int i = 0; i < k - winSize; i++)
	{
		mhMat matGenotypeA;
		ohResult.at(i) = matGenotypeA.groupWrapCommand(Genotype.cols(i, i + winSize), maxoh);
	}

	return (wrap(ohResult));

}

