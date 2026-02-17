/*
 * FixBothStrand.cpp
 *
 *  Created on: 11/04/2014
 *      Author: mhf
 */

#include "FixBothStrand.h"

SEXP fixBothStrand(SEXP groupsR)
{
	NumericMatrix groups(groupsR);
	int n = groups.nrow(), k = groups.ncol();
	arma::mat Groups(as<arma::mat>(groups).begin(), n, k, TRUE);
	arma::mat result = arma::zeros<arma::mat>(n, k);

	for (int i = 0; i < n ; i = i + 2)
	{
		for (int j = 0; j < k; j++)
		{
		if (Groups(i, j) != 0 && Groups(i + 1, j) != 0)
			{
				result(i, j) = 100;
				result(i + 1, j) = 100;
			}
			else
			{
				result(i, j) = Groups(i, j);
				result(i + 1, j) = Groups(i + 1, j);
			}
		}
	}
	return wrap(result);
}
