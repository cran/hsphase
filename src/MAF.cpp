/*
 * MAF.cpp
 *
 *  Created on: 16/01/2014
 *      Author: mhf
 */

#include "MAF.h"
SEXP MAFC(SEXP snp)
{
	NumericVector SNPs(snp);

	double z = 0, o = 0, t = 0, result = 0;
	for (int i = 0; i < SNPs.length(); i++)
	{
		if (SNPs[i] == 0)
		{
			z = z + 1;
		}
		if (SNPs[i] == 1)
		{
			o = o + 1;
		}
		if (SNPs[i] == 2)
		{
			t = t + 1;
		}
	}
	result = (z * 2 + o) / ((z + o + t) * 2);
	if (result > 0.5)
		result = 1 - result;
	return wrap(result);
}

