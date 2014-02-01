/*
 * ohp.cpp
 *
 *  Created on: 07/09/2013
 *      Author: mhf
 */

#include "ohp.h"
#include <iostream>
#include <R.h>

int ohpFunction(int const * matrix, const int* nrow, const int* ncol, int* result)
{
	//For matrix

	int const**pRowsMat = new int const*[*nrow];

	for (int i = 0; i < *nrow; i++)
	{
		pRowsMat[i] = matrix + (*ncol) * i;
	}

	// For result
	int **pRowsRes = new int*[*nrow];

	for (int i = 0; i < *nrow; i++)
	{
		pRowsRes[i] = result + (*nrow) * i;
	}

	int frq = 0;
    #pragma omp parallel for private(frq) schedule(dynamic) num_threads(2)
	for (int i = 0; i < *nrow; i++)
	{
		for (int j = i; j < *nrow; j++)
		{
			for (int k = 0; k < *ncol; k++)
			{
				if ((pRowsMat[i][k] == 2 && pRowsMat[j][k] == 0) || (pRowsMat[i][k] == 0 && pRowsMat[j][k] == 2))
				{
					frq = frq + 1;
				}
			}
			pRowsRes[i][j] = frq;
			frq = 0;
		}
	}


	delete[] pRowsMat;
	delete[] pRowsRes;
}

extern "C"
{

void ohp(int const * matrix, int *nrow, int *ncol, int* result)
{

	ohpFunction(matrix, nrow, ncol, result);

}
}

