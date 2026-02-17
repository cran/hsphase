/*
 * fixRotation.cpp
 *
 *  Created on: 04/04/2014
 *      Author: mhf
 */

#include "fixRotation.h"

#define mat imat

using namespace Rcpp;
NumericVector contrast(NumericVector& blockVec)
{
	NumericVector tempforwardVector(blockVec.size(), 20.0);
	for (int i = 0; i < blockVec.size(); i++)
	{
		if (blockVec[i] == 2)
		{
			tempforwardVector[i] = 1;
		}
		if (blockVec[i] == 1)
		{
			tempforwardVector[i] = 2;
		}
		if (blockVec[i] == 0)
		{
			tempforwardVector[i] = 0;
		}

	}
	return (tempforwardVector);
}

int check(NumericVector& A, NumericVector& B)
{
	int temp = 0;
	for (int j = 0; j < A.size(); j++)
	{
		if (A[j] == B[j] && A[j] != 0 && B[j] != 0)
			temp = temp + 1;
	}
	return temp;
}
SEXP fixRotation(SEXP blockStructure)
{
	NumericMatrix bk(blockStructure);
	NumericMatrix res(bk.nrow(), bk.ncol());

	for (int i = 0; i < bk.nrow(); i++)
	{
		for (int j = 0; j < bk.ncol(); j++)
		{
			res(i, j) = bk(i, j);
		}
	}

	NumericVector tempforwardVector(bk.nrow(), 20.0);
	NumericVector MainforwardVector(bk.nrow(), 20.0);
	int checksumMain = 0;
	int checksumCotrast = 0;

	for (int i = 0; i < (bk.ncol() - 1); i++)
	{

				tempforwardVector = res(_, i + 1);
				MainforwardVector = res(_, i);


		checksumMain = check(tempforwardVector, MainforwardVector);
		tempforwardVector = contrast(tempforwardVector);
		checksumCotrast = check(tempforwardVector, MainforwardVector);
		if (checksumMain < checksumCotrast)
		{

			for (int j = 0; j < bk.nrow(); j++)
			{
				if (res(j, i + 1) == 1)
					res(j, i + 1) = 3;
				if (res(j, i + 1) == 2)
					res(j, i + 1) = 1;
				if (res(j, i + 1) == 3)
					res(j, i + 1) = 2;
			}
		}
		/*if (checksumMain == checksumCotrast)
		 {
		 for (int j = 0; j < bk.nrow(); j++)
		 {
		 tempforwardVector = res(_, i);
		 }
		 }*/

	}

	return wrap(res);
}
