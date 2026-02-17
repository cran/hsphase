/*
 * mhMat.cpp
 *
 *  Created on: 02/01/2014
 *      Author: mhf
 */

#include "mhMat.h"
/**
 * Calculate the opposing homozygote matrix
 */
void mhMat::ohGenotype()
{
	arma::mat cm = mainGenotypeMat - (floor(mainGenotypeMat / 9) * 8);
	arma::mat fpart = floor(cm / 2);
	arma::mat lPart = arma::trans(ceil(((cm) - 2) / 2));
	arma::mat result = (fpart * lPart) * (-1);
	ohMat = arma::trans(result) + result;
}
/**
 * Simplify the opposing homozygote matrix the \v ohMat
 * @param maxOH maximum number of oh in each group that is acceptable
 */
void mhMat::ohSimplify(int maxOH)
{
	if (ohMat.is_empty())
	{
		throw "The oh must be calculated!";
	}
	if (ohMatStatus == "simple")
	{
		throw "The oh must not be simplified!";
	}
	for (uint j = 0; j < ohMat.n_cols; j++)
	{
		for (uint i = 0; i < ohMat.n_rows; i++)
		{
			if (ohMat(i, j) > maxOH)
			{
				ohMat(i, j) = 10;
			}
			else
			{
				ohMat(i, j) = 0;
			}
		}
	}
	ohMatStatus = "simple";
}

/**
 * Group the individuals based on the simplified \v  maxOH
 */
void mhMat::validateohSimplify()
{
	if (ohMatStatus != "simple")
	{
		throw "The ohMat must be simplified first";
	}

	resultContainer result(ohMat.n_rows);

	vector<int>::iterator iElemnt;
	int lastGroup = 1;
	bool newGroup = FALSE;
	for (uint j = 0; j < ohMat.n_cols; j++)
	{
		for (uint i = 0; i < ohMat.n_rows; i++)
		{
			if (ohMat(i, j) == 0 && result.getCondition(i) == FALSE && result.getIndex(i) != 0 /*&& i != j*/)
			{
				result.setCondition(TRUE, i);
				result.setIndex(lastGroup, i);
				newGroup = TRUE;
			}
			// Check other columns of ohMat for inconsistency
			if (ohMat(i, j) != 0 && result.getCondition(i) == TRUE && result.getIndex(i) == result.getIndex(j))
			{
				result.setIndex(0, i);
				result.setIndex(0, j);

				/* Rprintf("i = %d", i);
				 Rprintf("j = %d", j);
				 */
			}
		}
		if (newGroup)
		{
			lastGroup = lastGroup + 1;
			newGroup = FALSE;
		}
	}
	for (uint i = 0; i < ohMat.n_rows; i++)
	{
		groupResult.push_back(result.getIndex(i));
	}
}

void mhMat::keep2MostFrequent()
{
	CounterMap counts;
	for (int j = 0; j < groupResult.size(); ++j)
	{
		CounterMap::iterator i(counts.find(groupResult[j]));
		if (i != counts.end())
		{
			i->second++;
		}
		else
		{
			counts[groupResult[j]] = 1;
		}
	}

	saveMap sorted;
	for (CounterMap::iterator it = counts.begin(); it != counts.end(); ++it)
	{

		sorted.insert(make_pair(it->second, it->first));
	}

	saveMap::iterator itr = sorted.end();
	saveMap::iterator itr2 = sorted.end();
	--itr;
	--itr2;
	--itr2;

	if (itr->second == 0)
	{
		--itr;
		--itr2;
	}
	if (itr2->second == 0)
	{
		--itr2;
	}
	//cout << itr->first << " " << itr2->first << endl;

	for (vector<int>::iterator it = groupResult.begin(); it != groupResult.end(); ++it)
	{

		if (sorted.size() > 1)
		{
			if (itr->first > 1 && itr2->first > 1)
			{
				if (*it != itr->second && *it != itr2->second)
				{
					*it = 0;
				}
				else if (*it == itr->second)
				{
					*it = 1;
				}
				else if (*it == itr2->second)
				{
					*it = 2;
				}
			}
			else if (itr->first > 1)
			{
				if (*it != itr->second)
				{
					*it = 0;
				}
				else if (*it == itr->second)
				{
					*it = 1;
				}
			}
			else
			{
				*it = 0;
			}

		}
		else
		{
			*it = 1;
		}
	}
}
