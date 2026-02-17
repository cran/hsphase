/*
 * FillGap.cpp
 *
 *  Created on: 26 Sep 2014
 *      Author: mhf
 */

#include "FillGap.h"
int fillGapFunction(int *vector, int * length, int *result)
{
	/*for (int i = 0; i < *length; i++)
	 {
	 std::cout << vector[i] << "\n";
	 }*/

	for (int i = 0; i < *length; i++)
	{
		result[i] = vector[i];
	}

	int gapLength = 0;
	if (result[0] == 0)
	{
		for (size_t j = 0; j < *length; j++)
		{
			if (result[j] != 0)
				break;
			gapLength += 1;
		}

		int secondVal = result[gapLength];
		if (gapLength != *length)
			for (size_t j = 0; j < gapLength; j++)
			{
				result[j] = secondVal;
			}
	}

	gapLength = 0;
	if (result[*length - 1] == 0)
	{
		for (size_t j = *length - 1; j > 0; j--)
		{
			if (result[j] != 0)
				break;
			gapLength += 1;
		}
//		std::cout << gapLength << " " << *length << std::endl;
		int secondVal = result[*length - gapLength - 1];
		if (gapLength != *length - 1)
			for (size_t j = *length - 1; j > (*length-gapLength - 1); j--)
			{
				result[j] = secondVal;
			}
	}

	gapLength = 0;
	bool gap = false;
	for (int i = 1; i < (*length - 1); i++)
	{

		if (result[i] != 0)
		{
			gap = false;
		}
		else
		{
			gap = true;
		}
		if (gap)
		{
			for (size_t j = i; j < (*length - 1); j++)
			{
				if (result[j] != 0)
					break;
				gapLength += 1;

			}

			int firstVal = result[i - 1];
			for (size_t j = i; j < i + (gapLength / 2); j++)
			{
				result[j] = firstVal;
			}

			int secondVal = result[i + gapLength];
			for (size_t j = i + (gapLength / 2); j < i + gapLength; j++)
			{
				result[j] = secondVal;
			}
			gapLength = 0;
		}
	}

	return (0);

}

