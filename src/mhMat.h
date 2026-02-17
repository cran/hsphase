/*
 * mhMat.h
 *
 *  Created on: 02/01/2014
 *      Author: mhf
 */

#ifndef MHMAT_H_
#define MHMAT_H_
#include <math.h>
#include <RcppArmadillo.h>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <map>

#define uint unsigned int
typedef std::multimap<int, unsigned int> saveMap;
typedef std::unordered_map<int, int> CounterMap;

using namespace std;
class mhMat
{
public:

	void ohGenotype();
	void ohSimplify(int maxOH);
	void validateohSimplify();
	/**
	 * Group the subMatrix of genotype
	 * @param subGenotype
	 * @param maxoh
	 * @return the results will be saved in  \v groupResult variable
	 */
	const vector<int>& groupWrapCommand(arma::mat subGenotype, int maxoh)
	{
		setMainMat(subGenotype);
		ohGenotype();
		ohSimplify(maxoh);
		validateohSimplify();
		keep2MostFrequent();
        return groupResult;
	}
	mhMat() :
			ohMatStatus("NA")
	{
	}

	const arma::mat& getMainMat() const
	{
		return mainGenotypeMat;
	}

	void setMainMat(const arma::mat& mainMat)
	{
		this->mainGenotypeMat = mainMat;
	}

	const arma::mat& getOhMat() const
	{
		return ohMat;
	}

	void setOhMat(const arma::mat& ohMat)
	{
		this->ohMat = ohMat;
	}

	const vector<int>& getGroupResult() const
	{
		return groupResult;
	}

	const string& getOhMatStatus() const
	{
		return ohMatStatus;
	}

private:
	arma::mat mainGenotypeMat;
	arma::mat ohMat;
	vector<int> groupResult;
	string ohMatStatus;
	void keep2MostFrequent();
};

class resultContainer
{
public:
	resultContainer(int size)
	{
		itsSize_ = size;
		index_ = new int[size];
		condition_ = new bool[size];
		for (int i = 0; i < size; i++)
		{
			condition_[i] = FALSE;
			index_[i] = 20;
		}
	}
	virtual ~resultContainer()
	{
		delete[] index_;
		delete[] condition_;
	}
	bool getCondition(int i) const
	{
		return condition_[i];
	}

	void setCondition(bool cond, int i)
	{
		*(condition_ + i) = cond;
	}

	int getIndex(int i) const
	{
		return index_[i];
	}

	void setIndex(int indexA, int i)
	{
		*(index_ + i) = indexA;
	}

private:
	int* index_;
	bool* condition_;
	int itsSize_;

};


#endif /* MHMAT_H_ */
