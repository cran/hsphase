/*

 * ibdcluster.h
 *
 *  Created on: 30/12/2013
 *      Author: mhf
*/


#ifndef IBDCLUSTER_H_
#define IBDCLUSTER_H_

#include <RcppArmadillo.h>
#include "mhMat.h"

using namespace Rcpp;
RcppExport SEXP ibdCluster(SEXP genotype, SEXP maxOH, SEXP windowsSize, SEXP oh) ;
/**
 * Compare two integer vector and change them to the new vectors based on their common elements
 * @param first output of groupWrapCommand in mhMat class
 * @param second output of groupWrapCommand in mhMat class
 * @return
 */
int compareTwoGroups(vector<int>& first, vector<int>& second);
#endif
