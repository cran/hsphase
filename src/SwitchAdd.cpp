// Copyright (C) 2014 Mohammad H. Ferdosi
//
// HSPhase is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// HSPhase program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
/*
 * SwitchAdd.cpp
 *
 *  Created on: 09/06/2014
 *      Author: mhf
 */

#include "SwitchAdd.h"
#include <algorithm>   // std::swap

SEXP switchAdd(SEXP haplotype, SEXP switchpoints, SEXP minlength)
{
    NumericMatrix Haplotypemain(haplotype);
    const int nrow    = Haplotypemain.nrow();
    const int ncol    = Haplotypemain.ncol();
    const int minLength = as<int>(minlength);

    List SwitchPoints(switchpoints);
    const int nInd = SwitchPoints.size();

  
    if (2 * nInd > nrow)
        stop("switchAdd: number of switch-point vectors (%d) exceeds nrow/2 (%d).",
             nInd, nrow / 2);

   
    std::vector<std::vector<int>> Haplotype(nrow, std::vector<int>(ncol));
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < ncol; j++)
        {
            const double v = Haplotypemain(i, j);
            Haplotype[i][j] = ISNAN(v) ? NA_INTEGER : static_cast<int>(v);
        }
    }

  
    std::vector<std::vector<int>> allInd(nInd);
    for (int i = 0; i < nInd; i++)
        allInd[i] = as<std::vector<int>>(SwitchPoints[i]);

    #pragma omp parallel for
    for (int i = 0; i < nInd; i++)
    {
        const std::vector<int>& ind = allInd[i];
        std::vector<int>& h1 = Haplotype[2 * i];
        std::vector<int>& h2 = Haplotype[2 * i + 1];

   
        for (size_t j = 1; j + 1 < ind.size(); j++)
        {
            const int p    = ind[j];
            const int prev = ind[j - 1];
            const int next = ind[j + 1];


            if (p < 0 || p >= ncol)
                continue;

            if ((p - prev) > minLength && (next - p) > minLength)
            {
                for (int k = p; k < ncol; k++)
                    std::swap(h1[k], h2[k]);
            }
        }
    }

  
    NumericMatrix output(nrow, ncol);
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < ncol; j++)
        {
            output(i, j) = (Haplotype[i][j] == NA_INTEGER)
                               ? NA_REAL
                               : static_cast<double>(Haplotype[i][j]);
        }
    }

    return wrap(output);
}
