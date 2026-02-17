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

SEXP switchAdd(SEXP haplotype, SEXP switchpoints, SEXP minlength)
{
    NumericMatrix Haplotypemain(haplotype);
    int nrow = Haplotypemain.nrow();
    int ncol = Haplotypemain.ncol();
    int minLength = as<int>(minlength);

    std::vector<std::vector<int>> Haplotype(nrow, std::vector<int>(ncol));
    std::vector<std::vector<int>> newHaplotype(nrow, std::vector<int>(ncol));


    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < ncol; j++)
        {
            Haplotype[i][j] = static_cast<int>(Haplotypemain(i, j));
        }
    }

    List SwitchPoints(switchpoints);

//#pragma omp parallel for
    for (int i = 0; i < SwitchPoints.size(); i++)
    {
        std::vector<int> ind = SwitchPoints[i];
        #pragma omp parallel for
        for (size_t j = 1; j < ind.size() - 1; j++)
        {
            //if ((ind[j] - ind[j - 1]) > minLength && (ind[j + 1] - ind[j] > minLength))
             //{
                std::vector<int> haplotype1 = Haplotype[i * 2];
                std::vector<int> haplotype2 = Haplotype[i * 2 + 1];

                for (int k = ind[j]; k < ncol; k++)
                {
                    haplotype1[k] = Haplotype[i * 2 + 1][k];
                    haplotype2[k] = Haplotype[i * 2][k];
                }

                Haplotype[i * 2] = haplotype1;
                Haplotype[i * 2 + 1] = haplotype2;
            //}
        }

        newHaplotype[i * 2] = Haplotype[i * 2];
        newHaplotype[i * 2 + 1] = Haplotype[i * 2 + 1];
    }

    NumericMatrix output(nrow, ncol);


    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < ncol; j++)
        {
            output(i, j) = newHaplotype[i][j];// ? 1.0 : 0.0;
        }
    }

    return wrap(output);
}

