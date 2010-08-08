// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#include <cassert>

#include "pvalue.hpp"

using namespace std;
using namespace Permory;

double Permory::quisq(double t, double df)
{
    return 1 - gsl_cdf_chisq_P(t, df);
}

std::vector<double> Permory::single_step_pvalues(
        const std::vector<double>& t, 
        const std::vector<double>& tperm, 
        bool cc) 
{
    if (!vector_is_sorted(t, true)) { //sorted in decreasing order?
        cerr << "Bad values in function 'Permory::single_step_pvalues': " <<
            "the test statistics are not monotonically decreasing!\n";
        exit(-1);
    }

    std::vector<double> p(t.size(), 0);
    for (size_t i=0; i<tperm.size(); ++i) {     //for each permutation
        for (size_t j=0; j<t.size(); ++j) {     //for each test statistic
            if (tperm[i] >= t[j]) p[j]++;
        }
    }
    if (!cc) {
        for (size_t i=0; i<p.size(); ++i)
            p[i] = (p[i] + 1)/(double(tperm.size() + 1)); //p-values from counts
    }
    return p;
}

std::vector<double> Permory::step_down_pvalues(
        const std::vector<double>& t, 
        const matrix_type& m, 
        const std::vector<double>& v, 
        bool cc)       
{
    if (!vector_is_sorted(t, true)) { //sorted in decreasing order?
        cerr << "Bad values in function 'Permory::single_step_pvalues': " <<
            "the test statistics are not monotonically decreasing!\n";
        exit(-1);
    }
    assert (t.size() == m.size());
    size_t X = m.size();
    size_t nn = m.front().size(); //number of permutations

    // The "step-down matrix" 'm' contains the null distribution obtained via 
    // permutation for the top X markers (X == #rows of m):
    //           |  1   2  ... nPerm 
    //         ----------------------
    //         1 |  .   .  ...   .    
    //         2 |  .   .  ...   .    
    //         . |  .   .  ...   .    
    //         X |  .   .  ...   .   
    //
    // The rest of the null distr is condensed in 'v'.
    // In order to determine step-down adjusted p-values, we examine the test 
    // statistics for each permutation step, i.e., we go column by column 
    // through this matrix.

    std::vector<double> p(t.size(), 0);
    for (size_t j=0; j<nn; ++j) {   //for each permutation
        double d = v[j];            //for now consider only permutation j
        // In each column determine successive maxima: 
        // Let k = X and start with rank k, followed by k-1, k-2, ..., 2, 1. 
        // That is, in the 1st step we only consider the test stat belonging 
        // to rank X, in the 2nd step the one belonging to ranks X and X-1, 
        // and so on ...
        size_t k = X;
        while (--k > 0) {
            d = max(d, m[k][j]);
            if(d >= t[k]) 
                p[k]++;
        }
    }
    // Enforce monotonicity using succesive maximization
    for(size_t i=1; i<p.size(); ++i) 
        p[i] = max(p[i-1], p[i]);

    if (!cc) {
        for (size_t i=0; i<p.size(); ++i)
            p[i] = (p[i] + 1)/(double(nn + 1)); //p-values from counts
    }
    return p;
}
