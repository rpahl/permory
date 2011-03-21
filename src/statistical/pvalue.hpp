// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_pvalue_hpp
#define permory_pvalue_hpp

#include <vector>
#include <deque>
#include <gsl/gsl_cdf.h>

#include "detail/config.hpp"
#include "detail/vector.hpp"

namespace Permory { namespace statistic {
    typedef std::vector<std::vector<double> > matrix_type;

    // Computes single step counts
    std::deque<size_t> single_step_counts(
            const std::deque<double>& t,        //test statistics
            std::deque<double>* tperm)          //max test statistic per permutation
    {
        using namespace std;
        sort(tperm->begin(), tperm->end());
        deque<size_t> cnts(t.size(), 0);
        for (size_t j=0; j<t.size(); ++j) {             //for each test statistic
            deque<double>::iterator it = lower_bound(tperm->begin(), tperm->end(), t[j]);
            cnts[j] = tperm->end() - it;
        }
        return cnts;
    }

    // Computes single step p-values 
    std::deque<double> single_step_pvalues(
            const std::deque<size_t>& counts,   //see function above
            size_t nperm)                       //number of permutations
    {
        using namespace std;
        deque<double> p(counts.size());
        for (size_t i=0; i<p.size(); ++i) {
            p[i] = double(counts[i] + 1)/(nperm + 1); 
        }
        return p;
    }

    // Computes step down counts
    std::deque<size_t> step_down_counts(
            const std::deque<double>& t,//top X test stats sorted _decreasingly_
            const matrix_type& m,//max teststats per permutation for top X markers
            const std::deque<double>& v)//condensed max test stats of the non-top markers
    {
        using namespace std;
        if (not detail::sequence_is_sorted(t, true)) { //sorted in decreasing order?
            throw invalid_argument("Test statistics must be monotonically decreasing.");
        }
        assert (t.size() == m.size());
        size_t X = m.size();
        size_t nperm = m.front().size(); //number of permutations

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

        //deque<double> p(t.size(), 0);
        deque<size_t> cc(t.size(), 0);
        for (size_t j=0; j<nperm; ++j) {   //for each permutation
            double d = v[j];            //for now consider only permutation j
            // In each column determine successive maxima: 
            // Let k = X and start with rank k, followed by k-1, k-2, ..., 2, 1. 
            // That is, in the 1st step we only consider the test stat belonging 
            // to rank X, in the 2nd step the one belonging to ranks X and X-1, 
            // and so on ...
            size_t k = X;
            while (--k > 0) {
                d = max(d, m[k][j]);
                if(d >= t[k]) {
                    cc[k]++;
                }
            }
        }
        // Enforce monotonicity using succesive maximization
        for(size_t i=1; i<cc.size(); ++i) {
            cc[i] = max(cc[i-1], cc[i]);
        }
        return cc;
    }

    // Computes step down p-values
    std::deque<double> step_down_pvalues(
            const std::deque<double>& t,
            const matrix_type& m,
            const std::deque<double>& v)
    {
        using namespace std;
        deque<double> p(t.size());
        deque<size_t> cc = step_down_counts(t, m, v);
        size_t nperm = m.front().size(); //number of permutations
        for (size_t i=0; i<p.size(); ++i) {
            p[i] = double(cc[i] + 1)/(nperm + 1); 
        }
        return p;
    }
} // namespace statistic
} // namespace Permory

#endif // include guard
