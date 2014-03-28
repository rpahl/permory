// Copyright (c) 2010-2011 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_statistic_hpp
#define permory_statistic_hpp

//#include <algorithm>
//#include <numeric>
//#include <set>
//#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>

#include "detail/config.hpp"
#include "detail/parameter.hpp"
#include "gwas/locusdata.hpp"
#include "permutation/fast_count.hpp"  
#include "statistical/testpool.hpp"

namespace Permory { namespace statistic {
    using namespace permutation;
    using namespace Permory::detail;

    //
    // Base class for all classes analyzing genotype data.
    //
    template<class T> class Statistic {
        public:
            // Iterator pass through
            typedef typename std::vector<double>::const_iterator const_iterator;
            const_iterator tmax_begin() const { return tMax_.begin(); }
            const_iterator tmax_end() const { return tMax_.end(); }


        protected:
            // This function does the "permutation work"
            template<class D> void do_permutation(const gwas::Locus_data<D>&);

            T marginal_sum_;       // sum of all nomdenom_buf_ elements
            boost::ptr_vector<Fast_count<T> > boosters_;

            // contains the intermediate result as contingency table in
            // dichotom and extension of the nominator and denominator in
            // quantitative.
            Matrix<T> res_;   //intermediate results

            std::vector<double> tMax_;  //max test statistics
    };
    // ========================================================================
    // Statistic implementations

    template<class T> template<class D> inline void
        Statistic<T>::do_permutation(const gwas::Locus_data<D>& data)
    {
        size_t card = data.domain_cardinality();
        std::vector<uint> boost_index(card);    //indices of applied booster
        std::vector<Bitset_with_count> dummy_codes(card);

        // the unique_iterator is defined in discretedata.hpp:
        // std::map<elem_type, count_type> unique_;//unique elements with counts
        typename gwas::Locus_data<D>::unique_iterator it = data.unique_begin();
        size_t worst_dist = 0;
        size_t worst_idx = 0;
        for (uint i=0; i < card; i++) {
            D val = it->first;         
            dummy_codes[i] = dummy_code<D>(data.begin(), data.end(), val);

            // Search for the most similar bitset in the buffer with respect to 
            // hamming distance. The lower bound of the distance is the bit 
            // count of the dummy code.
            size_t cnt = dummy_codes[i].count();
            size_t dist = boosters_[i].find_similar_bitset_in_buffer(dummy_codes[i], cnt);

            // Keep track of the worst boostable element, which is the one with
            // highest occurences of the code and/or the least similarity to
            // codes in the buffer
            if (dist > worst_dist) {
                worst_dist = dist;
                worst_idx = i;
            }
            it++;
        }

        // Permute for each genotype code, except the one that can be least
        // boosted for permutation, and for which thus the result will be 
        // derived using the marginal sum.
        res_[worst_idx] = marginal_sum_; //init with marginal sum
        for (uint i=0; i < card; i++) {
            if (i != worst_idx) {
                res_[i] = boosters_[i](dummy_codes[i]);
                res_[worst_idx] -= res_[i];
                boosters_[i].add_to_buffer(dummy_codes[i], res_[i]);
            }
        }
        // Finally update booster's buffer of the "worst index" 
        boosters_[worst_idx].add_to_buffer(dummy_codes[worst_idx], res_[worst_idx]);
    }


} // namespace statistic
} // namespace Permory

#endif // include guard

