// Copyright (c) 2010 Roman Pahl
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
//#include "contab.hpp"
//#include "detail/functors.hpp" //pair_comp_2nd
#include "detail/parameter.hpp"
#include "gwas/locusdata.hpp"
//#include "gwas/gwas.hpp"
#include "permutation/booster.hpp"  //Bitset_with_count,
#include "permutation/perm.hpp"
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
            boost::ptr_vector<Perm_boost<T> > boosters_;

            // contains the intermediate result as contingency table in
            // dichotom and extension of the nominator and denominator in
            // quantitative.
            Matrix<T> extension_;

            std::vector<double> tMax_;  //max test statistics
    };
    // ========================================================================
    // Statistic implementations

    template<class T> template<class D> inline void
        Statistic<T>::do_permutation(const gwas::Locus_data<D>& data)
    {
        size_t maxcnt = 0;
        size_t worst_idx = 0;
        size_t card = data.domain_cardinality();
        std::vector<uint> boost_index(card);    //indices of applied booster
        std::vector<std::vector<int> > index_codes(card);
        std::vector<Bitset_with_count> dummy_codes(card);

        // the unique_iterator is defined in discretedata.hpp:
        // std::map<elem_type, count_type> unique_;//unique elements with counts
        typename gwas::Locus_data<D>::unique_iterator it = data.unique_begin();
        for (uint i=0; i < card; i++) {
            // Prepare the raw data in different formats (index code and dummy 
            // code), which are later used for boosting the permutation 
            D code_value = it->first;         
            size_t cnt = (size_t) it->second;   //#occurences of the code value
            index_codes[i] = index_code<D>(data.begin(), data.end(), code_value);
            dummy_codes[i] = dummy_code<D>(data.begin(), data.end(), code_value);

            // Determine index of most similar dummy code in the
            // booster's buffer. If the smallest distance between this
            // dummy code and the most similar one is smaller than 'cnt',
            // 'cnt' will contain this distance after the call.
            boost_index[i] = boosters_[i].find_most_similar(dummy_codes[i], cnt);

            // Keep track of the worst boostable element, which is the one
            // that shows both highest occurences of the code and the
            // highest distance to "neighboured" dummy codes
            if (cnt > maxcnt) {
                maxcnt = cnt;
                worst_idx = i;
            }
            it++;
        }

        // Permute for each allelic code, except the one that can be least
        // optimized/boosted for permutation. For this, the frequency is
        // simply derived via the marginal sum.
        extension_[worst_idx] = marginal_sum_; //init with marginal sum
        for (uint i=0; i < card; i++) {
            if (i != worst_idx) {
                boosters_[i].permute(
                        index_codes[i],     //index coded data
                        dummy_codes[i],     //dummy coded data
                        boost_index[i],     //index into booster's buffer
                        &extension_[i]);    //resulting sums
                extension_[worst_idx] -= extension_[i];
            }
        }
        // Since it was left out, the dummy code and the resulting case 
        // frequencies of the code belonging to the "worst index" must be added 
        // post hoc "by hand" to the buffer of the booster.
        boosters_[worst_idx].add_to_buffer(dummy_codes[worst_idx]);
        boosters_[worst_idx].add_to_buffer(extension_[worst_idx]);
    }


} // namespace statistic
} // namespace Permory

#endif // include guard

