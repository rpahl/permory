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
    template<class T, class E, uint L> class Statistic {
        public:
            // Iterator pass through
            typedef typename std::vector<double>::const_iterator const_iterator;
            const_iterator tmax_begin() const { return tMax_.begin(); }
            const_iterator tmax_end() const { return tMax_.end(); }

            // Inspection
            size_t size() const { return testPool_.size(); }

        protected:
            // This function does the "permutation work"
            template<class D> void do_permutation(const gwas::Locus_data<D>&);

            Test_pool<T> testPool_;
            E marginal_sum_;       // sum of all nomdenom_buf_ elements
            boost::ptr_vector<Perm_boost<E> > boosters_;

            // contains the intermediate result as contingency table in
            // dichotom and extension of the nominator and denominator in
            // quantitative.
            Matrix<E> extension_;

            std::vector<double> tMax_;  //max test statistics

            // For caching purpose
            std::vector<T> tabs_;
            std::vector<int> index_[L+1];   //L+1 integer vectors
            Bitset_with_count dummy_[L+1];  //L+1 bitsets with cached bit counts
    };
    // ========================================================================
    // Statistic implementations

    template<class T, class E, uint L> template<class D> inline void
        Statistic<T, E, L>::do_permutation(const gwas::Locus_data<D>& data)
    {
        if (not (data.domain_cardinality() == L+1)) {
            throw std::runtime_error("Bad domain cardinality in permutation test.");
        }
        uint boost_index[L+1];  // see below
        size_t worst_idx = 0;
        size_t maxcnt = 0;

        typename gwas::Locus_data<D>::unique_iterator
            it = data.unique_begin();
        for (uint i=0; i < L+1; i++) {
            size_t cnt = (size_t) it->second; //#occurences of the code
            D allelic_code = it->first;
            index_[i].clear();
            index_[i].reserve(cnt);
            index_code<D>(index_[i], data.begin(), data.end(), allelic_code);
            dummy_[i] = dummy_code<D>(data.begin(), data.end(), allelic_code);

            // Determine index of most similar dummy code in the
            // booster's buffer. If the smallest distance between this
            // dummy code and the most similar one is smaller than 'cnt',
            // 'cnt' will contain this distance after the call.
            boost_index[i] = boosters_[i].find_most_similar(dummy_[i], cnt);

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
        for (uint i=0; i < L+1; i++) {
            if (i != worst_idx) {
                boosters_[i].permute(
                        index_[i],          //index coded data
                        dummy_[i],          //dummy coded data
                        boost_index[i],     //index into booster's buffer
                        &extension_[i]);    //resulting sums
                extension_[worst_idx] -= extension_[i];
            }
        }
        // Since it was left out, the dummy code and the resulting case
        // frequencies of the skipped code are added post hoc "by hand"
        // to the buffer of the booster.
        boosters_[worst_idx].add_to_buffer(dummy_[worst_idx]);
        boosters_[worst_idx].add_to_buffer(extension_[worst_idx]);
    }


} // namespace statistic
} // namespace Permory

#endif // include guard

