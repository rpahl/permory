// Copyright (c) 2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_quantitative_hpp
#define permory_quantitative_hpp

#include <algorithm>
#include <numeric>
#include <set>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "detail/config.hpp"
#include "contab.hpp"
#include "detail/functors.hpp" //pair_comp_2nd
#include "detail/parameter.hpp"
#include "detail/pair.hpp"
#include "gwas/locusdata.hpp"
#include "gwas/gwas.hpp"
#include "permutation/booster.hpp"  //Bitset_with_count,
#include "permutation/perm.hpp"
#include "statistical/testpool.hpp"

namespace Permory { namespace statistic {
    using namespace permutation;
    using namespace Permory::detail;

    //
    // Analyze genotype data with binary/dichotomous trait
    //
    template<uint L, class T = double> class Quantitative {
        public:

            typedef Pair<T> element_t;

            // Iterator pass through
            typedef typename std::vector<double>::const_iterator const_iterator;
            const_iterator tmax_begin() const { return tMax_.begin(); }
            const_iterator tmax_end() const { return tMax_.end(); }

            Quantitative(
                    const Parameter&,
                    gwas::Gwas::const_inderator ind_begin,//individuals to create
                    gwas::Gwas::const_inderator ind_end,  // dichotomous trait from
                    const Permutation* pp=0);       //pre-stored permutations

            // Inspection

            // Modification
            void renew_permutations(
                    const Permutation* pp,  //creates the permutations
                    size_t nperm,           //number of permutations
                    size_t tail_size);      //parameter of the permutation booster

            // Conversion
            // Compute test statistics for the data
            template<class D> std::vector<double> test(const gwas::Locus_data<D>&);
            template<class D> void
                make_table(typename std::vector<element_t>::const_iterator,
                           typename std::vector<element_t>::const_iterator,
                           typename std::vector<D>::const_iterator,
                           typename std::vector<D>::const_iterator,
                           std::vector<element_t>& tab);
            // Compute permutation test statistics
            template<class D> void permutation_test(const gwas::Locus_data<D>&);

        private:
            // This function does the "permutation work"
            template<class D> void do_permutation(const gwas::Locus_data<D>&);

            std::vector<T> prepare_trait(gwas::Gwas::const_inderator begin,
                    gwas::Gwas::const_inderator end);

            // Corrects wrong sum tables which appear if genotypes are missing
            // in markers.
            template<class D> void correct_tab(const Locus_data<D>&,
                    std::vector<element_t>&);

            std::vector<T> trait_;
            Trend_continuous<T> *test_stat_;
            Test_pool<std::vector<element_t> > testPool_;
            const std::vector<element_t>& nomdenom_buf_;
                                        // nominator-denominator buffer:
                                        //   first  = \sum_{i=1}^{N} (Y_i - \mu_y)
                                        //   second = \sum_{i=1}^{N} (Y_i - \mu_y)^2
            element_t sum_;             // sum of all nomdenom_buf_ elements

            std::vector<double> tMax_;  //max test statistics

            Matrix<element_t> extension_;
                                        // contains the intermediate
                                        // result (extension) of the nominator
                                        // and denominator

            boost::ptr_vector<Perm_boost<element_t> > boosters_;

            // For caching purpose
            std::vector<std::vector<element_t> > tabs_;
            bool useBitarithmetic_;
            std::vector<int> index_[L+1];   //L+1 integer vectors
            Bitset_with_count dummy_[L+1];  //L+1 bitsets with cached bit counts
    };
    // ========================================================================
    // Quantitative implementations
    template<uint L, class T> inline
        Quantitative<L, T>::Quantitative(
                const Parameter& par,
                gwas::Gwas::const_inderator ind_begin,
                gwas::Gwas::const_inderator ind_end,
                const Permutation* pp)
        : trait_(prepare_trait(ind_begin, ind_end)),
          test_stat_(new Trend_continuous<T>(trait_)),
          nomdenom_buf_(test_stat_->get_buffer()),
          useBitarithmetic_(par.useBar)
        {
            if (useBitarithmetic_) {
                throw std::invalid_argument(
                    "Bit arithmetics not allowed for quantitative phenotypes!");
            }

            testPool_.add(test_stat_); // NOTE: Pointer will be deleted by
                                       // Test_pool.
            sum_ = std::accumulate(nomdenom_buf_.begin(), nomdenom_buf_.end(),
                                    make_pair(0.,0.));

            bool yesPermutation = (pp != 0);
            if (yesPermutation) {
                renew_permutations(pp, par.nperm_block, par.tail_size);
            }
        }

    template<uint L, class T> inline
        void Quantitative<L, T>::renew_permutations(const Permutation* pp, size_t nperm,
                size_t tail_size)
    {
        tabs_.resize(nperm);
        for (size_t i = 0; i < tabs_.size(); ++i) {
            tabs_[i].resize(L+1);
        }

        tMax_.clear();
        tMax_.resize(nperm);

        boost::shared_ptr<Perm_matrix<element_t> > pmat(
                new Perm_matrix<element_t>(nperm, *pp, nomdenom_buf_, useBitarithmetic_));

        // Prepare permutation booster
        boosters_.clear();
        boosters_.reserve(L+1);
        for (uint i=0; i<L+1; i++) {
            boosters_.push_back(new Perm_boost<element_t>(pmat, tail_size));
        }
    }

    template<uint L, class T> template<class D> inline void
        Quantitative<L, T>::make_table(
            typename std::vector<element_t>::const_iterator it_buf,
            typename std::vector<element_t>::const_iterator buf_end,
            typename std::vector<D>::const_iterator it_data,
            typename std::vector<D>::const_iterator data_end,
            std::vector<element_t>& result)
    {
        using namespace std;

        result.resize(L);

        set<D> unique(it_data, data_end);
        vector<D> indexed(unique.begin(), unique.end());

        assert(indexed.size() <= L);
        while (it_buf != buf_end && it_data != data_end) {
            uint index = distance(indexed.begin(),
                            find(indexed.begin(), indexed.end(), *it_data++));
            assert(index < result.size());
            result[index] += *it_buf++;
        }
    }

    template<uint L, class T> template<class D> inline
        std::vector<double> Quantitative<L, T>::test(const gwas::Locus_data<D>& data)
        {
            std::vector<element_t> tab;
            if (data.hasMissings()) {
                // TODO
                throw std::runtime_error("Not implementet yet.");
            }
            else {
                make_table<D>(nomdenom_buf_.begin(), nomdenom_buf_.end(),
                              data.begin(), data.end(), tab);
            }
            bool correction_needed =
                    data.domain_cardinality() != data.data_cardinality();
            if (correction_needed) {
                correct_tab<D>(data, tab);
            }
            test_stat_->update(data);
            std::vector<double> v(testPool_.size());
            for_each_test(tab, testPool_.begin(), testPool_.end(), v.begin());
            return v;
        }

    template<uint L, class T> template<class D> inline void
        Quantitative<L, T>::permutation_test(const gwas::Locus_data<D>& data)
        {
            test_stat_->update(data);
            extension_.resize(L+1, tMax_.size());  // one extra row for missings

            do_permutation(data);

            // Fill tables
            uint j = 0; //row index of matrix with case frequency results
            uint c = 0; //column index of contingency table

            // the unique_iterator is defined in discretedata.hpp:
            // std::map<elem_type, count_type> unique_;//unique elements with counts
            typename gwas::Locus_data<D>::unique_iterator uniques = data.unique_begin();
            for (; uniques!=data.unique_end(); uniques++) {
                bool ok = not (uniques->first == data.get_undef());
                if (ok) {
                    for (uint t=0; t < tabs_.size(); ++t) {
                        tabs_[t][c] = extension_[j][t];
                    }
                    c++;
                }
                j++;
            }
            for_each_test_and_tab(tabs_, testPool_, tMax_.begin());
        }

    template<uint L, class T> template<class D> inline void
        Quantitative<L, T>::do_permutation(const gwas::Locus_data<D>& data)
        {
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
            extension_[worst_idx] = sum_; //init with marginal sum
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

    template<uint L, class T> std::vector<T>
        Quantitative<L, T>::prepare_trait(
                    gwas::Gwas::const_inderator ind_begin,
                    gwas::Gwas::const_inderator ind_end)
    {
        std::vector<T> result;
        result.reserve(ind_end - ind_begin);
        for (gwas::Gwas::const_inderator i = ind_begin;
                i != ind_end;
                ++i) {
            result.push_back(i->begin()->val);
        }
        return result;
    }
    template<uint L, class T> template<class D> inline void
        Quantitative<L, T>::correct_tab(const Locus_data<D>& data,
                    std::vector<element_t>& tab)
    {
        using std::vector;

        vector<element_t> result;
        typename vector<element_t>::const_iterator it_tab = tab.begin();
        for (typename Locus_data<D>::unique_iterator it = data.unique_begin();
             it != data.unique_end();
             ++it) {
            if (it->second == 0) {
                result.push_back(make_pair<T>(0, 0));
            }
            else {
                result.push_back(*it_tab++);
            }
        }
        tab.swap(result);
    }


} // namespace statistic
} // namespace Permory

#endif // include guard

