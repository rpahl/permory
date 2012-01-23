// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
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
#include "permutation/permutation.hpp"
#include "statistical/testpool.hpp"
#include "statistical/statistic.hpp"

namespace Permory { namespace statistic {
    using namespace permutation;
    using namespace Permory::detail;

    //
    // Analyze genotype data with continuous/quantitative trait.
    //
    template<uint L> class Quantitative
        : public Statistic<Pair<double> > {
            public:
                typedef detail::Pair<double> pair_t;

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
                    make_table(typename std::vector<pair_t>::const_iterator,
                            typename std::vector<pair_t>::const_iterator,
                            typename std::vector<D>::const_iterator,
                            typename std::vector<D>::const_iterator,
                            std::vector<pair_t>& tab);
                // Compute permutation test statistics
                template<class D> void permutation_test(const gwas::Locus_data<D>&);

            private:
                Test_pool<std::vector<pair_t> > testPool_; 
                std::vector<std::vector<pair_t> > pairs_;   

                std::vector<double> prepare_trait(gwas::Gwas::const_inderator begin,
                        gwas::Gwas::const_inderator end);

                // Corrects wrong sum tables which appear if genotypes are missing
                // in markers.
                template<class D> void correct_tab(const Locus_data<D>&,
                        std::vector<pair_t>&);

                std::vector<double> trait_;
                Trend_continuous* test_stat_;
                const std::vector<pair_t>& nomdenom_buf_;
                // nominator-denominator buffer:
                //   first  = \sum_{i=1}^{N} (Y_i - \mu_y)
                //   second = \sum_{i=1}^{N} (Y_i - \mu_y)^2
        };
    // ========================================================================
    // Quantitative implementations
    template<uint L> inline Quantitative<L>::Quantitative(
                const Parameter& par,
                gwas::Gwas::const_inderator ind_begin,
                gwas::Gwas::const_inderator ind_end,
                const Permutation* pp)
        : trait_(prepare_trait(ind_begin, ind_end)),
        test_stat_(new Trend_continuous(trait_)),
        nomdenom_buf_(test_stat_->get_buffer())
    {
        this->marginal_sum_ = test_stat_->get_sum();

        this->testPool_.add(test_stat_); // NOTE: Pointer will be deleted by
        // Test_pool.

        bool yesPermutation = (pp != 0);
        if (yesPermutation) {
            renew_permutations(pp, par.nperm_block, par.tail_size);
        }
    }

    template<uint L> inline
        void Quantitative<L>::renew_permutations(const Permutation* pp, size_t nperm,
                size_t tail_size)
        {
            this->pairs_.resize(nperm);
            for (size_t i = 0; i < this->pairs_.size(); ++i) {
                this->pairs_[i].resize(L+1);
            }

            this->tMax_.clear();
            this->tMax_.resize(nperm);

            boost::shared_ptr<Perm_matrix<pair_t> > pmat(
                    new Perm_matrix<pair_t>(nperm, *pp, nomdenom_buf_, false));

            // Prepare permutation booster
            this->boosters_.clear();
            this->boosters_.reserve(L+1);
            for (uint i=0; i<L+1; i++) {
                this->boosters_.push_back(new Fast_count<pair_t>(pmat, tail_size));
            }
        }

    template<uint L> template<class D> inline void
        Quantitative<L>::make_table(
                typename std::vector<pair_t>::const_iterator it_buf,
                typename std::vector<pair_t>::const_iterator buf_end,
                typename std::vector<D>::const_iterator it_data,
                typename std::vector<D>::const_iterator data_end,
                std::vector<pair_t>& result)
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

    template<uint L> template<class D> inline
        std::vector<double> Quantitative<L>::test(const gwas::Locus_data<D>& data)
        {
            std::vector<pair_t> tab;
            if (data.hasMissings()) {
                std::vector<D> dd; dd.reserve(data.size());
                std::vector<pair_t> tt; tt.reserve(data.size());

                typename gwas::Locus_data<D>::const_iterator it = data.begin();
                BOOST_FOREACH(pair_t t, nomdenom_buf_) {
                    if (*it != data.get_undef()) { //only keep the valid values
                        dd.push_back(*it);
                        tt.push_back(t);
                    }
                    it++;
                }
                make_table<D>(tt.begin(), tt.end(), dd.begin(), dd.end(), tab);
            }
            else {
                make_table<D>(nomdenom_buf_.begin(), nomdenom_buf_.end(),
                        data.begin(), data.end(), tab);
            }
            bool correction_needed =
                data.domain_cardinality() !=
                (data.data_cardinality() + (data.hasMissings() ? 0 : 1));
            if (correction_needed) {
                correct_tab<D>(data, tab);
            }
            test_stat_->update(data);
            std::vector<double> v(this->testPool_.size());
            for_each_test(tab, this->testPool_.begin(), this->testPool_.end(), v.begin());
            return v;
        }

    template<uint L> template<class D> inline void
        Quantitative<L>::permutation_test(const gwas::Locus_data<D>& data)
        {
            if (not (data.domain_cardinality() == L+1)) { 
                throw std::runtime_error("Bad domain cardinality in permutation test.");
            }
            test_stat_->update(data);
            this->extension_.resize(L+1, this->tMax_.size());  // one extra row for missings

            bool useBooster = (not this->boosters_.empty());
            if (useBooster) {
                this->do_permutation(data);
            }

            // Fill tables
            uint j = 0; //row index of matrix with case frequency results
            uint c = 0; //column index of contingency table

            // the unique_iterator is defined in discretedata.hpp:
            // std::map<elem_type, count_type> unique_;//unique elements with counts
            typename gwas::Locus_data<D>::unique_iterator uniques = data.unique_begin();
            for (; uniques!=data.unique_end(); uniques++) {
                bool ok = not (uniques->first == data.get_undef());
                if (ok) {
                    for (uint t=0; t < this->pairs_.size(); ++t) {
                        this->pairs_[t][c] = this->extension_[j][t];
                    }
                    c++;
                }
                j++;
            }
            // For each permutation i (i.e. for each obtained contingency
            // table) compute the max over all test statistics, say max(i), and
            // then update tMax_[i] = max(tMax_[i], max(i))
            each_test_for_each_element(this->pairs_, this->testPool_, this->tMax_.begin());
        }

    template<uint L> std::vector<double>
        Quantitative<L>::prepare_trait(
                gwas::Gwas::const_inderator ind_begin,
                gwas::Gwas::const_inderator ind_end)
        {
            std::vector<double> result;
            result.reserve(ind_end - ind_begin);
            for (gwas::Gwas::const_inderator i = ind_begin;
                    i != ind_end;
                    ++i) {
                result.push_back(i->begin()->val);
            }
            return result;
        }
    template<uint L> template<class D> inline void
        Quantitative<L>::correct_tab(const Locus_data<D>& data,
                std::vector<pair_t>& tab)
        {
            using std::vector;

            vector<pair_t> result;
            typename vector<pair_t>::const_iterator it_tab = tab.begin();
            for (typename Locus_data<D>::unique_iterator it = data.unique_begin();
                    it != data.unique_end();
                    ++it) {
                if (it->second == 0) {
                    result.push_back(make_pair<double>(0, 0));
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

