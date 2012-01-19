// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_dichotom_hpp
#define permory_dichotom_hpp

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
#include "gwas/locusdata.hpp"
#include "permutation/booster.hpp"  //Bitset_with_count,
#include "permutation/perm.hpp"
#include "statistical/testpool.hpp"
#include "statistical/statistic.hpp"

namespace Permory { namespace statistic {
    using namespace permutation;

    //
    // Analyze genotype data with binary/dichotomous trait
    //
    template<uint K, uint L, class T = unsigned short int> class Dichotom
            : public Statistic<Con_tab<K,L>, T, L> {
        public:
            typedef Statistic<Con_tab<K,L>, T, L> S;

            Dichotom(
                    const detail::Parameter&,
                    gwas::Gwas::const_inderator ind_begin,//individuals to create
                    gwas::Gwas::const_inderator ind_end,  // dichotomous trait from
                    const Permutation* pp=0);       // pre-stored permutations

            // Modification
            void renew_permutations(
                    const Permutation* pp,  //creates the permutations
                    size_t nperm,           //number of permutations
                    size_t tail_size);      //parameter of the permutation booster

            // Conversion
            // Compute test statistics for the data
            template<class D> std::vector<double> test(const gwas::Locus_data<D>&);
            // Compute permutation test statistics
            template<class D> void permutation_test(const gwas::Locus_data<D>&);

        private:
            // Corrects wrong sums in contingency table which appear if
            // not all genotype values are present in marker.
            template<class D> void correct_tab(const gwas::Locus_data<D>&,
                    Con_tab<K, L>*&);

            std::vector<T> prepare_trait(gwas::Gwas::const_inderator begin,
                    gwas::Gwas::const_inderator end);

            std::vector<T> trait_;

            // For caching purpose
            bool useBitarithmetic_;
    };
    // ========================================================================
    // Dichotom implementations
    template<uint K, uint L, class T> inline
        Dichotom<K, L, T>::Dichotom(
                const detail::Parameter& par,
                gwas::Gwas::const_inderator ind_begin,
                gwas::Gwas::const_inderator ind_end,
                const Permutation* pp)
        : trait_(prepare_trait(ind_begin, ind_end)),
          useBitarithmetic_(par.useBar)
        {
            this->testPool_.add(par);
            this->marginal_sum_ = std::accumulate(trait_.begin(), trait_.end(), 0);
            bool yesPermutation = (pp != 0);
            if (yesPermutation) {
                renew_permutations(pp, par.nperm_block, par.tail_size);
            }
        }

    template<uint K, uint L, class T> inline
        void Dichotom<K, L, T>::renew_permutations(const Permutation* pp, size_t nperm,
                size_t tail_size)
    {
        this->tabs_.resize(nperm);
        this->tMax_.clear();
        this->tMax_.resize(nperm);
        this->extension_.resize(L+1, nperm); //one extra row to account for missings

        // Create and store permutations in matrix
        boost::shared_ptr<Perm_matrix<T> > pmat(
                new Perm_matrix<T>(nperm, *pp, trait_, useBitarithmetic_));

        // Prepare permutation booster
        this->boosters_.clear();
        this->boosters_.reserve(L+1);
        for (uint i=0; i<L+1; i++) {
            this->boosters_.push_back(new Perm_boost<T>(pmat, tail_size));
        }
    }

    template<uint K, uint L, class T> template<class D> inline
        std::vector<double> Dichotom<K, L, T>::test(const gwas::Locus_data<D>& data)
        {
            // Create contingency tab and analyze it with each test of the test pool
            Con_tab<K, L>* tab;
            if (data.hasMissings()) { //requires extra work to get rid of the missings
                std::vector<D> dd; dd.reserve(data.size());
                std::vector<T> tt; tt.reserve(data.size());

                typename gwas::Locus_data<D>::const_iterator it = data.begin();
                BOOST_FOREACH(T t, trait_) {
                    if (*it != data.get_undef()) { //only keep the valid values
                        dd.push_back(*it);
                        tt.push_back(t);
                    }
                    it++;
                }
                tab = make_Con_tab<T,D,K,L>(tt.begin(), tt.end(), dd.begin(), dd.end());
            }
            else {
                tab = make_Con_tab<T,D,K,L>(
                        trait_.begin(), trait_.end(), data.begin(), data.end());
            }
            bool correction_needed =
                    data.domain_cardinality() !=
                    (data.data_cardinality() + (data.hasMissings() ? 0 : 1));
            if (correction_needed) {
                correct_tab<D>(data, tab);
            }

            std::vector<double> v(this->testPool_.size());
            for_each_test(*tab, this->testPool_.begin(), this->testPool_.end(), v.begin());
            delete tab;
            return v;
        }

    template<uint K, uint L, class T> template<class D> inline void
        Dichotom<K, L, T>::correct_tab(const gwas::Locus_data<D>& data,
                Con_tab<K, L>*& tab)
        {
            using gwas::Locus_data;

            Con_tab<K, L> *result = new Con_tab<K, L>();
            size_t col = 0;
            size_t col_result = 0;
            typename Locus_data<D>::unique_iterator it = data.unique_begin();
            do {
                while (col_result < L && it->second == 0) {
                    ++col_result;
                    ++it;
                }
                if (col_result < L) {
                    for (uint row = 0; row < K; ++row) {
                        result->assign(row, col_result, tab->at(row, col));
                    }
                    ++it;
                    ++col_result;
                    ++col;
                }
            } while (col_result < L);

            delete tab;
            tab = result;
        }

    template<uint K, uint L, class T> template<class D> inline void
        Dichotom<K, L, T>::permutation_test(const gwas::Locus_data<D>& data)
        {
            assert (trait_.size() == data.size());
            bool yesPermutation = (not this->boosters_.empty());
            if (yesPermutation) {
                this->do_permutation(data);
            }

            // Fill contingency tables
            uint j = 0; //row index of matrix with case frequency results
            uint c = 0; //column index of contingency table

            // the unique_iterator is defined in discretedata.hpp:
            // std::map<elem_type, count_type> unique_;//unique elements with counts
            typename gwas::Locus_data<D>::unique_iterator uniques = data.unique_begin();
            for (; uniques!=data.unique_end(); uniques++) {
                bool ok = not (uniques->first == data.get_undef());
                if (ok) {
                    uint n = uniques->second; //frequency of both (cases + controls)
                    for (uint t=0; t<this->tabs_.size(); ++t) {
                        this->tabs_[t][0][c] = this->extension_[j][t];     //cases r[j]
                        this->tabs_[t][1][c] = n - this->extension_[j][t]; //controls s[j]
                    }
                    c++;
                }
                j++;
            }
            // For each permutation i (i.e. for each obtained contingency
            // table) compute the max over all test statistics, say max(i), and
            // then update tMax_[i] = max(tMax_[i], max(i))
            for_each_test_and_tab(this->tabs_, this->testPool_, this->tMax_.begin());
        }

    template<uint K, uint L, class T> std::vector<T>
        Dichotom<K, L, T>::prepare_trait(
                    gwas::Gwas::const_inderator ind_begin,
                    gwas::Gwas::const_inderator ind_end)
    {
        std::vector<T> result(ind_end - ind_begin);
        transform(ind_begin, ind_end, result.begin(),
                    std::mem_fun_ref(&Individual::isAffected));
        return result;
    }

} // namespace statistic
} // namespace Permory

#endif // include guard

