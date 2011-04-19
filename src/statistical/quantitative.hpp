// Copyright (c) 2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_quantitative_hpp
#define permory_quantitative_hpp

#include <algorithm>
#include <numeric>
#include <set>
#include <vector>
#include <utility> // pair

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "detail/config.hpp"
#include "contab.hpp"
#include "detail/functors.hpp" //pair_comp_2nd
#include "detail/parameter.hpp"
#include "detail/pair.hpp"
#include "gwas/locusdata.hpp"
#include "permutation/booster.hpp"  //Bitset_with_count,
#include "permutation/perm.hpp"
#include "statistical/testpool.hpp"

namespace Permory { namespace statistic {
    using namespace permutation;

    //
    // Analyze genotype data with binary/dichotomous trait
    //
    template<uint L, class T = double> class Quantitative {
        public:

            typedef std::pair<T, T> element_t;

            // Iterator pass through
            typedef typename std::vector<double>::const_iterator const_iterator;
            const_iterator tmax_begin() const { return tMax_.begin(); }
            const_iterator tmax_end() const { return tMax_.end(); }

            Quantitative(
                    const detail::Parameter&,
                    const std::vector<T>& trait,    //quantitative trait
                    const Permutation* pp=0);       //pre-stored permutations

            // Inspection

            // Modification
            void renew_permutations(
                    const Permutation* pp,  //creates the permutations
                    size_t nperm,           //number of permutations
                    size_t tail_size);      //parameter of the permutation booster
            template<class D> void calculate_mu_j(const gwas::Locus_data<D>& data);
            void calculate_denom_invariant();

            // Conversion
            // Compute test statistics for the data
            template<class D> std::vector<double> test(const gwas::Locus_data<D>&);
            template<class D> std::vector<element_t>
                make_table(const gwas::Locus_data<D>&);
            void test_for_tab(std::vector<element_t>&, std::vector<double>*);
            // Compute permutation test statistics
            template<class D> void permutation_test(const gwas::Locus_data<D>&);

        private:
            // This function does the "permutation work"
            template<class D> void do_permutation(const gwas::Locus_data<D>&);

            std::vector<T> trait_;
            double mu_y_;               // mean of phenotypes:
                                        //   1/N * \sum_{i=1}^{N} Y_i
            double mu_j_;               // mean of genotypes of marker j:
                                        //   1/N * \sum_{i=1}^{M} X_ji
            double denom_invariant_;    // invariant summand of denominator in
                                        // marker j:
                                        //   \mu_j^2 * \sum_{i=1}^{N} (Y_i - \mu_y)^2
            std::vector<element_t> nomdenom_buf_;
                                        // nominator-denominator buffer:
                                        //   first  = \sum_{i=1}^{N} (Y_i - \mu_y)
                                        //   second = \sum_{i=1}^{N} (Y_i - \mu_y)^2

            std::vector<double> tMax_;  //max test statistics

            boost::shared_ptr<Perm_matrix<element_t> > permMatrix_;

            detail::Matrix<element_t> extension_;
                                        // contains the intermediate
                                        // result (extension) of the nominator
                                        // and denominator

            // For caching purpose
            bool useBitarithmetic_;
    };
    // ========================================================================
    // Quantitative implementations
    template<uint L, class T> inline
        Quantitative<L, T>::Quantitative(
                const detail::Parameter& par,
                const std::vector<T>& trait,
                const Permutation* pp)
        : trait_(trait.begin(), trait.end()), useBitarithmetic_(par.useBar)
        {
            if (useBitarithmetic_) {
                throw std::invalid_argument(
                    "Bit arithmetics not allowed for quantitative phenotypes!");
            }

            mu_y_ = std::accumulate(trait_.begin(), trait_.end(), 0.)
                  / trait_.size();
            nomdenom_buf_.resize(trait_.size());
            for (size_t i = 0; i < trait_.size(); ++i) {
                const T tmp = trait_[i] - mu_y_;
                nomdenom_buf_[i] = make_pair(tmp, tmp*tmp);
            }

            bool yesPermutation = (pp != 0);
            if (yesPermutation) {
                renew_permutations(pp, par.nperm_block, par.tail_size);
            }
        }

    template<uint L, class T> inline
        void Quantitative<L, T>::renew_permutations(const Permutation* pp, size_t nperm,
                size_t tail_size)
    {
        tMax_.clear();
        tMax_.resize(nperm);

        boost::shared_ptr<Perm_matrix<element_t> > pmat(
                new Perm_matrix<element_t>(nperm, *pp, nomdenom_buf_, useBitarithmetic_));
        permMatrix_.swap(pmat);
    }

    template<uint L, class T>
    template<class D> void Quantitative<L, T>::calculate_mu_j(
                const gwas::Locus_data<D>& data)
    {
        mu_j_ = std::accumulate(data.begin(), data.end(), 0.) / data.size();
    }

    template<uint L, class T>
    void Quantitative<L, T>::calculate_denom_invariant()
    {
        typedef element_t P;
        denom_invariant_ = 0;
        BOOST_FOREACH(P x, nomdenom_buf_) {
            denom_invariant_ += x.second;
        }
        denom_invariant_ *= mu_j_ * mu_j_;
    }

    template<uint L, class T> template<class D> inline
        std::vector<typename Quantitative<L, T>::element_t>
        Quantitative<L, T>::make_table(const gwas::Locus_data<D>& data)
        {
            std::vector<element_t> result;
            result.resize(L);
            typename std::vector<element_t>::const_iterator
                it = nomdenom_buf_.begin();
            typename gwas::Locus_data<D>::const_iterator it_data = data.begin();
            while (it != nomdenom_buf_.end() && it_data != data.end()) {
                assert(*it_data < L);
                result[*it_data++] += *it++;
            }
            return result;
        }

    template<uint L, class T>
        void Quantitative<L, T>::test_for_tab(
                std::vector<element_t>& tab,
                std::vector<double> *result)
    {
        double nominator = 0.;
        double denominator = denom_invariant_;
        for (size_t i = 0; i < tab.size(); ++i) {
            element_t x = tab[i];
            nominator += i * x.first;
            denominator += (i * i - 2 * i * mu_j_) * x.second;
        }
        nominator *= nominator;
        result->push_back(nominator / denominator);
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
                tab = make_table(data);
            }
            calculate_mu_j(data);
            calculate_denom_invariant();
            std::vector<double> v;
            test_for_tab(tab, &v);
            return v;
        }

    template<uint L, class T> template<class D> inline void
        Quantitative<L, T>::permutation_test(const gwas::Locus_data<D>& data)
        {
            calculate_mu_j(data);
            calculate_denom_invariant();
            extension_.resize(L+1, tMax_.size());  // one extra row for missings

            do_permutation(data);

            std::vector<double> result;
            result.reserve(extension_.ncol());
            for (size_t i = 0; i < extension_.ncol(); ++i) {
                std::vector<element_t> tab(extension_.nrow());
                for (size_t j = 0; j < tab.size(); ++j) {
                    tab[j] = extension_[j][i];
                }
                test_for_tab(tab, &result);
            }
            for (size_t i = 0; i < result.size(); ++i) {
                tMax_[i] = std::max(tMax_[i], result[i]);
            }
        }

    template<uint L, class T> template<class D> inline void
        Quantitative<L, T>::do_permutation(const gwas::Locus_data<D>& data)
        {
            std::vector<int> index_[L+1];   //L+1 integer vectors

            typename gwas::Locus_data<D>::unique_iterator it = data.unique_begin();
            for (uint i=0; i < L+1; i++) {
                size_t cnt = (size_t) it->second; //#occurences of the code
                D allelic_code = it->first;
                index_[i].clear();
                index_[i].reserve(cnt);
                index_code<D>(index_[i], data.begin(), data.end(), allelic_code);
                ++it;

                git(*permMatrix_, index_[i], extension_[i]);
            }
        }
} // namespace statistic
} // namespace Permory

#endif // include guard

