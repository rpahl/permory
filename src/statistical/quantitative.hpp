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

            // Conversion
            // Compute test statistics for the data
            template<class D> std::vector<double> test(const gwas::Locus_data<D>&);
            template<class D> std::vector<std::pair<T, T> >
                make_table(const gwas::Locus_data<D>&);
            // Compute permutation test statistics 
            template<class D> void permutation_test(const gwas::Locus_data<D>&);

        private:
            // This function does the "permutation work"
            template<class D> void do_permutation(const gwas::Locus_data<D>&);

            std::vector<T> trait_;
            double mu_y_;               // mean of phenotypes:
                                        //   1/N * \sum_{i=1}^{N} Y_i
            std::vector<std::pair<T,T> > nomdenom_buf_;
                                        // nominator-denominator buffer:
                                        //   first  = \sum_{i=1}^{N} (Y_i - \mu_y)
                                        //   second = \sum_{i=1}^{N} (Y_i - \mu_y)^2

            std::vector<double> tMax_; //max test statistics 

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
                T tmp = trait_[i] - mu_y_;
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
    }

    template<uint L, class T> template<class D> inline 
        std::vector<std::pair<T, T> > Quantitative<L, T>::make_table(const gwas::Locus_data<D>& data)
        {
            std::vector<std::pair<T, T> > result;
            result.resize(L);
            typename std::vector<std::pair<T, T> >::const_iterator
                it = nomdenom_buf_.begin();
            typename gwas::Locus_data<D>::const_iterator it_data = data.begin();
            while (it != nomdenom_buf_.end() && it_data != data.end()) {
                assert(*it_data < L);
                result[*it_data++] += *it++;
            }
            return result;
        }

    template<uint L, class T> template<class D> inline 
        std::vector<double> Quantitative<L, T>::test(const gwas::Locus_data<D>& data)
        {
            std::vector<std::pair<T, T> > tab;
            if (data.hasMissings()) {
            }
            else {
                tab = make_table(data);
            }
            std::vector<double> v;
            return v;
        }

    template<uint L, class T> template<class D> inline void 
        Quantitative<L, T>::permutation_test(const gwas::Locus_data<D>& data) 
        {
        }

    template<uint L, class T> template<class D> inline void 
        Quantitative<L, T>::do_permutation(const gwas::Locus_data<D>& data)
        {
        }
} // namespace statistic
} // namespace Permory

#endif // include guard

