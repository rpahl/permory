// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_permutation_perm_matrix_hpp
#define permory_permutation_perm_matrix_hpp

#include <valarray>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "detail/config.hpp"
#include "permutation/perm.hpp"

namespace Permory { namespace permutation {

    // Stores permutations and provides fast methods for counting
    // allele/genotype frequencies in these permutations
    template<class T> class Perm_matrix {
        public:
            typedef boost::dynamic_bitset<> bitset_t;

            // Ctor
            Perm_matrix(
                    const size_t nperm,
                    const Permutation& p, 
                    const std::vector<T>& trait,
                    bool useBitmat=true); 
            // Inspection
            size_t nperm() const { return tpermMat_.ncol(); }
            size_t nsubject() const { return tpermMat_.nrow(); }
            bool hasBitmat() const { return hasBitmat_; }

            // Modification
            void reshuffle(const size_t, const Permutation&);

            // For definitions of these functions see booster.hpp
            template<class T2> friend void rem(const Perm_matrix<T2>&, 
                    const bitset_t& b, const bitset_t&, std::valarray<T2>&); 
            template<class T2> friend void bar(const Perm_matrix<T2>&,
                    const bitset_t&, std::valarray<T2>&);
            template<class T2> friend void git(const Perm_matrix<T2>&,
                const std::vector<int>&, std::valarray<T2>&); 

        private:
            detail::Matrix<T> tpermMat_;    //transposed permutations
            std::vector<bitset_t> bitMat_;  //bit-coded permutations
            bool hasBitmat_;
    };

    template<class T> inline Perm_matrix<T>::Perm_matrix( 
            const size_t nperm,
            const Permutation& p, 
            const std::vector<T>& trait, 
            bool useBitmat
            ) 
        : tpermMat_(trait.size(), nperm), hasBitmat_(useBitmat)
    {
        if (useBitmat) {
            bitset_t bs(trait.size());
            bitMat_.resize(nperm, bs);
            //bitMat_.resize(nperm, bitset_t(trait.size()));
        }
        std::vector<T> v = trait;                
        for (size_t i=0; i<nperm; ++i) {  
            p.shuffle(&v[0], v.size()); //next permutation
            for (size_t j=0; j<v.size(); ++j) {
                tpermMat_[j][i] = v[j];     //fill by column 
            }
            if (useBitmat) {
                bitMat_[i] = detail::vector_to_bitset(v);
            }
        }
    }
    template<class T> inline void Perm_matrix<T>::reshuffle(
            const size_t nperm,
            const Permutation& p)
    {
        assert(tpermMat_.ncol() > 0);
        // Copy trait out of the first column
        std::vector<T> some_trait(tpermMat_.nrow());
        for (size_t i=0; i<tpermMat_.nrow(); i++) {
            some_trait[i] = tpermMat_[i][0];
        }
        *this = Perm_matrix(nperm, p, some_trait, hasBitmat_);
    }

} // namespace permutation
} // namespace Permory

#endif // include guard
