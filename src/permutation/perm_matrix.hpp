// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_permutation_perm_matrix_hpp
#define permory_permutation_perm_matrix_hpp

#include <valarray>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "config.hpp"
#include "permutation/perm.hpp"

namespace Permory { namespace permutation {
    typedef boost::dynamic_bitset<> bitset_t;

    // Stores permutations and provides fast methods for counting
    // allele/genotype frequencies in these permutations
    template<class T> class Perm_matrix {
        public:
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

            // Conversion
            // bit arithmetics
            void bar(const bitset_t&, std::valarray<T>&);
            // genotype indexing
            void git(const std::vector<int>&, std::valarray<T>&);
            // reconstruction memoization
            void rem(
                    const bitset_t& b,   //dummy-coded data
                    const bitset_t& bb,  //a bitset similar to b
                    std::valarray<T>& res); //results are written into res
            void print(size_t i) { //XXX debugging
                for (size_t j=0; j<tpermMat_.nrow(); j++)
                    cout << tpermMat_[j][i];
                cout << endl;
                print_vec(bitMat_[i]);
            }
        private:
            Matrix<T> tpermMat_;             //transposed permutations
            std::vector<bitset_t> bitMat_;   //bit-coded permutations
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
                bitMat_[i] = vector_to_bitset(v);
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

    template<class T> inline void Perm_matrix<T>::bar(
            const bitset_t& b, 
            std::valarray<T>& res) 
    {
        assert (b.size() == bitMat_.front().size());
        assert (res.size() == bitMat_.size());
        for (size_t i=0; i<res.size(); i++) {
            res[i] = (b & bitMat_[i]).count();
        }
    }

    template<class T> inline void Perm_matrix<T>::git(
            const std::vector<int>& v, 
            std::valarray<T>& res) 
    {
        BOOST_FOREACH(int i, v) {
            assert (i < tpermMat_.size());
            res += tpermMat_[i]; 
        }
    }

    template<class T> inline void Perm_matrix<T>::rem(
            const bitset_t& b,   //dummy-coded data
            const bitset_t& bb,  //a bitset similar to b
            std::valarray<T>& res) 
    {
        assert (b.size() == bb.size());
        assert (b.size() == tpermMat_.nrow());
        bitset_t b1 = b ^ bb;   
        bitset_t b2 = b & b1;   //b[i]==1 vs bb[i]==0 :"positions" are added
        b1 &= bb;               //b[i]==0 vs bb[i]==1 :"..." are subtracted 

        // "Jump" to each different position and add/subtract the 
        // corresponding transposed permutations
        size_t pos = b2.find_first(); 
        while (pos < b2.size()) {   
            res += tpermMat_[pos]; //adding
            pos = b2.find_next(pos);
        }
        pos = b1.find_first();
        while (pos < b1.size()) {   
            res -= tpermMat_[pos]; //subtracting
            pos = b1.find_next(pos);
        }
    }
} // namespace permutation
} // namespace Permory

#endif // include guard
