// Copyright (c) 2010-2011 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_permutation_boost_algorithms_hpp
#define permory_permutation_boost_algorithms_hpp

#include <valarray>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "detail/config.hpp"
#include "perm_matrix.hpp"

namespace Permory { namespace permutation {
    typedef boost::dynamic_bitset<> bitset_t;

    // The following three functions (bar, git, and rem) are the core elements 
    // to accelerate the overall permutation process. For each data vector, 
    // which basically is the genotype data corresponding to some marker, and for 
    // each genotype (e.g. 0,1,2) of this vector, that of the three methods is 
    // chosen at runtime, which for each permutation can determine the
    // corresponding genotype frequencies the fastest.
    // For more details, we refer to the publication:
    // PERMORY: an LD-exploiting permutation test algorithm for powerful genome-wide 
    // association testing. Bioinformatics. 2010 Sep 1;26(17):2093-100. Epub 2010 Jul 6.
    
    //
    // *b*it *ar*ithmetics (BAR)
    //
    template<class T> inline void bar(   
            const Perm_matrix<T>& pmat, //matrix of predefined permutations
            const bitset_t& b,          //bitset aka dummy coded data
            std::valarray<T>& res)      //results are written into res
    {
        assert (b.size() == pmat.bitMat_.front().size());
        assert (res.size() == pmat.bitMat_.size());
        for (size_t i=0; i<res.size(); i++) {
            res[i] = (b & pmat.bitMat_[i]).count();
        }
    }

    //
    // *g*enotype *i*ndexing using *t*ransposed permutations (GIT)
    //
    template<class T> inline void git(   
            const Perm_matrix<T>& pmat, //matrix of predefined permutations
            const std::vector<int>& idx,//genotype index vector
            std::valarray<T>& res)      //results are written into res
    {
        BOOST_FOREACH(int i, idx) {
            assert (size_t(i) < pmat.tpermMat_.size());
            res += pmat.tpermMat_[i]; 
        }
    }

    //
    // *re*construction *m*emoization (REM)
    //
    template<class T> inline void rem(   
            const Perm_matrix<T>& pmat, //matrix of predefined permutations
            const bitset_t& b,          //dummy-coded data
            const bitset_t& bb,         //a bitset similar to b
            std::valarray<T>& res)      //results are written into res
    {
        assert (b.size() == bb.size());
        assert (b.size() == pmat.tpermMat_.nrow());
        bitset_t b1 = b ^ bb;   
        bitset_t b2 = b & b1; //b[i]==1 vs bb[i]==0 - "positions" will be added
        b1 &= bb;             //b[i]==0 vs bb[i]==1 - ... will be subtracted 

        // "Jump" to each differing position and add/subtract the 
        // corresponding transposed permutations
        size_t pos = b2.find_first(); 
        while (pos < b2.size()) {   
            res += pmat.tpermMat_[pos]; //adding
            pos = b2.find_next(pos);
        }
        pos = b1.find_first();
        while (pos < b1.size()) {   
            res -= pmat.tpermMat_[pos]; //subtracting
            pos = b1.find_next(pos);
        }
    }

} // namespace permutation
} // namespace Permory

#endif // include guard


