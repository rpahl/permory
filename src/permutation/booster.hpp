// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_permutation_booster_hpp
#define permory_permutation_booster_hpp

#include <valarray>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/circular_buffer.hpp>

#include "detail/config.hpp"
#include "detail/matrix.hpp" 
#include "detail/vector.hpp" 
#include "detail/pair.hpp"
#include "permutation/perm_matrix.hpp"
#include "recode.hpp"

namespace Permory { namespace permutation {
    typedef boost::dynamic_bitset<> bitset_t;

    // The following three functions (bar, git, and rem) are the core boosting
    // elements of the overall permutation boosting process. For each data
    // vector, which basically is the data corresponding to some marker, and for
    // each genotype within this data vector the fastest of the three methods 
    // is chosen at runtime. For more details, we refer to the publication:
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

    // Helper class
    class Bitset_with_count { 
        public:
            // Ctors
            Bitset_with_count() : cnt_(0) {}
            Bitset_with_count(const boost::dynamic_bitset<>& bs) 
                : bs_(bs), cnt_(bs.count())
            {}

            // Inspection
            const boost::dynamic_bitset<>& get() const { return bs_; }
            size_t count() const { return cnt_; }
            size_t size() const { return bs_.size(); }
            bool operator<(const Bitset_with_count& b) const { 
                return cnt_ < b.count(); 
            }

            // Modification
            Bitset_with_count& operator=(const boost::dynamic_bitset<>& bs) {
                bs_ = bs;
                cnt_ = bs.count();
                return *this;
            }
        private:
            boost::dynamic_bitset<> bs_;
            size_t cnt_;
    };

    // Helper function
    typedef Bitset_with_count Bitset2;
    inline size_t hamming_dist(
            const Bitset2& b1, 
            const Bitset2& b2) {
        return (b1.get() ^ b2.get()).count();
    }

    // Permutation booster
    template<class T> class Perm_boost {
        public:
            // Ctor + Dtor
            Perm_boost(
                    boost::shared_ptr<Perm_matrix<T> >, 
                    size_t buffer_sz=100); 
            ~Perm_boost(){}

            // Inspector
            size_t get_tradeOff() const { return tradeOff_; }

            // Modifier
            void add_to_buffer(const Bitset2&);
            void add_to_buffer(const std::valarray<T>&);
            void permute(
                    const std::vector<int>&,    //index-coded data
                    const Bitset2&,             //dummy-coded data
                    int,                        //index into buffers
                    std::valarray<T>* result);

            // Conversion
            int find_most_similar(const Bitset2&, size_t&) const; //searches in dataBuf_
        private:
            boost::shared_ptr<Perm_matrix<T> > permMatrix_;
            size_t tradeOff_;
            //partial memoization uses data and results buffer:
            boost::circular_buffer<Bitset2> dataBuf_;        
            boost::circular_buffer<std::valarray<T> > resultBuf_;
    };

    // =====================================================================
    // Perm_boost implementation
    template<class T> inline Perm_boost<T>::Perm_boost(
            boost::shared_ptr<Perm_matrix<T> > pmat, 
            size_t buffer_sz) 
        : permMatrix_(pmat), dataBuf_(buffer_sz), resultBuf_(buffer_sz)
    { 
        tradeOff_ = (size_t)(double(pmat->nsubject())/6.0); 
    }

    template<class T> inline void Perm_boost<T>::add_to_buffer(
            const Bitset2& b) {
        dataBuf_.push_back(b);                  
    }

    template<class T> inline void Perm_boost<T>::add_to_buffer(
            const std::valarray<T>& v) {
        resultBuf_.push_back(v);            
    }

    template<class T> inline int Perm_boost<T>::find_most_similar(
            const Bitset2& b, size_t& toBeat) const
    {
        typedef boost::circular_buffer<Bitset2> bitset_buffer;
        typedef bitset_buffer::const_reverse_iterator bufferator;
        if (permMatrix_->hasBitmat()) {
            toBeat = std::min(toBeat, this->tradeOff_);
        }
        bool haveImproved = false;
        bufferator itmin = dataBuf_.rbegin();
        for (bufferator i=dataBuf_.rbegin(); i!=dataBuf_.rend(); i++) {
            if (toBeat == 0) break; //cannot improve anymore 
            //bit count is a lower bound for the hamming distance, thus 
            //providing a first guess in constant time
            size_t x = b.count();   
            size_t y = (*i).count();
            bool isPromising = (std::max (x, y) - std::min (x, y)) < toBeat;
            if (isPromising) { 
                size_t hd = hamming_dist(b, *i);
                if (hd < toBeat) { //more similar bitset found?
                    haveImproved = true;
                    toBeat = hd;
                    itmin = i;
                }
            }
        }
        // return imin such that minimum can be acessed via dataBuf_[imin]
        int imin = -1;
        if (haveImproved) {
            imin = distance(itmin, dataBuf_.rend()) - 1;
        }
        return imin;
    }

    template<class T> inline void Perm_boost<T>::permute(
            const std::vector<int>& dix,    //index-coded data
            const Bitset2& dummy,           //dummy-coded data
            int imin,                       //index into result buffer
            std::valarray<T>* result)
    {
        // tradeoff threshold regarding bit arithmetics vs indexing methods
        bool noBAR = !(permMatrix_->hasBitmat()); //no bit arithmetics anyway?
        bool usePAM = imin >= 0;                  //partial memoization
        bool useGIT = !usePAM && (dix.size() < tradeOff_); //genotype indexing
        bool isBufferEmpty = resultBuf_.size() == 0;

        if (usePAM || useGIT || noBAR) {
            if (usePAM && !isBufferEmpty) { 
                *result = resultBuf_[imin];   //memorize previous results
                rem(*permMatrix_, dummy.get(), dataBuf_[imin].get(), *result);
            }
            else { //useGIT
                *result = 0;
                git(*permMatrix_, dix, *result); 
            }
        }
        else {
            bar(*permMatrix_, dummy.get(), *result); 
        }
        // update buffers
        dataBuf_.push_back(dummy);                  
        resultBuf_.push_back(*result);            
    }

} // namespace permutation
} // namespace Permory

#endif // include guard


