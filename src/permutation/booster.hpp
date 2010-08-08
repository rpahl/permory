// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_permutation_booster_hpp
#define permory_permutation_booster_hpp

#include <valarray>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/circular_buffer.hpp>

#include "config.hpp"
#include "detail/bitset.hpp" 
#include "detail/matrix.hpp" 
#include "detail/recode.hpp"
#include "detail/vector.hpp" 
#include "locusdata.hpp"
#include "permutation/perm_matrix.hpp"

namespace Permory { namespace permutation {
    // Helper class
    class Bitset_with_count { 
        public:
            boost::dynamic_bitset<> bs;
            size_t cnt;

            // Ctors
            Bitset_with_count() { cnt = 0; } 
            Bitset_with_count(const boost::dynamic_bitset<>& bs) 
                : bs(bs) { cnt = bs.count(); }

            // Inspectors
            size_t count() const { return cnt; }
            size_t size() const { return bs.size(); }
            bool operator<(const Bitset_with_count& b) const { 
                return cnt < b.cnt; 
            }
    };
    // Helper function
    typedef Bitset_with_count Bitset2;
    inline size_t hamming_dist(
            const Bitset2& b1, 
            const Bitset2& b2) {
        return (b1.bs ^ b2.bs).count();
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
            std::pair<int, size_t> find_most_similar( //searches in dataBuf_
                    const Bitset2&, size_t) const;
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
        tradeOff_ = (size_t)(double(pmat->nsubject())/7.0); 
    }

    template<class T> inline void Perm_boost<T>::add_to_buffer(
            const Bitset2& b) {
        dataBuf_.push_back(b);                  
    }

    template<class T> inline void Perm_boost<T>::add_to_buffer(
            const std::valarray<T>& v) {
        resultBuf_.push_back(v);            
    }

    template<class T> inline std::pair<int, size_t> Perm_boost<T>::
        find_most_similar(const Bitset2& b, size_t toBeat) const
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
            // if we have improved, toBeat was also updated to a smaller value
            return std::make_pair(imin, toBeat); 
        }

    template<class T> inline void Perm_boost<T>::permute(
            const std::vector<int>& dix,    //index-coded data
            const Bitset2& dummy,           //dummy-coded data
            int imin,                       //index into result buffer
            std::valarray<T>* result       
            )
    {
        // tradeoff threshold regarding bit arithmetics vs indexing methods
        bool noBAR = !(permMatrix_->hasBitmat()); //no bit arithmetics anyway?
        bool usePAM = imin >= 0;                  //partial memoization
        bool useGIT = dix.size() < tradeOff_;     //genotype indexing
        bool isBufferEmpty = resultBuf_.size() == 0;

        /* XXX
           if (imin > 0) 
           cout << hamming_dist(dummy, dataBuf_[imin]) << endl;
           else
           cout << tradeOff_ << endl;
           */
        //print_vec(dummy.bs); std::cout << std::endl; 
        if (usePAM || useGIT || noBAR) {
            if (usePAM && !isBufferEmpty) { 
                *result = resultBuf_[imin];   //memorize previous results
                permMatrix_->rem(dummy.bs, dataBuf_[imin].bs, *result);
            }
            else {
                *result = 0;
                permMatrix_->git(dix, *result); 
            }
        }
        else {
            permMatrix_->bar(dummy.bs, *result); 
        }
        // update buffers
        dataBuf_.push_back(dummy);                  
        resultBuf_.push_back(*result);            
    }

} // namespace permutation
} // namespace Permory

#endif // include guard


