// Copyright (c) 2010-2011 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_permutation_fast_count_hpp
#define permory_permutation_fast_count_hpp

#include <valarray>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/circular_buffer.hpp>

#include "detail/config.hpp"
#include "detail/matrix.hpp" 
#include "detail/vector.hpp" 
#include "detail/pair.hpp"
#include "gwas/locusdata.hpp"
#include "perm_matrix.hpp"
#include "boost_algorithms.hpp"
#include "recode.hpp"

namespace Permory { namespace permutation {
    typedef boost::dynamic_bitset<> bitset_t;

    //
    // Helper class: we use this class in the buffer of class Fast_count in 
    // order to save repeatedly counting bits during search in the buffer
    //
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
    typedef Bitset_with_count Bitset_wc;
    inline size_t hamming_dist(
            const Bitset_wc& b1, 
            const Bitset_wc& b2) {
        return (b1.get() ^ b2.get()).count();
    }

    // Permutation booster
    template<class T> class Fast_count {
        public:
            typedef std::pair<Bitset_with_count, std::valarray<T> > elem_t;
            typedef boost::circular_buffer<elem_t> buffer_t;
            typedef typename buffer_t::const_reverse_iterator buf_iterator;

            // Ctor + Dtor
            Fast_count(
                    boost::shared_ptr<Perm_matrix<T> >, 
                    size_t buffer_sz=100); 
            ~Fast_count(){}

            // Inspector
            size_t get_tradeOff() const { return tradeOff_; }
            bool empty_buffer() const { return buf_.empty(); }
            bool hasMemoized() const { return (!empty_buffer()) && (itMem_ != buf_.rend()); }

            // Modifier
            void add_to_buffer(const Bitset_wc& b, const std::valarray<T>& v) {
                buf_.push_back(std::make_pair(b, v)); 
            }

            // Conversion
            std::valarray<T> operator()(const Bitset_with_count&);
            size_t find_similar_bitset_in_buffer(const Bitset_with_count&, size_t toBeat);

        private:
            boost::shared_ptr<Perm_matrix<T> > permMatrix_;
            size_t tradeOff_;

            // reconstruction memoization requires buffering input with results 
            buffer_t buf_;
            buf_iterator itMem_;    //Memorize position in buffer
    };

    // =====================================================================
    // Fast_count implementation
    template<class T> inline Fast_count<T>::Fast_count(
            boost::shared_ptr<Perm_matrix<T> > pmat, 
            size_t buffer_sz) 
        : permMatrix_(pmat), buf_(buffer_sz)
    { 
        tradeOff_ = (size_t)(double(pmat->nsubject())/6.0); 
        itMem_ = buf_.rend();
    }

    //
    // Takes a bitset b and searches for a bitset x in the buffer that shows the
    // smallest hamming distance between b and all bitsets in the buffer. The
    // function basically returns min(b.count(), "minimal hamming distance")
    //
    template<class T> inline size_t Fast_count<T>::find_similar_bitset_in_buffer(
            const Bitset_with_count& b, size_t toBeat)
    {
        if (permMatrix_->hasBitmat()) { toBeat = std::min(toBeat, this->tradeOff_); } //TODO auslagern mittels get_tradeOff()
        itMem_ = buf_.rend();
        for (buf_iterator itBuf=buf_.rbegin(); itBuf!=buf_.rend(); itBuf++) {
            if (toBeat == 0) break; //cannot improve anymore 

            // bit count is always a lower bound for the hamming distance, thus 
            // providing a quick estimate whether it is worth to compute the 
            // hamming distance at all
            bool isWorth = (itBuf->first).count() < toBeat;

            if (isWorth) { 
                size_t hd = hamming_dist(b, itBuf->first);
                if (hd < toBeat) {  //more similar bitset found?
                    toBeat = hd;
                    itMem_ = itBuf; //remember this position
                }
            }
        }
        return toBeat;
    }

    template<class T> inline std::valarray<T> Fast_count<T>::operator()(
            const Bitset_with_count& dummy_code) 
                    //boost::shared_ptr<Perm_matrix<T> > permMat) //TODO
    {
        std::vector<int> indices = index_code(dummy_code.get()); 

        // First determine, which of the accelerating methods are available
        bool useREM = this->hasMemoized();          //reconstruction memoization
        bool noBAR = !(permMatrix_->hasBitmat());   //bit arithmetics available?
        bool useGIT = indices.size() < tradeOff_;   //genotype indexing

        std::valarray<T> res(T(0), permMatrix_->nperm());
        if (useREM) {
            // recall/memoize previous results and update them using rem method
            res = itMem_->second;
            rem(*permMatrix_, dummy_code.get(), (itMem_->first).get(), res);
        }
        else if(useGIT 
                || noBAR) { //if BAR method NOT available, we use GIT anyway 
            res = 0;        //init *all* valarray entries with 0
            git(*permMatrix_, indices, res); 
        }
        else {
            bar(*permMatrix_, dummy_code.get(), res); 
        }
        return res;
    }

} // namespace permutation
} // namespace Permory

#endif // include guard


