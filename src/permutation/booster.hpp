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
#include "gwas/locusdata.hpp"
#include "permutation/perm_matrix.hpp"
#include "recode.hpp"

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
            size_t current_hamming_distance(const Bitset_wc& b) { 
                return empty_buffer() ? b.count() : hamming_dist(b, buf_->first); 
            }

            // Modifier
            void add_to_buffer(const Bitset_wc& b, const std::valarray<T>& v) {
                buf_.push_back(std::make_pair(b, v)); 
            }
            std::valarray<T> permute(
                    const std::vector<int>&,    //index-coded data
                    const Bitset_wc&,           //dummy-coded data
                    int);                       //index into buffer
            template<class D> std::valarray<T> count_fast(
                    const D val, 
                    typename std::vector<D>::const_iterator start, 
                    typename std::vector<D>::const_iterator end);

            // Conversion
            int find_most_similar(const Bitset_wc&, size_t&) const; //searches in buffer
            void find_similar_bitset_in_buffer(const Bitset_with_count&, size_t);

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
    // Takes a bitset b and searches for a bitset x in the buffer such that 
    // the hamming distance between b and x is smaller than the value specified
    // by 'toBeat'
    template<class T> inline int Fast_count<T>::find_most_similar(
            const Bitset_wc& b, size_t& toBeat) const
    {
        if (permMatrix_->hasBitmat()) {
            toBeat = std::min(toBeat, this->tradeOff_);
        }
        bool hasFound = false;
        buf_iterator itmin = buf_.rend();
        for (buf_iterator itBuf=buf_.rbegin(); itBuf!=buf_.rend(); itBuf++) {
            if (toBeat == 0) break; //cannot improve anymore 

            // bit count is a lower bound for the hamming distance, thus 
            // providing a quick estimate whether it is worth to compute the 
            // hamming distance 
            size_t x = b.count();   
            size_t y = (itBuf->first).count();
            bool isWorth = (std::max (x, y) - std::min (x, y)) < toBeat;

            if (isWorth) { 
                size_t hd = hamming_dist(b, itBuf->first);
                if (hd < toBeat) { //more similar bitset found?
                    hasFound = true;
                    toBeat = hd;
                    itmin = itBuf;
                }
            }
        }
        // If we have found a very similar bitset in buffer, return its index
        return hasFound ? std::distance(itmin, buf_.rend()) - 1 : -1;
    }

    /*
    template<class T> inline void Fast_count<T>::find_similar_bitset_in_buffer(
            const Bitset_with_count& b, size_t toBeat)
    {
        //if (permMatrix_->hasBitmat()) { toBeat = std::min(toBeat, this->tradeOff_); } //TODO auslagern mittels get_tradeOff()
        itMem_ = buf_.rend();
        for (buf_iterator itBuf=buf_.rbegin(); itBuf!=buf_.rend(); itBuf++) {
            if (toBeat == 0) break; //cannot improve anymore 

            // bit count is a lower bound for the hamming distance, thus 
            // providing a quick estimate whether it is worth to compute the 
            // hamming distance 
            size_t x = b.count();   
            size_t y = (itBuf->first).count();
            bool isWorth = (std::max (x, y) - std::min (x, y)) < toBeat;

            if (isWorth) { 
                size_t hd = hamming_dist(b, itBuf->first);
                if (hd < toBeat) {  //more similar bitset found?
                    toBeat = hd;
                    itMem_ = itBuf; //remember this position
                }
            }
        }
    }
    */

    template<class T> inline std::valarray<T> Fast_count<T>::permute(
            const std::vector<int>& index_code, //index-coded data
            const Bitset_wc& dummy,               //dummy-coded data
            int imin)                           //index into result buffer
    {
        // First determine, which of the boosting methods are available
        bool useREM = imin >= 0 && hasMemoized();  //reconstruction memoization
        bool noBAR = !(permMatrix_->hasBitmat());   //bit arithmetics available?
        bool useGIT = index_code.size() < tradeOff_;//genotype indexing

        std::valarray<T> res(T(0), permMatrix_->nperm());
        if (useREM) {
            // recall/memoize previous results and update them using rem method
            res = buf_[imin].second; 
            rem(*permMatrix_, dummy.get(), buf_[imin].first.get(), res);
        }
        // if BAR method is NOT used, we will use GIT in any case
        else if(useGIT || noBAR) {    
            res = 0;        //init *all* valarray entries with 0
            git(*permMatrix_, index_code, res); 
        }
        else {
            bar(*permMatrix_, dummy.get(), res); 
        }

        // update buffers
        this->add_to_buffer(dummy, res);
        return res;
    }

    /*
    template<class T> template<class D> inline 
        std::valarray<T> Fast_count<T>::count_fast(
                    const D val, 
                    typename std::vector<D>::const_iterator start, 
                    typename std::vector<D>::const_iterator end)
                    //boost::shared_ptr<Perm_matrix<T> > permMat) //TODO
    {
        std::vector<int> indices = index_code(start, end, val);
        Bitset_with_count dummy = dummy_code(start, end, val);

        // First determine, which of the accelerating methods are available
        bool useREM = this->hasMemoized();          //reconstruction memoization
        bool noBAR = !(permMatrix_->hasBitmat());   //bit arithmetics available?
        bool useGIT = indices.size() < tradeOff_;   //genotype indexing

        std::valarray<T> res(T(0), permMatrix_->nperm());
        if (useREM) {
            // recall/memoize previous results and update them using rem method
            res = *itRes_;
            rem(*permMatrix_, dummy.get(), itDat_->get(), res);
        }
        else if(useGIT 
                || noBAR) { //if BAR method NOT available, we use GIT anyway 
            res = 0;        //init *all* valarray entries with 0
            git(*permMatrix_, indices, res); 
        }
        else {
            bar(*permMatrix_, dummy.get(), res); 
        }

        // update buffers
        return res;
    }
    */

} // namespace permutation
} // namespace Permory

#endif // include guard


