// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_locus_hpp
#define permory_locus_hpp

#include <string>
#include <vector>

#include "config.hpp"
#include "individual.hpp"

namespace Permory 
{
    class Individual; // forward declaration

    class Locus {
        public: 
            enum Chr {
                none=0, chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, 
                chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, 
                chr17, chr18, chr19, chr20, chr21, chr22, X, Y, XY, MT};
            std::vector<Permory::Individual>* individuals_; 

            // Ctor
            Locus(
                    size_t id,
                    std::string name="", 
                    Chr chr=none,
                    size_t bp = 0,
                    double cm = 0.0
                 ) 
                : id_(id), name_(name), chr_(chr), bp_(bp), cm_(cm),
                individuals_(0)
            { }

            // Inspection 
            size_t id() const { return id_; }
            std::string name() const { return name_; }
            Chr chr() const { return chr_; }
            size_t bp() const { return bp_; }
            double cm() const { return cm_; }
            bool operator<(size_t i) const { return id_ < i; }
            bool operator<(double d) const { return tsMax_ < d; }
            bool operator<(const Locus& x) const { 
                return (chr_ < x.chr() || (chr_ == x.chr() && bp_ < x.bp()) );
            }
            bool operator==(const Locus& x) const { return ( id_ == x.id() ); }
            std::vector<double> test_stats() const { return ts_; }
            double tmax() const { return tsMax_; }

            // Modification
            template<int K, int L> void add_test_stat(
                    std::vector<double>::const_iterator start,
                    std::vector<double>::const_iterator end);
        private:
            size_t id_;     //unique id
            string name_;     //rs-id or other Locus identifier
            Chr chr_;            //chr1-chr22, X, Y or na if undefined
            size_t bp_;          //base pair position in bp units
            double cm_;          //cM map position

            std::vector<double> ts_;    //test statistics for the locus
            double tsMax_;              //max of ts_
    };

    // ========================================================================
    // Locus implementations
    template<int K, int L> inline void Locus::add_test_stat(
            std::vector<double>::const_iterator start,
            std::vector<double>::const_iterator end)
    {
        copy (start, end, back_inserter(ts_));
        tsMax_ = std::max(tsMax_, *std::max_element(ts_.begin(), ts_.end()));
    }


} // namespace Permory

#endif // include guard
