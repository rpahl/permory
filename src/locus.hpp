// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_locus_hpp
#define permory_locus_hpp

#include <sstream>
#include <string>
#include <vector>

#include "detail/config.hpp"

namespace Permory 
{
    // Modelling a genetic locus, for which test statistics are computed
    class Locus {
        public: 
            enum Chr {
                none=0, chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, 
                chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, 
                chr17, chr18, chr19, chr20, chr21, chr22, X, Y, XY, MT};

            // Ctor
            Locus(
                    size_t id,              //unique id
                    std::string rs="",      //rs-id or other Locus identifier
                    std::string gene="",    //gene the locus belongs to
                    Chr chr=none,           //chr1-chr22, X, Y or na if undefined
                    size_t bp = 0,          //base pair position in bp units
                    double cm = 0.0,        //cM map position
                    bool isPolymorph=true   //polymorph yes/no
                 ) 
                : id_(id), rs_(rs), gene_(gene), chr_(chr), bp_(bp), cm_(cm),
                isPolymorph_(isPolymorph), tsMax_(0.0), flag_(false)
            { }

            // Inspection 
            size_t id() const { return id_; }
            bool isPolymorph() const { return isPolymorph_; }
            bool isFlagged() const { return flag_; }
            bool hasTeststat() const { return (not ts_.empty()); }
            const std::string& rs() const { return rs_; }
            const std::string& gene() const { return gene_; }
            Chr chr() const { return chr_; }
            size_t bp() const { return bp_; }
            double cm() const { return cm_; }
            double tmax() const { return tsMax_; }
            bool operator<(const Locus& x) const;
            bool operator==(const Locus& x) const; 
            std::vector<double> test_stats() const { return ts_; }

            // Modification
            void set_flag(bool x) { flag_ = x; }
            void set_polymorph(bool x) { isPolymorph_ = x; }
            void set_gene(const std::string& s) { gene_ = s; }
            void add_test_stats(const std::vector<double>&);
        private:
            size_t id_;     //unique id
            std::string rs_;     //rs-id or other Locus identifier
            std::string gene_;   //gene the locus belongs to
            Chr chr_;       //chr1-chr22, X, Y or na if undefined
            size_t bp_;     //base pair position in bp units
            double cm_;     //cM map position

            std::vector<double> ts_;    //test statistics for the locus
            double tsMax_;              //max of ts_
            bool isPolymorph_;          //polymorph yes/no
            bool flag_;                 //flag this locus
    };

    // Locus implementation
    // ========================================================================
    inline bool Locus::operator<(const Locus& x) const
    {
        if (chr_ != none && bp_ > 0) {
            if (chr_ == x.chr())
                return bp_ < x.bp();
            else
                return chr_ < x.chr();
        }
        else
            return id_ < x.id();
    }
    inline bool Locus::operator==(const Locus& x) const
    {
        if (chr_ != none && bp_ > 0) {
            if (chr_ == x.chr())
                return bp_ == x.bp();
            else
                return chr_ == x.chr();
        }
        else
            return id_ == x.id();
    }
    inline void Locus::add_test_stats(const std::vector<double>& v)
    {
        copy (v.begin(), v.end(), back_inserter(ts_));
        tsMax_ = std::max(tsMax_, *std::max_element(ts_.begin(), ts_.end()));
    }


    Locus::Chr string2chr(const std::string& s)
    {
        if(s.empty())
            return Locus::none;
        std::istringstream iss(s);
        int i;
        iss >> i;
        if (i >= 1 && i <= 22) {
            return Locus::Chr(i);
        }
        if (s == "X")
            return Locus::X;
        else if (s == "Y")
            return Locus::Y;
        else if (s == "XY")
            return Locus::XY;
        else if (s == "MT")
            return Locus::MT;
        else
            return Locus::none;
    }
} // namespace Permory

#endif // include guard
