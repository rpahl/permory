// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_gwas_hpp
#define permory_gwas_hpp

#include <deque>
#include <vector>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/deque.hpp>

#include "detail/config.hpp"
#include "individual.hpp"
#include "locus.hpp"

namespace Permory { namespace gwas {
    //using namespace ::boost::multi_index;
    //namespace bmi = boost::multi_index;

    // Genome-wide association study
    class Gwas {
        public:
            //typedef typename boost::ptr_deque<data_type>::iterator iterator;
            typedef std::deque<Locus>::iterator iterator;
            typedef std::deque<Locus>::const_iterator const_iterator;
            typedef std::vector<Individual>::const_iterator const_inderator;

            // Iterator pass through
            iterator begin() { return loci_.begin(); }
            iterator end() { return loci_.end(); }
            const_iterator begin() const { return loci_.begin(); }
            const_iterator end() const { return loci_.end(); }
            const_inderator ind_begin() const { return ind_.begin(); }
            const_inderator ind_end() const { return ind_.end(); }

            // Ctor
            Gwas(const std::vector<Individual>& ind);   
            Gwas(const_inderator start, const_inderator end);   

            // Inspection
            size_t m() const { return loci_.size(); }
            size_t meff(size_t maxval) const { return std::min(maxval, meff_); }
            size_t ncase() const;
            size_t ncontrol() const { return (sample_size() - ncase()); }
            size_t sample_size() const { return ind_.size(); }
            bool has_unique_loci() const;

            // Modification
            std::deque<Locus>* pointer_to_loci() { return &loci_; }
            void add_loci(const_iterator start, const_iterator end);
            void resize_loci(size_t n) { if (n < this->m()) loci_.resize(n); }
            void set_meff(double, double alpha=0.05, bool Bonf=true);

        private:
            std::vector<Individual> ind_;   //recruited individuals
            std::deque<Locus> loci_;
            size_t meff_;                   //effective number of tests

            // serialization stuff
            friend class boost::serialization::access;
            template<class Archive>
            void serialize(Archive & ar, const unsigned int version)
            {
                ar & ind_;
                ar & loci_;
            }

    };

    // Gwas implementation
    // ========================================================================
    inline Gwas::Gwas(
            const std::vector<Individual>& ind)
       : ind_(ind)
    {
        if (ind_.size() < 2) {
            throw std::invalid_argument("Gwas must have at least 2 individuals.");
        }
    }

    inline Gwas::Gwas(
            const_inderator start, const_inderator end)
        : ind_(start, end)
    {
        if (ind_.size() < 2) {
            throw std::invalid_argument("Gwas must have at least 2 individuals.");
        }
    }

    //
    // Effective number of tests (or markers, or hypotheses)
    inline void Gwas::set_meff(
            double p_adjusted,  //adjusted p-value
            double alpha,       //significance threshold
            bool Bonf)          //Bonferroni (true) or Sidak (false) correction
    {
        if (Bonf) { //Bonferroni
            meff_ = alpha/p_adjusted;
        }
        else {      // Sidak
            meff_ = log(1-alpha)/log(1-p_adjusted);
        }
    }

    inline size_t Gwas::ncase() const
    {
        return std::count_if(ind_.begin(), ind_.end(), 
                    std::mem_fun_ref(&Individual::isAffected));
    }

    inline bool Gwas::has_unique_loci() const
    {
        std::deque<Locus> loc(loci_);
        std::sort(loc.begin(), loc.end());
        return (std::adjacent_find(loc.begin(), loc.end()) == loc.end());
    }

    inline void Gwas::add_loci(
            const_iterator start, const_iterator end) 
    {
        std::copy(start, end, std::back_inserter(loci_));
    }
} // namespace gwas
} // namespace Permory

// Boost Serialization API Version Information
// ========================================================================
BOOST_CLASS_VERSION(Permory::gwas::Gwas, 1)

#endif // include guard

/*
   bmi::multi_index_container<
   Locus,
   indexed_by<
   ordered_unique<identity<Locus> >,
   ordered_non_unique<
   BOOST_MULTI_INDEX_CONST_MEM_FUN(Locus, double, tmax)
   >
   > 
   > loci_;
   */


