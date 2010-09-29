// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_snp_data_hpp
#define permory_snp_data_hpp

#include <numeric>

#include<boost/bind.hpp>
#include<boost/lexical_cast.hpp>

#include "detail/config.hpp"
#include "detail/enums.hpp" //Genetic_type
#include "discretedata.hpp"
#include "locus.hpp"

namespace Permory {
    using namespace detail;

    template<class T> class Locus_data : public Discrete_data<T> {
        public:
            typedef T elem_t; //internal type of the data (e.g. int, char)
            typedef typename Discrete_data<T>::count_t count_t;

            // Ctor
            explicit Locus_data(
                    const std::vector<T>&,  //the data
                    const T&);              //code of undefined/missing value
            explicit Locus_data(
                    typename std::vector<T>::const_iterator, //start
                    typename std::vector<T>::const_iterator, //end
                    const T&);

            // Inspectors
            T get_major() const { return major_; }
            T get_minor() const { return minor_; }
            T get_target() const { return target_; }
            T get_undef() const { return undef_; }
            bool hasMissings() const { return count_elem(undef_) > 0; }
            bool isPolymorph() const;
            size_t nValid() const { return this->size() - count_elem(undef_); }
            size_t nMiss() const { return count_elem(undef_); }
            double maf(Genetic_type type) const; //minor allele frequency

            // Modifiers
            void set_target(const T&); 
            void set_minor(const T&); 
            void set_major(const T&); 

            // Conversions
            // @parameter 'a' specifies how many alleles to merge
            Locus_data<uint>* merge_alleles_to_genotypes(uint a=2) const;
            Locus_data<uint>* as_numeric() const;

        private:
            void init();
            T major_;   //e.g. major marker
            T minor_;   //e.g. minor marker
            T target_;  //e.g. target/risk marker
            T undef_;   //e.g. undefined/missing code
    };


    // Locus_data<T> implementation
    // ========================================================================
    template<class T> inline Locus_data<T>::Locus_data(
            const std::vector<T>& v, const T& u) 
        : Discrete_data<T>(v), undef_(u) 
    {
        init();
    }
    template<class T> inline Locus_data<T>::Locus_data(
            typename std::vector<T>::const_iterator start, 
            typename std::vector<T>::const_iterator end, const T& u) 
        : Discrete_data<T>(start, end), undef_(u) 
    {
        init();
    }
    template<class T> inline void Locus_data<T>::init()
    {
        // determine minor and major allele:
        std::map<elem_t, count_t> m = this->unique_with_counts();
        m.erase(undef_); //undefined is not allowed
        minor_ =  (*min_element(m.begin(), m.end(), 
                    detail::comp_second<std::pair<T, int> >())).first;
        major_ =  (*max_element(m.begin(), m.end(), 
                    detail::comp_second<std::pair<T, int> >())).first;
        target_ = minor_; //default target is minor allele

        this->add_to_domain(undef_); //make undef_ always a part of domain
    }
    template<class T> inline bool Locus_data<T>::isPolymorph() const
    {
        return this->data_cardinality() > ( 1 + (size_t) hasMissings() );
    }

    template<class T> inline double Locus_data<T>::maf(Genetic_type type) const
    {
        if (nValid() == 0) {
            return 0.0;
        }
        if (type == genotype) {
            // derive maf via genotype
            std::map<elem_t, count_t> m = this->unique_with_counts();
            m.erase(undef_); //undefined does not count

            // First we need to transform the elements, which could be of type
            // int but also of type string, char, etc... into numeric type
            std::vector<count_t> vv; 
            vv.reserve(m.size());
            typename std::map<elem_t, count_t>::const_iterator itMap;
            for (itMap=m.begin(); itMap != m.end(); ++itMap) {
                vv.push_back(boost::lexical_cast<count_t>(itMap->first));
            }

            count_t max_sum = nValid()*(*max_element(vv.begin(), vv.end()));
            count_t sum = 0;
            typename Discrete_data<T>::unique_iterator it;
            for (it = this->unique_begin(); it != this->unique_end(); ++it) {
                elem_t g = it->first;    //the genotype
                bool isValid = g != undef_;
                if (isValid) {
                    //#alleles = genotype*(#occurrence)
                    sum += boost::lexical_cast<count_t>(g)*it->second; 
                }
            }
            return double(sum)/double(max_sum);
        }
        else { //or via haplotype, which is simple to compute
            return double(count_elem(minor_)) / double(nValid());
        }
    }

    template<class T> inline void Locus_data<T>::set_target(const T& x) 
    {
        if (x != undef_ && this->isInDomain(x)) {
            target_ = x;
        }
        else {
            throw(std::domain_error("Bad target allele/marker."));
        }
    }
    template<class T> inline void Locus_data<T>::set_minor(const T& x) 
    {
        if (x != undef_ && this->isInDomain(x)) {
            minor_ = x;
        }
        else {
            throw(std::domain_error("Bad minor allele/marker."));
        }
    }
    template<class T> inline void Locus_data<T>::set_major(const T& x) 
    {
        if (x != undef_ && this->isInDomain(x)) {
            major_ = x;
        }
        else {
            throw(std::domain_error("Bad major allele/marker."));
        }
    }

    template<class T> inline Locus_data<uint>* 
        Locus_data<T>::merge_alleles_to_genotypes(uint a) const
    {
        std::vector<uint> v; 
        // a=2 (default) means standard 2-allelic genotype
        v.reserve(this->size()/a);

        for (typename std::vector<T>::const_iterator it = this->begin(); 
                it+a-1 < this->end(); it+=a) {
            bool ok = (find(it, it+a, this->undef_) == it+a);
            if (ok) { 
                // genotype is the number of target/risk markers 
                v.push_back((uint) count(it, it+a, this->target_));
            }
            else {
                // Found one or more undefined alleles => undefined genotype 
                v.push_back(-1);
            }
        }
        //Locus_data<uint> ld(v, -1); //new Locus_data with undef set to -1 
        //return ld;
        return new Locus_data<uint>(v, -1); //new Locus_data with undef set to -1 
    }
    template<class T> inline Locus_data<uint>* Locus_data<T>::as_numeric() const
    {
        // We transform the values into integers by their position and hence
        // first need to put the values into a sequence such as vector
        std::vector<T> uu;
        uu.reserve(this->domain_cardinality()); 
        typedef Discrete_data<T> data_t;
        typedef typename std::map<
            typename data_t::elem_t, typename data_t::count_t > map_t;
        std::transform(this->unique_begin(), this->unique_end(), 
                std::back_inserter(uu), boost::bind(&map_t::value_type::first,_1) );

        std::vector<uint> v; 
        v.reserve(this->size());

        typedef typename std::vector<T>::const_iterator vec_iter;
        vec_iter start = uu.begin();
        vec_iter end = uu.end();
        for (vec_iter itData = this->begin(); itData != this->end(); ++itData) {
            bool ok = *itData != this->undef_;
            if (ok) {
                uint dx = std::distance(start, std::find(start, end, *itData));
                v.push_back(dx);
            }
            else {
                v.push_back(-1);    //undefined
            }
        }

        //Locus_data<uint> ld(v, -1); //new Locus_data with undef set to -1 
        //return ld;
        return new Locus_data<uint>(v, -1); //new Locus_data with undef set to -1 
    }
} // namespace Permory

#endif // include guard


