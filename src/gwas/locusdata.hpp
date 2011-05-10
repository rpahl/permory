// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_snp_data_hpp
#define permory_snp_data_hpp

#include <numeric>

#include<boost/bind.hpp>
#include<boost/lexical_cast.hpp>

#include "detail/config.hpp"
#include "detail/enums.hpp" //Marker_type
#include "discretedata.hpp"
#include "locus.hpp"

namespace Permory { namespace gwas {

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
            double maf(detail::Marker_type) const; //*m*inor *a*llele *f*requency
            void get_data_without_missings(std::vector<T>*) const;

            // Modifiers
            void set_target(const T&); 
            void set_minor(const T&); 
            void set_major(const T&); 

            // Conversions
            // @parameter 'a' specifies how many alleles to merge
            Locus_data<T> condense_alleles_to_genotypes(uint a=2) const;
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
        minor_ =  (*std::min_element(m.begin(), m.end(), 
                    detail::comp_second<std::pair<elem_t, count_t> >())).first;
        major_ =  (*std::max_element(m.begin(), m.end(), 
                    detail::comp_second<std::pair<elem_t, count_t> >())).first;
        target_ = minor_; //default target is minor allele

        this->add_to_domain(undef_); //make undef_ always a part of domain
    }
    template<class T> inline bool Locus_data<T>::isPolymorph() const
    {
        return this->data_cardinality() > ( 1 + (size_t) hasMissings() );
    }

    template<class T> inline double Locus_data<T>::maf(detail::Marker_type mt) const
    {
        if (nValid() == 0) {
            return 0.0;
        }
        if (mt == detail::genotype) {
            // derive maf via genotype
            std::map<elem_t, count_t> m = this->unique_with_counts();
            m.erase(undef_); //undefined does not count

            // First we need to transform the domain, which could be of type
            // int but also of type string, char, etc... into countable type
            std::vector<count_t> vv; 
            vv.reserve(m.size());
            typename std::map<elem_t, count_t>::const_iterator itMap;
            for (itMap=m.begin(); itMap != m.end(); ++itMap) {
                try {
                    vv.push_back(boost::lexical_cast<count_t>(itMap->first));
                }
                catch (const boost::bad_lexical_cast& e) {
                    std::string s = "Bad data entry '";
                    s.append(boost::lexical_cast<std::string>(itMap->first));
                    s.append("' because it is not part of the genotype domain.\n");
                    throw std::domain_error(s);
                }
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
            double maf = double(sum)/double(max_sum);
            return maf > 0.5 ? 1.0 - maf : maf;
        }
        else { //or via alleles, which is simple to compute
            double maf = double(count_elem(minor_)) / double(nValid());
            return maf > 0.5 ? 1.0 - maf : maf;
        }
    }

    template<class T> inline void
    Locus_data<T>::get_data_without_missings(std::vector<T>* result) const
    {
        result->reserve(nValid());
        for (typename Locus_data<T>::const_iterator it = this->begin();
             it != this->end();
             ++it) {
            if (not (*it == undef_)) {
                result->push_back(*it);
            }
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

    template<class T> inline Locus_data<T> 
        Locus_data<T>::condense_alleles_to_genotypes(uint a) const
    {
        T undef_condense = boost::lexical_cast<T>('?');
        std::vector<T> v; 
        // a=2 (default) means standard 2-allelic genotype
        v.reserve(this->size()/a);

        for (typename std::vector<T>::const_iterator it = this->begin(); 
                it+a-1 < this->end(); it+=a) {
            bool ok = (find(it, it+a, this->undef_) == it+a);
            if (ok) { 
                // genotype is the number of target/risk markers 
                size_t cnt = std::count(it, it+a, this->target_);
                v.push_back(boost::lexical_cast<T>(cnt));
            }
            else {
                // Found one or more undefined alleles => undefined genotype 
                v.push_back(undef_condense);
            }
        }
        return Locus_data<T>(v, undef_condense);
    }

    template<> inline Locus_data<uint>* Locus_data<char>::as_numeric() const
    {
        std::vector<uint> v; 
        v.reserve(this->size());
        typedef std::vector<char>::const_iterator vec_iter;
        for (std::vector<char>::const_iterator itChar = this->begin(); 
                itChar != this->end(); ++itChar) {
            v.push_back(uint(*itChar) - 48);
        }
        char c = uint(this->undef_) - 48;
        return new Locus_data<uint>(v, c); 
    }

    /* not used
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
        return new Locus_data<uint>(v, -1); //new Locus_data with undef set to -1 
    }
    */
} // namespace gwas
} // namespace Permory

#endif // include guard


