/**
 * @author Roman Pahl
 */

#ifndef permory_snp_data_hpp
#define permory_snp_data_hpp

#include <numeric>

#include "config.hpp"
#include "discretedata.hpp"
#include "helper/tokenizer.hpp" //copy_token
#include "helper/functors.hpp" //pair_comp_2nd
#include "locus.hpp"

namespace Permory 
{
    template<class T> class Locus_data : public Discrete_data<T> {
        public:
            typedef typename Discrete_data<T>::elem_type elem_type;
            typedef typename Discrete_data<T>::count_type count_type;
            enum Data_type {haplotype, genotype};

            // Ctor
            explicit Locus_data(
                    const std::vector<T>&, 
                    const T&, boost::shared_ptr<Locus>, Data_type);
            explicit Locus_data(
                    typename std::vector<T>::const_iterator, //start
                    typename std::vector<T>::const_iterator, //end
                    const T&, boost::shared_ptr<Locus>, Data_type);
            explicit Locus_data(
                    Tokenizer&, 
                    const T&, boost::shared_ptr<Locus>, Data_type);

            // Inspectors
            T get_major() const { return major_; }
            T get_minor() const { return minor_; }
            T get_target() const { return target_; }
            T get_undef() const { return undef_; }
            Data_type get_type() const { return type_; }
            bool hasMissings() const { return count_elem(undef_) > 0; }
            bool isPolymorphic() const;
            int nValid() const { return this->size() - count_elem(undef_); }
            int nMiss() const { return count_elem(undef_); }
            double maf() const; //minor allele frequency

            // Modifiers
            void set_target(const T&); 
            void set_minor(const T&); 
            void set_major(const T&); 
            boost::shared_ptr<Locus> get_locus() { return loc_; }
            
            // Conversions
            Locus_data<int> make_genotypic(int a=2) const;

        private:
            void init();
            T major_;    //e.g. major marker
            T minor_;    //e.g. minor marker
            T target_;   //e.g. target/risk marker
            T undef_;    //e.g. undefined/missing code
            Data_type type_; //haplo- or genotype
            boost::shared_ptr<Locus> loc_;  
    };


    // ========================================================================
    // Locus_data<T> implementation
    template<class T> inline Locus_data<T>::Locus_data(
            const std::vector<T>& v, 
            const T& undef, boost::shared_ptr<Locus> loc, Data_type dt)
        : Discrete_data<T>(v), undef_(undef), loc_(loc), type_(dt)
    {
        init();
    }
    template<class T> inline Locus_data<T>::Locus_data(
            typename std::vector<T>::const_iterator start, 
            typename std::vector<T>::const_iterator end,
            const T& undef, boost::shared_ptr<Locus> loc, Data_type dt)
        : Discrete_data<T>(start, end), undef_(undef), loc_(loc), type_(dt)
    {
        init();
    }
    template<class T> inline Locus_data<T>::Locus_data(
            Tokenizer& tok, 
            const T& undef, boost::shared_ptr<Locus> loc, Data_type dt)
        : Discrete_data<T>(tok), undef_(undef), loc_(loc), type_(dt)
    {
        init();
    }
    template<class T> inline void Locus_data<T>::init()
    {
        // determine minor and major allele:
        std::map<elem_type, count_type> m = this->unique_with_counts();
        m.erase(undef_); //undefined is not allowed
        minor_ =  (*min_element(m.begin(), m.end(), 
                    comp_second<std::pair<T, int> >())).first;
        major_ =  (*max_element(m.begin(), m.end(), 
                    comp_second<std::pair<T, int> >())).first;
        target_ = minor_; //default target is minor allele

        this->add_to_domain(undef_); //make undef_ always a part of domain
    }
    template<class T> inline bool Locus_data<T>::isPolymorphic() const
    {
        return this->data_cardinality() > ( 1 + (size_t) hasMissings() );
    }

    template<class T> inline double Locus_data<T>::maf() const
    {
        if (nValid() == 0)
            return 0.0;
        if (type_ == genotype) {
            std::map<elem_type, count_type> m = this->unique_with_counts();
            m.erase(undef_); //undefined does not count
            elem_type max_sum = 
                nValid()*(*max_element(m.begin(), m.end())).first;

            elem_type sum = 0;
            for(typename Discrete_data<T>::counts_iterator 
                    i=this->counts_begin(); i!=this->counts_end(); i++) {
                elem_type e = i->second; //the genotype
                bool isValid = e != undef_;
                if (isValid)
                    sum += e*(i->first); //#alleles = genotype*(#occurrence)
            }
            return double(sum)/double(max_sum);
        }
        else { //haplotype
            return double(count_elem(minor_)) / double(nValid());
        }
    }

    template<class T> inline void Locus_data<T>::set_target(const T& x) 
    {
        //bool isElement = m[x] > 0;
        if (x != undef_ && this->isInDomain(x)) {
            target_ = x;
        }
        else {
            std::cerr << "Error: target allele is undefined.\n";
            exit(-1);
        }
    }
    template<class T> inline void Locus_data<T>::set_minor(const T& x) 
    {
        if (x != undef_ && this->isInDomain(x)) {
            minor_ = x;
        }
        else {
            std::cerr << "Error: minor allele is undefined.\n";
            exit(-1);
        }
    }
    template<class T> inline void Locus_data<T>::set_major(const T& x) 
    {
        if (x != undef_ && this->isInDomain(x)) {
            major_ = x;
        }
        else {
            std::cerr << "Error: major allele is undefined.\n";
            exit(-1);
        }
    }

    template<class T> inline Locus_data<int> Locus_data<T>::make_genotypic(int a) const
    {
        std::vector<int> v; 
        // a=2 (default) means standard 2-allelic genotype
        v.reserve(this->size()/a);

        for (typename std::vector<T>::const_iterator it = this->begin(); 
                it+a-1 < this->end(); it+=a) {
            if (find(it, it+a, this->undef_) == it+a) { //check for missing
                // genotype is the number of target/risk markers 
                v.push_back((count_type) count(it, it+a, this->target_));
            }
            else {
                // One or more undefined alleles => undefined genotype 
                v.push_back(-1);
            }
        }
        Locus_data<int> s(v, -1, genotype); //new Locus_data with undef set to -1 
        return s;
    }
} // namespace Permory

#endif // include guard


