/**
 * @author Roman Pahl
 */

#ifndef permory_discretedata_hpp
#define permory_discretedata_hpp

#include <iostream>
#include <map>
#include <set>
#include <vector>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

#include "helper/algorithm.hpp" //copy_token
#include "helper/functors.hpp" //predicates
#include "helper/vector.hpp"

namespace Permory 
{
    using namespace ::boost::multi_index;

    template<class T> class Discrete_data { 
        public:
            typedef T elem_type;
            typedef int count_type;

            // Iterator pass through
            typedef typename std::vector<T>::const_iterator const_iterator;
            typedef typename std::map<T, int>::const_iterator unique_iterator;
            typedef typename std::map<int, T>::const_iterator counts_iterator;
            const_iterator begin() const { return data_.begin(); }
            const_iterator end() const { return data_.end(); }
            unique_iterator unique_begin() const { return unique_.begin(); } 
            unique_iterator unique_end() const { return unique_.end(); }
            counts_iterator counts_begin() const { return counts_.begin(); } 
            counts_iterator counts_end() const { return counts_.end(); }

            // Ctor
            explicit Discrete_data(const std::vector<T>&);  //direct copy
            explicit Discrete_data(const_iterator, const_iterator);//direct copy
            explicit Discrete_data(Tokenizer&);//tokenized string

            // Inspector
            const T& operator[](const size_t pos) const { return data_[pos]; }
            size_t size() const { return data_.size(); }
            size_t domain_cardinality() const { return unique_.size(); }
            size_t data_cardinality() const; 
            std::map<T, count_type> unique_with_counts() const { return unique_; } 
            bool isInDomain(const T&) const;
            count_type count_elem(const T&) const;
            void print(); //for debugging

            // Modifier
            void assign(Tokenizer&);
            template<class Compare> void regroup(const std::vector<int> v);
            void add_to_domain(const std::set<T>& s);
            void add_to_domain(const T& x);

            // Conversion
            Discrete_data<T> mask(const std::vector<bool>&);
        private:
            void init();
            std::vector<T> data_;        
            std::map<elem_type, count_type> unique_;//unique elements with counts
            std::multimap<count_type, elem_type> counts_;//vice versa
    };

    // ========================================================================
    // Discrete_data<T> implementation
    template<class T> inline Discrete_data<T>::Discrete_data(const std::vector<T>& d) 
        : data_(d) 
    {
        init();
    }
    template<class T> inline Discrete_data<T>::Discrete_data(const_iterator start, const_iterator end)
        :data_(start, end)
    {
        init();
    }
    template<class T> inline Discrete_data<T>::Discrete_data(Tokenizer& tok) 
    {
        copy_token<T>(tok, data_);
        //data_ = copy_token<T>(tok);
        init();
    }
    template<class T> inline void Discrete_data<T>::assign(Tokenizer& tok) 
    {
        unique_.clear();
        counts_.clear();
        data_.clear();
        copy_token<T>(tok, data_);
        init();
    }

    template<class T> inline void Discrete_data<T>::init() 
    {
        BOOST_FOREACH(T i, this->data_) {
            unique_[i]++;
        }
        for (typename std::map<T, int>::iterator i=unique_.begin(); 
                i!=unique_.end(); i++)
            counts_.insert(std::make_pair(i->second, i->first));
    }


    template<class T> template<class Compare> inline void 
        Discrete_data<T>::regroup(const std::vector<int> v) 
    {
        this-> data_ = regroup<T, Compare>(data_, v); //see helper/vector.hpp
    }

    template<class T> inline size_t Discrete_data<T>::data_cardinality() const
    {
        return count_if(unique_.begin(), unique_.end(), 
                greater_than_second<std::pair<T, int>, int>(0));
    }

    template<class T> inline bool Discrete_data<T>::isInDomain(const T& x) const
    {
        return unique_.find(x) != unique_.end();
    }

    template<class T> inline int Discrete_data<T>::count_elem(const T& x) const
    {
        typename std::map<T, int>::const_iterator it = unique_.find(x);
        return it != unique_.end() ? it->second : 0;
    }

    template<class T> inline void Discrete_data<T>::add_to_domain(
            const std::set<T>& s) 
    {
        BOOST_FOREACH(T x, s) {
            if (!isInDomain(x)) {
                unique_[x] = 0;
                counts_.insert(std::make_pair(0, x));
            }
        }
    }
    template<class T> inline void Discrete_data<T>::add_to_domain(const T& x) 
    {
        if (!isInDomain(x)) {
            unique_[x] = 0;
            counts_.insert(std::make_pair(0, x));
        }
    }

    template<class T> inline Discrete_data<T> Discrete_data<T>::mask(
            const std::vector<bool>& b)
    {
        std::vector<T> v;
        for (int i=0; i<std::min(data_.size(), b.size()); ++i)
            if (b[i]) v.push_back(data_[i]);
        return Discrete_data<T>(v);
    }

    template<class T> inline void Discrete_data<T>::print()
    {
        copy(data_.begin(), data_.end(), std::ostream_iterator<T>(std::cout," "));
        std::cout << std::endl;
    }
} // namespace Permory

#endif // include guard


