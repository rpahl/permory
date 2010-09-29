// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_discretedata_hpp
#define permory_discretedata_hpp

#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "detail/functors.hpp" //predicates
#include "detail/vector.hpp"

namespace Permory 
{
    template<class T> class Discrete_data { 
        public:
            typedef T elem_t;
            typedef uint count_t;

            // Iterator pass through
            typedef typename std::vector<T>::const_iterator const_iterator;
            typedef typename std::map<T, uint>::const_iterator unique_iterator;
            typedef typename std::map<uint, T>::const_iterator counts_iterator;
            const_iterator begin() const { return data_.begin(); }
            const_iterator end() const { return data_.end(); }
            unique_iterator unique_begin() const { return unique_.begin(); } 
            unique_iterator unique_end() const { return unique_.end(); }
            counts_iterator counts_begin() const { return counts_.begin(); } 
            counts_iterator counts_end() const { return counts_.end(); }

            // Ctor
            explicit Discrete_data(const std::vector<T>&);  
            explicit Discrete_data(const_iterator, const_iterator);

            // Inspector
            const T& operator[](const size_t pos) const { return data_[pos]; }
            size_t size() const { return data_.size(); }
            size_t domain_cardinality() const { return unique_.size(); }
            size_t data_cardinality() const; 
            std::map<T, count_t> unique_with_counts() const { return unique_; } 
            bool isInDomain(const T&) const;
            count_t count_elem(const T&) const;
            void print(); //for debugging

            // Modifier
            template<class Compare> void regroup(const std::vector<uint> v);
            void add_to_domain(const std::set<T>& s);
            void add_to_domain(const T& x);

            // Conversion
            Discrete_data<T> mask(const std::vector<bool>&);
        private:
            void init();
            std::vector<T> data_;        
            std::map<elem_t, count_t> unique_;//unique elements with counts
            std::multimap<count_t, elem_t> counts_;//and vice versa
    };

    // ========================================================================
    // Discrete_data<T> implementation
    template<class T> inline Discrete_data<T>::Discrete_data(const std::vector<T>& d) 
        : data_(d) 
    {
        init();
    }
    template<class T> inline Discrete_data<T>::Discrete_data(
            const_iterator start, const_iterator end)
        : data_(start, end)
    {
        init();
    }
    template<class T> inline void Discrete_data<T>::init() 
    {
        BOOST_FOREACH(T i, this->data_) {
            unique_[i]++;
        }
        for (typename std::map<T, uint>::iterator i=unique_.begin(); 
                i!=unique_.end(); i++)
            counts_.insert(std::make_pair(i->second, i->first));
    }

    template<class T> template<class Compare> inline void 
        Discrete_data<T>::regroup(const std::vector<uint> v) 
    {
        // Reorder data according to group indices contained in v. For more 
        // details see the 'regroup' function in detail/vector.hpp.
        this-> data_ = regroup<T, Compare>(data_, v); 
    }

    template<class T> inline size_t Discrete_data<T>::data_cardinality() const
    {
        // Count each unique element that appears at least once in the data
        return count_if(unique_.begin(), unique_.end(), 
                detail::greater_than_second<std::pair<T, uint>, uint>(0));
    }

    template<class T> inline bool Discrete_data<T>::isInDomain(const T& x) const
    {
        return unique_.find(x) != unique_.end();
    }

    template<class T> inline uint Discrete_data<T>::count_elem(const T& x) const
    {
        typename std::map<T, uint>::const_iterator it = unique_.find(x);
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
        for (uint i=0; i<std::min(data_.size(), b.size()); ++i)
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


