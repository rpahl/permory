// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_permutation_recode_hpp
#define permory_permutation_recode_hpp

#include <vector>

#include <boost/dynamic_bitset.hpp>

namespace Permory { namespace permutation {
    
    // Index code (return by value)
    template<class T> std::vector<int> index_code(
            typename std::vector<T>::const_iterator start, 
            typename std::vector<T>::const_iterator end, 
            const T& val)
    {
        std::vector<int> v;
        v.reserve(std::distance(start, end));
        int i=0;
        while (start != end) {
            if (*start++ == val) { 
                v.push_back(i); //store index
            }
            i++;
        }
        return v;
    }

    // Index code (by reference)
    template<class T> void index_code(
            std::vector<int>& v,
            typename std::vector<T>::const_iterator start, 
            typename std::vector<T>::const_iterator end, 
            const T& val)
    {
        int i=0;
        while (start != end) {
            if (*start++ == val) { 
                v.push_back(i); //store index
            }
            i++;
        }
    }

    // Bitset dummy code
    template<class T> boost::dynamic_bitset<> dummy_code(
            typename std::vector<T>::const_iterator start, 
            typename std::vector<T>::const_iterator end, 
            const T& val)
    {
        boost::dynamic_bitset<> b(distance(start, end));
        int i=0;
        while (start != end) {
            if (*start++ == val) { 
                b.set(i);      
            }
            i++;
        }
        return b;
    }

    /*
    boost::dynamic_bitset<> dummy_code(
            std::vector<int>::const_iterator start, 
            std::vector<int>::const_iterator end,
            size_t len)
    {
        boost::dynamic_bitset<> b(len);
        while (start != end) {
            b.set(*start++);      
        }
        return b;
    }

    template<class T, class U> std::vector<U> dummy_code(
            typename std::vector<T>::const_iterator start, 
            typename std::vector<T>::const_iterator end, 
            const T& val)
    {
        std::vector<U> v(distance(start, end), 0);
        for (int i=0; i<v.size(); i++) {
            if (*start++ == val) { 
                v[i] = (U) 1;
            }
        }
        return v;
    }

    template<class T, class U> std::vector<U> dummy_code(
            std::vector<U>& v,
            typename std::vector<T>::const_iterator start, 
            typename std::vector<T>::const_iterator end, 
            const T& val)
    {
        for (int i=0; i<v.size(); i++) {
            if (*start++ == val) { 
                v[i] = 1;
            }
        }
    }
    */
} // namespace detail
} // namespace Permory

#endif // include guard


