// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_permutation_recode_hpp
#define permory_permutation_recode_hpp

#include <vector>

#include <boost/dynamic_bitset.hpp>

namespace Permory { namespace permutation {
    
    //
    // For an input sequence s and some value x, this function returns a vector
    // v of indices such that for i=0,1,...,v.size()-1 : v[i] = x;
    // Example: 
    // s = 1 2 0 1 1 2 0 0 1 0
    // x = 0
    // v = index_code(s.start(), s.end(), x)
    // v = 2 6 7 9
    //
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

    //
    // Derive index code from a bitset, that is, return a vector of all indices
    // that are set in the bitset
    //
    std::vector<int> index_code(const boost::dynamic_bitset<>& b)
    {
        std::vector<int> v;
        v.reserve(b.count());
        for (int i=0; i<int(b.size()); ++i) {
            if (b.test(i)) {    //if (b[i] == 1)
                v.push_back(i);
            }
        }
        return v;
    }

    //
    // For an input sequence s and some value x, this function returns a bitset
    // b such that for i=0,1,...,b.size()-1 : b[i] = 1 iff v[i] = x;
    // Example: 
    // s = 1 2 0 1 1 2 0 0 1 0
    // x = 0
    // b = dummy_code(s.start(), s.end(), x)
    // b = 0 0 1 0 0 0 1 1 0 1
    //
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

    /* not used
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


