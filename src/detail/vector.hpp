// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_vector_hpp
#define permory_detail_vector_hpp

#include <map>
#include <valarray>
#include <vector>
#include <boost/dynamic_bitset.hpp>

#include "detail/functors.hpp"
#include "detail/pair.hpp"

namespace Permory { namespace detail {

    template<class T> boost::dynamic_bitset<>
        vector_to_bitset(const T& v)
    {
        boost::dynamic_bitset<> b(v.size());
        for (size_t i=0; i<v.size(); ++i) {
            b[i] = v[i];
        }
        return b;
    }

    // Specialization for Pair because of static typing. Needed to not produce
    // compilation error in constuctor of Perm_matrix (reason bitMat).
    template<class T>
    boost::dynamic_bitset<> vector_to_bitset(const std::vector<Pair<T> >& v)
    {
        throw std::runtime_error(
            "Bit arithmetics not allowed with quantitative phenotypes!");
    }

    template<class T> std::vector<T>
        bitset_to_vector(const boost::dynamic_bitset<>& b)
    {
        std::vector<T> v(b.size());
        for (size_t i=0; i<v.size(); ++i) {
            if (b[i])
                v[i] = 1;
        }
        return v;
    }

    template<class T> bool 
        sequence_is_sorted(const T& seq, bool decreasing=false)
    {
        typename T::const_iterator i;
        if (decreasing) {
            //if no value is less than its successor, it is decreasingly sorted 
            i = adjacent_find(seq.begin(), seq.end(), std::less<typename T::value_type>());
        }
        else { //increasing 
            i = adjacent_find(seq.begin(), seq.end(), std::greater<typename T::value_type>());
        }
        return (i == seq.end());
    }

    template<class T> std::valarray<T> 
        vector_to_valarray(const std::vector<T>& v)
    {
        typename std::valarray<T> va(&v[0], v.size());
        return va;
    }
    template<class T> std::vector<T> 
        valarray_to_vector(const std::valarray<T>& va)
    {
        typename std::vector<T> v(&va[0], va.size());
        return v;
    }

    template<class T> inline typename std::vector<T>::iterator max_pairwise(
            const std::vector<T>& v1, std::vector<T>& v2)
    {
            transform(v1.begin(), v1.end(), v2.begin(), v2.begin(), std::max<T>()); 
            return v2;
    }

    template<class T, class Compare> std::vector<T> regroup(
            const std::vector<T> v,     //will be reordered according to groups
            const std::vector<int> g)   //groups
    {
        // Two examples:
        // =========
        // vector<int> v = 1 2 3 4 5 6
        // vector<int> g = 1 1 0 0 4 4
        // 1) regroup<int, less<int> >(v, g):
        //                 3 4 1 2 5 6   
        // 2) regroup<int, greater<int> >(v, g): 
        //                 5 6 1 2 3 4
        typedef std::map<int, std::vector<int>, Compare> index_map;
        index_map m;
        for (int i=0; i<g.size(); i++) {
            m[g[i]].push_back(i);
        }
        std::vector<T> vv; 
        vv.reserve(v.size());
        for (typename index_map::const_iterator i=m.begin(); i!=m.end(); i++) {
            BOOST_FOREACH(int ii, i->second) { 
                vv.push_back(v[ii]); 
            }
        }
        return vv;
    }

} //namespace detail
} //namespace Permory

#endif

