

/**
 * @author Roman Pahl
 */

#ifndef permory_helper_vector_hpp
#define permory_helper_vector_hpp

#include <map>
#include <valarray>
#include <vector>

#include "helper/functors.hpp"

namespace Permory {

    template<class T> bool 
        vector_is_sorted(
            const std::vector<T>& v, 
            bool decreasing=false)
    {
        typename std::vector<T>::const_iterator i;
        if (decreasing) {
            //if no value is less than its successor, it is decreasingly sorted 
            i = adjacent_find(v.begin(), v.end(), std::less<T>());
        }
        else { //increasing 
            i = adjacent_find(v.begin(), v.end(), std::greater<T>());
        }
        return (i == v.end());
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
            const std::vector<T> v,
            const std::vector<int> groups) 
    {
        // Examples:
        // =========
        // vector<int> v = 1 2 3 4 5 6
        // vector<int> g = 1 1 0 0 4 4
        // 1) regroup<int, less<int> >(v, g):
        //                 3 4 1 2 5 6   
        // 2) regroup<int, greater<int> >(v, g): 
        //                 5 6 1 2 3 4
        typedef std::map<int, std::vector<int>, Compare> index_map;
        index_map m;
        for (int i=0; i<groups.size(); i++) {
            m[groups[i]].push_back(i);
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

} //namespace Permory

#endif

