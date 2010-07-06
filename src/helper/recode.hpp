/**
 * @author Roman Pahl
 */

#ifndef permory_helper_recode_hpp
#define permory_helper_recode_hpp

#include <boost/dynamic_bitset.hpp>
#include <vector>

namespace Permory 
{
    // Build index code
    template<class T> std::vector<int> index_code(
            typename std::vector<T>::const_iterator start, 
            typename std::vector<T>::const_iterator end, 
            const T& val)
    {
        std::vector<int> vi;
        int i=0;
        while (start != end) {
            if (*start++ == val) { 
                vi.push_back(i); //store index
            }
            i++;
        }
        return vi;
    }

    // Build dummy code
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
    template<class T, class U> std::vector<U> dummy_code(
            typename std::vector<T>::const_iterator start, 
            typename std::vector<T>::const_iterator end, 
            const T& val)
    {
        std::vector<U> v(distance(start, end), 0);
        BOOST_FOREACH(U i, v) {
            if (*start++ == val) { 
                i = 1;
            }
        }
        return v;
    }
} // namespace Permory

#endif // include guard


