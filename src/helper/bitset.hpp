/**
 * @author Roman Pahl
 */

#ifndef permory_helper_bitset_hpp
#define permory_helper_bitset_hpp

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "config.hpp"

namespace Permory 
{

    template<class T> boost::dynamic_bitset<>
        vector_to_bitset(const T& v)
    {
        boost::dynamic_bitset<> b(v.size());
        for (size_t i=0; i<v.size(); ++i) {
            b[i] = v[i];
        }
        return b;
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

} // namespace Permory

#endif // include guard


