// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_bitset_hpp
#define permory_detail_bitset_hpp

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "detail/config.hpp"

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

} // namespace detail
} // namespace Permory

#endif // include guard


