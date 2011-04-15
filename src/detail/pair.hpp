// Copyright (c) 2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_pair_hpp
#define permory_detail_pair_hpp

#include <utility>
#include <ostream>

namespace Permory { namespace detail {

    template<class T1, class T2>
    std::pair<T1, T2> operator+(
                const std::pair<T1, T2>& a,
                const std::pair<T1, T2>& b) {
        return std::make_pair(a.first + b.first, a.second + b.second);
    }

    template<class T1, class T2>
    std::pair<T1, T2>& operator+=(
                std::pair<T1, T2>& a,
                const std::pair<T1, T2>& b) {
        a.first += b.first;
        a.second += b.second;
        return a;
    }

    template<class T1, class T2>
    std::ostream& operator<<(std::ostream& o, const std::pair<T1, T2>& x) {
        o << "(" << x.first << "," << x.second << ")";
        return o;
    }
} //namespace detail
} //namespace Permory

#endif

