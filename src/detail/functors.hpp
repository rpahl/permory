// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_functors_hpp
#define permory_detail_functors_hpp

#include <functional>
#include <deque>

#include "detail/config.hpp"

namespace Permory { namespace detail {
    // Predicates
    // ========================================================================
    // Compares two pairs only based on their second argument
    template<class T> class comp_second : 
        public std::binary_function<T, T, bool>
    {
        public:
            bool operator()(const T& x, const T& y) const {
                return x.second < y.second;
            }
    };

    // True if the second argument of a pair is greater than some value x
    template<class T1, class T2> class greater_than_second : 
        public std::unary_function<T2, bool> 
    {
        public:
            explicit greater_than_second(const T2& x) : arg2(x) {}
            bool operator()(const T1& x) const { return x.second > arg2; }
        private:
            T2 arg2;
    };

    // Operators
    // =======================================================================

    // Concatenates two deques into one.
    template<class T> struct deque_concat :
        public std::binary_function< std::deque<T>, std::deque<T>, std::deque<T> >
    {
            std::deque<T> operator()(const std::deque<T>& v1, const std::deque<T>& v2) {
                std::deque<T> result;
                copy(v1.begin(), v1.end(), back_inserter(result));
                copy(v2.begin(), v2.end(), back_inserter(result));
                return result;
            }
    };

} // namespace detail
} // namespace Permory

#endif // include guard
