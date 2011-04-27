// Copyright (c) 2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_pair_hpp
#define permory_detail_pair_hpp

#include <utility>
#include <ostream>

namespace Permory { namespace detail {

    //
    // Specialization of std::pair for usage in PERMORY.
    //
    template<class T> struct Pair : public std::pair<T, T>
    {
        public:
            // Ctors
            Pair() : std::pair<T, T>() { }
            Pair(T x) : std::pair<T, T>(x, x) { }
            Pair(T f, T s) : std::pair<T, T>(f, s) { }
            Pair(const std::pair<T, T>& p) : std::pair<T, T>(p) { }

            // Operators
            Pair<T> operator+(const Pair<T>& other) const;
            Pair<T>& operator+=(const Pair<T>& other);
            Pair<T>& operator-=(const Pair<T>& other);
            Pair<T>& operator=(const Pair<T>& r);
            Pair<T>& operator=(const T& r);
    };

    // ========================================================================
    // Pair implementations
    template<class T> Pair<T> Pair<T>::operator+(const Pair<T>& other) const
    {
        return Pair(this->first + other.first, this->second + other.second);
    }

    template<class T> Pair<T>& Pair<T>::operator+=(const Pair<T>& other)
    {
        this->first += other.first;
        this->second += other.second;
        return *this;
    }

    template<class T> Pair<T>& Pair<T>::operator-=(const Pair<T>& other)
    {
        this->first -= other.first;
        this->second -= other.second;
        return *this;
    }

    template<class T> Pair<T>& Pair<T>::operator=(const Pair<T>& r)
    {
        this->first = r.first;
        this->second = r.second;
        return *this;
    }

    template<class T> Pair<T>& Pair<T>::operator=(const T& r)
    {
        this->first = r;
        this->second = r;
        return *this;
    }

    // None-Class Functions
    // ========================================================================

    //
    // Helper Funktion to create instances of Pair.
    template<class T> Pair<T> make_pair(T a, T b)
    {
        return Pair<T>(a, b);
    }

    //
    // Make instances of Pair printable to streams.
    template<class T>
    std::ostream& operator<<(std::ostream& o, const Pair<T>& x)
    {
        o << "(" << x.first << "," << x.second << ")";
        return o;
    }

} // namespace detail
} // namespace Permory
#endif

