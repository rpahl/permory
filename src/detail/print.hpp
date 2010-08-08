// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_print_hpp
#define permory_detail_print_hpp

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include <map>
#include <valarray>
#include <vector>

#include "config.hpp"
#include "detail/exception.hpp"

using std::cout;
using std::endl;
using std::string;

namespace Permory { namespace detail {

    class Print {
        public:
            std::ostream& out;
            string l, r;
            Print(std::ostream& os, string l="", string r="") 
                : out(os), l(l), r(r) {}
            template <class T> void operator()(T const& x) { 
                out << l << x << r; }
    };

    template<class T> void print_seq( 
            const T& v, const string l="", const string r=" ", int len=-1)
    {
        if(len < 0 || len > v.size()) len = v.size();
        //std::ostream_iterator<class T::value_type> output(cout, del.c_str());
        //std::copy( v.begin(), v.begin()+len, output );
        std::for_each(v.begin(), v.begin()+len, Print(cout, l, r));
    }
    template<class InIt> void print_seq( 
            InIt first, InIt last, 
            const string l="", const string r=" ")
    {
        std::for_each(first, last, Print(cout, l, r));
    }

    template<class T, class U> void print_map( 
            const std::map<T, U>& m,
            const string del="", int len=-1)
    {
        if(len < 0 || len > m.size()) len = m.size();
        class std::map<T, U>::const_iterator it = m.begin();
        for (it; it!=m.end(); it++) 
            cout << "[" << (*it).first << "]" << (*it).second << " ";
    }

    template<class C> void print_vec( 
            const C& v,
            const string del="", 
            int len=-1)
    {
        if (v.size() == 0) 
            throw Dimension_error("Error in print_vec: vector is empty!\n");
        if(len < 0 || len > v.size()) len = v.size();
        for (size_t i=0; i<v.size(); ++i)
            cout << v[i] << del;
    }

    template<class C> void print_mat( 
            const C& m,
            const string del="", 
            int r=-1, 
            int c=-1)
    {
        if (m.size() == 0) 
            throw Dimension_error("Error in print_mat: matrix is empty!\n");

        if(r < 0 || r > m.nrow()) r = m.nrow();
        if(c < 0 || c > m.ncol()) c = m.ncol();
        for (int i=0; i<r; ++i) {
            for (int j=0; j<c; ++j) 
                cout << m[i][j] << del;
            cout << endl;
        }
    }


} //namespace detail
} //namespace Permory

#endif


