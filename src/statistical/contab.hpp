// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_contab_hpp
#define permory_contab_hpp

#include <iostream>
#include <set>
#include <valarray>
#include <vector>

#include "detail/exception.hpp"
namespace Permory { namespace statistic {
    // Contingency table with R and C being the number of rows and columns, resp
    template<uint R, uint C> class Con_tab {
        public:
            typedef std::vector<uint>::const_iterator uint_iter;
            // Ctors
            Con_tab(); 
            Con_tab(uint_iter, uint_iter, uint_iter, uint_iter);

            // Inspectors
            uint nrow() const { return R; }
            uint ncol() const { return C; }
            uint rowsum(uint x) const;
            uint colsum(uint x) const;
            const uint& operator()(uint r, uint c) const { return tab_[r][c]; }
            uint* operator[](uint r) { return tab_[r]; }
            uint at(const uint r, const uint c) const; //range checked

            // Modifier
            void assign(const uint r, const uint c, const uint val); //range checked
            // Conversions
            void print() const; 
        private:
            uint tab_ [R] [C]; 
    };

    template<uint R, uint C> inline Con_tab<R, C>::Con_tab()
    {
        for (uint i=0; i<R; ++i) {
            for (uint j=0; j<C; ++j)
                tab_[i][j] = 0;
        }
    }
    template<uint R, uint C> inline Con_tab<R,C>::Con_tab(
            uint_iter start1, uint_iter end1, 
            uint_iter start2, uint_iter end2)
    {
        this->Con_tab(); //need to initialize with 0
#ifdef NDEBUG
        while (start1 != end1 && start2 != end2) {
            tab_[*start1++][*start2++]++;
        }
#else
        while (start1 != end1 && start2 != end2) {
            uint a = this->at(*start1, *start2);
            tab_.assign(*start1++, *start2++, a+1);
        }
#endif
    }

    /*
       template<class T1, class T2> template<uint R, uint C> inline 
       Con_tab<R,C>::Con_tab(
       std::vector<T1>::const_iterator start1,
       std::vector<T1>::const_iterator end1,
       std::vector<T2>::const_iterator start2,
       std::vector<T2>::const_iterator end2,
       )
       */

    template<uint R, uint C> inline uint Con_tab<R, C>::rowsum(uint x) const 
    {
        uint sum = 0;
        for (uint j=0; j<C; ++j) {
            sum += tab_[x][j];
        }
        return sum;
    }

    template<uint R, uint C> inline uint Con_tab<R, C>::colsum(uint x) const 
    {
        uint sum = 0;
        for (uint i=0; i<R; ++i) {
            sum += tab_[i][x];
        }
        return sum;
    }

    template<uint R, uint C> inline void Con_tab<R, C>::assign(
            const uint r, 
            const uint c, 
            const uint val) 
    {
        if (val < 0) {
            throw std::runtime_error("Con_tab: set negative value.");
        }
        if (!(r >= 0 && r < R && c >= 0 && c < C)) {
            throw std::runtime_error("Con_tab: Out of bounds.");
        }
        tab_[r][c] = val;
    }

    template<uint R, uint C> inline uint Con_tab<R, C>::at(
            const uint r, const uint c) const 
    {
        if (r >= 0 && r < R && c >= 0 && c < C) {
            return tab_[r][c];
        }
        else {
            throw std::runtime_error("Con_tab: Out of bounds.");
        }
    }

    template<uint R, uint C> inline void Con_tab<R, C>::print() const 
    {
        for (uint i=0; i<R; ++i) {
            for (uint j=0; j<C; ++j) {
                std::cout << tab_[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // ========================================================================
    // Non-member functions
    template<uint L> void fill_2xL_tabs(
            std::vector<Con_tab<2, L> >& tabs,
            const std::vector<std::valarray<uint> >& r,
            const std::vector<std::valarray<uint> >& s)
    {
#ifdef NDEBUG
        for (size_t j=0; j<L; ++j) {                //for each column
            for (size_t i=0; i<tabs.size(); ++i) {  //for each con tab
                tabs[i][0][j] = r[j][i]; //cases r[j]
                tabs[i][1][j] = s[j][i]; //controls s[j]
            }
        }
#else
        assert (r.size() == s.size());
        assert (r.size() == L);
        try {
            for (size_t j=0; j<L; ++j) {                //for each column
                for (size_t i=0; i<tabs.size(); ++i) {  //for each con tab
                    tabs[i].assign(0, j, r[j][i]); 
                    tabs[i].assign(1, j, s[j][i]);
                }
            }
        }
        catch (std::exception& e) {
            std::cerr << e.what() <<  std::endl;
            exit(-1);
        }
#endif
    }

    template<class T1, class T2, uint R, uint C> Con_tab<R,C>* make_Con_tab(
            typename std::vector<T1>::const_iterator start1, 
            typename std::vector<T1>::const_iterator end1, 
            typename std::vector<T2>::const_iterator start2, 
            typename std::vector<T2>::const_iterator end2)
    {
        // Filter out only the unique elements
        std::set<T1> s1(start1, end1);
        std::set<T2> s2(start2, end2);
        /* XXX
        PRINT(s1.size());
        PRINT(s2.size());
        detail::print_seq(s1.begin(), s1.end(), "", ",");
        detail::print_seq(s2.begin(), s2.end(), "", ",");
        */

        if (s1.size() > R || s2.size() > C) {
            throw std::out_of_range("Exceeding dimension in contingency table.");
        }
        std::vector<T1> v1(s1.begin(), s1.end());
        std::vector<T2> v2(s2.begin(), s2.end());

        Con_tab<R,C>* ctab = new Con_tab<R,C>(); //initialized with all entries to 0
        while (start1 != end1 && start2 != end2) {
            size_t i1 = std::distance(v1.begin(), 
                    find(v1.begin(), v1.end(), *start1++));
            size_t i2 = std::distance(v2.begin(), 
                    find(v2.begin(), v2.end(), *start2++));
            (*ctab)[i1][i2]++;
        }
        return ctab;
    }

} // namespace statistic
} // namespace Permory

#endif // include guard

