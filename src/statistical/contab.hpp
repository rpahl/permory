/**
 * @author Roman Pahl
 */

#ifndef permory_contab_hpp
#define permory_contab_hpp

#include <iostream>
#include <valarray>

#include "exception.hpp"

namespace Permory 
{
    // Contingency table with R and C being the number of rows and columns, resp
    template<int R, int C> class Con_tab {
        public:
            typedef std::vector<int>::const_iterator int_iter;
            // Ctors
            Con_tab(); 
            Con_tab(int_iter, int_iter, int_iter, int_iter);

            // Inspectors
            int nrow() const { return R; }
            int ncol() const { return C; }
            int rowsum(int x) const;
            int colsum(int x) const;
            const int& operator()(int r, int c) const { return tab_[r][c]; }
            int* operator[](int r) { return tab_[r]; }
            int at(const int r, const int c) const //range checked
                throw (std::runtime_error);

            // Modifier
            void assign(const int r, const int c, const int val) //range checked
                throw (std::runtime_error); 
            // Conversions
            void print() const; 
        private:
            int tab_ [R] [C]; 
    };

    template<int R, int C> inline Con_tab<R, C>::Con_tab()
    {
        for (size_t i=0; i<R; ++i) {
            for (size_t j=0; j<C; ++j)
                tab_[i][j] = 0;
        }
    }
    template<int R, int C> inline Con_tab<R, C>::Con_tab(
            int_iter start1, int_iter end1, 
            int_iter start2, int_iter end2)
    {
        this->Con_tab(); //need to initialize with 0
#ifdef NDEBUG
        while (start1 != end1 && start2 != end2) 
            tab_[*start1++][*start2++]++;
#else
        try {
            while (start1 != end1 && start2 != end2) {
                int a = this->at(*start1, *start2);
                tab_.assign(*start1++, *start2++, a+1);
            }
        }
        catch (std::exception& e) {
            std::cerr << e.what() <<  "\n";
            exit(-1);
        }
#endif
    }


    template<int R, int C> inline int Con_tab<R, C>::rowsum(int x) const 
    {
        int sum = 0;
        for (int j=0; j<C; ++j)
            sum += tab_[x][j];
        return sum;
    }

    template<int R, int C> inline int Con_tab<R, C>::colsum(int x) const 
    {
        int sum = 0;
        for (int i=0; i<R; ++i)
            sum += tab_[i][x];
        return sum;
    }

    template<int R, int C> inline void Con_tab<R, C>::assign(
            const int r, 
            const int c, 
            const int val) throw (std::runtime_error)
    {
        if (val < 0) 
            throw std::runtime_error("Con_tab: set negative value.");
        if (!(r >= 0 && r < R && c >= 0 && c < C))
            throw std::runtime_error("Con_tab: Out of bounds.");
        tab_[r][c] = val;
    }

    template<int R, int C> inline int Con_tab<R, C>::at(
            const int r, const int c) const throw (std::runtime_error)
    {
        if (r >= 0 && r < R && c >= 0 && c < C)
            return tab_[r][c];
        else
            throw std::runtime_error("Con_tab: Out of bounds.");
    }

    template<int R, int C> inline void Con_tab<R, C>::print() const 
    {
        for (int i=0; i<R; ++i) {
            for (int j=0; j<C; ++j)
                std::cout << tab_[i][j] << " ";
            std::cout << std::endl;
        }
    }

    // ========================================================================
    // Non-member functions
    template<int L> void fill_2xL_tabs(
            std::vector<Con_tab<2, L> >& tabs,
            const std::vector<std::valarray<int> >& r,
            const std::vector<std::valarray<int> >& s)
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
            std::cerr << e.what() <<  "\n";
            exit(-1);
        }
#endif
    }

} // namespace Permory

#endif // include guard

