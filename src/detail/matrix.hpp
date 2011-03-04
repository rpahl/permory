// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_matrix_hpp
#define permory_detail_matrix_hpp

#include <exception>
#include <iostream>
#include <map>
#include <valarray>
#include <vector>

namespace Permory { namespace detail {
    template<class T> class Matrix {
        public:
            typedef std::valarray<T> array_t;
            typedef std::valarray<array_t> matrix_t;

            // Ctors (throw if memory allocation fails)
            Matrix(size_t r=0, size_t c=0) throw(std::exception); 
            Matrix(size_t, const array_t&) throw(std::exception); 
            Matrix(Matrix<T>&) throw(std::exception); 

            // Inspection
            bool empty() const { return m_.size() == 0; }
            size_t size() const { return nrow()*ncol(); }
            size_t nrow() const { return m_.size(); }
            size_t ncol() const { return m_.size() > 0 ? m_[0].size() : 0; }
            const array_t& operator[](size_t i) const { return m_[i]; }

            // Conversion
            array_t row_apply(int (std::valarray<T>::*f)() const) const;
            array_t apply_sum(size_t) const;
            array_t apply_sum(size_t, const std::valarray<size_t>&) const;
            array_t apply_max(size_t) const;
            array_t apply_max(size_t, const std::valarray<size_t>&) const;
            array_t apply_min(size_t) const;
            array_t apply_min(size_t, const std::valarray<size_t>&) const;
            matrix_t transpose(); //not implemented

            // Modification
            array_t& operator[](size_t i) { return m_[i]; }
            Matrix<T>& operator=(const T&);
            Matrix<T>& operator=(const array_t&);
            Matrix<T>& operator=(const Matrix<T>&);
            Matrix<T>& resize(size_t, size_t); //data will be lost

        private:
            matrix_t m_;
    };

    // ========================================================================
    // Matrix implementation
    template<class T> inline Matrix<T>::Matrix(size_t r, size_t c)
        throw(std::exception)
        { 
            try {
                std::valarray<T> a(c);
                m_.resize(r, a); //mat with dim=[r][c]
            }
            catch(std::exception& e) {
                std::cerr << "Error: memory allocation for Matrix(" <<
                    r <<", "<< c <<") failed (" << e.what() << ").\n";
                exit(-1);
            }
        }

    template<class T> inline Matrix<T>::Matrix(
            size_t r, 
            const std::valarray<T>& a) throw(std::exception)
    { 
        try {
            m_.resize(r, a); //mat with dim=[r][a.size()]
        }
        catch(std::exception& e) {
            std::cerr << "Error: memory allocation for Matrix(" <<
                r <<", "<< a.size() <<") failed (" << e.what() << ").\n";
            exit(-1);
        }
    }

    template<class T> inline Matrix<T>::Matrix(Matrix<T>& m) throw(std::exception)
    { 
        try {
            m_.resize(m.nrow(), m.ncol()); 
        }
        catch(std::exception& e) {
            std::cerr << "Error: memory allocation for Matrix(" <<
                m.nrow() <<", "<< m.ncol() <<") failed (" << e.what() << ").\n";
            exit(-1);
        }
        for (size_t i=0; i<m_.size(); i++) {
            m_[i] = m[i]; //copying
        }
    }

    // e.g.: v = m.row_apply(&std::valarray<int>::sum); however, has 
    // slightly worse performance than directly calling apply_sum(1) 
    template<class T> inline std::valarray<T> Matrix<T>::row_apply(
            int (std::valarray<T>::*f)() const) const
    {
        if (this->empty()) {
            throw std::out_of_range("Operation not defined for empty matrix.");
        }

        std::valarray<T> v(nrow());
        for (size_t i=0; i<nrow(); ++i) {
            v[i] = (m_[i].*f)();
        }
        return v;
    }

    template<class T> inline std::valarray<T> Matrix<T>::apply_sum(size_t dim) const
    {
        std::valarray<T> v;
        switch (dim) {
            case 1 :
                v.resize(nrow());
                for (size_t i=0; i<nrow(); ++i) {
                    v[i] = m_[i].sum();
                }
                break;
            case 2 :
                v.resize(ncol());
                for (size_t i=0; i<nrow(); ++i) {
                    v += m_[i];
                }
        }
        return v;
    }

    template<class T> inline std::valarray<T> Matrix<T>::apply_sum(
            size_t dim,
            const std::valarray<size_t>& indices) const
    {
        std::valarray<T> v;
        switch (dim) {
            case 1 : //sum within rows
                v.resize(nrow());
                for (size_t i=0; i<nrow(); ++i) {
                    v[i] = m_[i][indices].sum();
                }
                break;
            case 2 : //sum within columns
                v.resize(ncol());
                for (size_t i=0; i<indices.size(); ++i) {
                    v += m_[indices[i]];
                }
        }
        return v;
    }

    template<class T> inline std::valarray<T> Matrix<T>::apply_max(size_t dim) const
    {
        std::valarray<T> v;
        switch (dim) {
            case 1 :
                v.resize(nrow());
                for (size_t i=0; i<nrow(); ++i) {
                    v[i] = m_[i].max();
                }
                break;
            case 2 :
                v.resize(ncol());
                std::valarray<bool> mask;
                for (size_t i=0; i<nrow(); ++i) {
                    mask = m_[i] > v;
                    v[mask] = m_[i][mask];
                }
        }
        return v;
    }

    template<class T> inline std::valarray<T> Matrix<T>::apply_max(
            size_t dim,
            const std::valarray<size_t>& indices) const
    {
        std::valarray<T> v;
        switch (dim) {
            case 1 : //sum within rows
                v.resize(nrow());
                for (size_t i=0; i<nrow(); ++i) {
                    v[i] = m_[i][indices].max();
                }
                break;
            case 2 : //sum within columns
                v.resize(ncol());
                std::valarray<bool> mask;
                for (size_t i=0; i<indices.size(); ++i) {
                    mask = m_[indices[i]] > v;
                    v[mask] = m_[indices[i]][mask];
                }
        }
        return v;
    }
    template<class T> inline std::valarray<T> Matrix<T>::apply_min(size_t dim) const
    {
        std::valarray<T> v;
        switch (dim) {
            case 1 :
                v.resize(nrow());
                for (size_t i=0; i<nrow(); ++i) {
                    v[i] = m_[i].min();
                }
                break;
            case 2 :
                v.resize(ncol());
                std::valarray<bool> mask;
                for (size_t i=0; i<nrow(); ++i) {
                    mask = m_[i] < v;
                    v[mask] = m_[i][mask];
                }
        }
        return v;
    }

    template<class T> inline std::valarray<T> Matrix<T>::apply_min(
            size_t dim,
            const std::valarray<size_t>& indices) const
    {
        std::valarray<T> v;
        switch (dim) {
            case 1 : //sum within rows
                v.resize(nrow());
                for (size_t i=0; i<nrow(); ++i) {
                    v[i] = m_[i][indices].min();
                }
                break;
            case 2 : //sum within columns
                v.resize(ncol());
                std::valarray<bool> mask;
                for (size_t i=0; i<indices.size(); ++i) {
                    mask = m_[indices[i]] < v;
                    v[mask] = m_[indices[i]][mask];
                }
        }
        return v;
    }

    template<class T> inline Matrix<T>& Matrix<T>::operator=(const T& val)
    {
        for (size_t i=0; i<m_.size(); ++i) {
            m_[i] = val;
        }
        return *this;
    }
    template<class T> inline Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
    {
        //assert(nrow() == m.nrow() && ncol() == m.ncol());
        for (size_t i=0; i<nrow(); ++i) {
            m_[i] = m[i];
        }
        return *this;
    }
    template<class T> inline Matrix<T>& Matrix<T>::operator=(
            const std::valarray<T>& v)
    {
        //assert(ncol() == v.size());
        for (size_t i=0; i<nrow(); ++i) {
            m_[i] = v;
        }
        return *this;
    }

    template<class T> inline Matrix<T>& Matrix<T>::resize(size_t d1, size_t d2)
    {
        m_.resize(d1);
        for (size_t i=0; i<nrow(); i++) {
            m_[i].resize(d2);
        }
        return *this;
    }
} //namespace detail
} //namespace Permory

#endif

