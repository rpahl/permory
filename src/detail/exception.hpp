// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_exception_hpp
#define permory_detail_exception_hpp

#include <exception>
#include <stdexcept>
#include "config.hpp"

namespace Permory { namespace detail {

    class File_exception : public std::exception {
        public:
            explicit File_exception(const char* s) : s_(s) {}
            virtual const char* what() const throw() { return s_; }
        private:
            const char* s_;
    };

    class Dimension_error : public std::exception { 
        public:
            explicit Dimension_error(const char* s) : s_(s) {}
            virtual const char* what() const throw() { return s_; }
        private: 
            const char* s_;
    };

    class Math_error {
    };

} //namespace detail
} //namespace Permory

#endif
