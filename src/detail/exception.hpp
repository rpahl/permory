// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_exception_hpp
#define permory_detail_exception_hpp

#include <exception>
#include <stdexcept>
#include "config.hpp"

namespace Permory { namespace detail {

    class File_exception : public std::runtime_error {
        public:
            explicit File_exception(const std::string& s) 
                : std::runtime_error(s)
            {}
    };

    class File_not_found : public File_exception {
        public:
            explicit File_not_found(const std::string& s) 
                : File_exception(s)
            {}
    };

    class Dimension_error : public std::out_of_range { 
        public:
            explicit Dimension_error(const std::string& s) 
                : std::out_of_range(s)
            {}
    };

    class Math_error {
    };

} //namespace detail
} //namespace Permory

#endif
