// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_exception_hpp
#define permory_detail_exception_hpp

#include <exception>
#include <stdexcept>

#include <boost/lexical_cast.hpp>

#include "detail/config.hpp"

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

    //
    // Error for mismatching length of phenotype and marker data.
    // Thrown in Analyzer::analyze_dichotom().
    class Data_length_mismatch_error : public std::runtime_error
    {
        public:
            // Ctor
            Data_length_mismatch_error(size_t id, size_t pheno_length, size_t marker_length)
                : runtime_error(
                            std::string("")
                            .append( "At marker no ")
                            .append(boost::lexical_cast<std::string>(id))
                            .append(": length of phenotype data (")
                            .append(boost::lexical_cast<std::string>(pheno_length))
                            .append(") does not match with length of marker data (")
                            .append(boost::lexical_cast<std::string>(marker_length))
                            .append(").\n")
                            .append(
                                (marker_length == 2 * pheno_length
                                    ? "Maybe you forgot to set option '--allelic'?\n"
                                    : ""))
                        )
                { }
    };

    class Missing_option : public std::runtime_error {
        public:
            explicit Missing_option(const std::string& s) : std::runtime_error(s) {}
    };
    
    class Ambigous_option : public std::runtime_error {
        public:
            explicit Ambigous_option(const std::string& s) : std::runtime_error(s) {}
    };
} //namespace detail
} //namespace Permory

#endif
