// Copyright (c) 2010-2014 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_exception_hpp
#define permory_detail_exception_hpp

#include <exception>
#include <stdexcept>

#include <boost/lexical_cast.hpp>

#include "detail/config.hpp"
#include "detail/functors.hpp"  //comp_second

namespace Permory { namespace detail {

    class File_exception : public std::runtime_error {
        public:
            explicit File_exception(const std::string& s) 
                : std::runtime_error(s)
            {}
    };

    class Missing_option : public std::runtime_error {
        public:
            explicit Missing_option(const std::string& s) : std::runtime_error(s) {}
    };
    
    class Ambigous_option : public std::runtime_error {
        public:
            explicit Ambigous_option(const std::string& s) : std::runtime_error(s) {}
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
                            .append("At marker no ")
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

    class Wrong_missing_value_error : public std::runtime_error
    {
        public:
            // Ctor
            Wrong_missing_value_error(size_t id, char na_set, char na_real)
                : runtime_error(
                            std::string("\n")
                            .append("At marker no ")
                            .append(boost::lexical_cast<std::string>(id))
                            .append(": found more than 2 alleles.\nIf '")
                            .append(std::string(1, na_real))
                            .append("' represents missing values in your ")
                            .append("data set, then it differs from what ")
                            .append("Permory assumes ('")
                            .append(std::string(1, na_set))
                            .append("'), and you probably forgot to set ")
                            .append("option'--missing ")
                            .append(std::string(1, na_real))
                            .append("'.\n")
                        )
                { }
    };

} //namespace detail
} //namespace Permory

#endif
