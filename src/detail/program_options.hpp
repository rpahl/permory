// Copyright (c) 2011 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_program_options_hpp
#define permory_detail_program_options_hpp

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "detail/config.hpp"
#include "detail/exception.hpp"

namespace Permory { namespace detail {

    //
    //  Extends boost program_options library by user-chosen names for the
    //  arguments types like INT, REAL, STRING ... With the original library 
    //  one can only specify a single global name defaulting to "arg".
    //
    template<class T> 
        class My_typed_value : public boost::program_options::typed_value<T>
    {
        public:
            My_typed_value(T* val, std::string name="ARG")
                : boost::program_options::typed_value<T>(val), name_(name)
            { }
            My_typed_value(std::string name="ARG")
                : boost::program_options::typed_value<T>(0), name_(name)
            { }

            // Inspection
            std::string name() const { return name_; }

            // Modification
            My_typed_value* my_default_value(const T& val)
            {
                name_.append(" (=");
                std::ostringstream oss; 
                oss << val;
                name_.append(oss.str());
                name_.append(")");
                this->default_value(val);
                return this;
            }
            My_typed_value* my_default_value(const T& val, const std::string& textual)
            {
                name_.append(" ");
                name_.append(textual);
                this->default_value(val);
                return this;
            }
        private:
            std::string name_;
    };

    //
    //  Functions providing shorter syntax similar to program_options library
    template<class T>
        My_typed_value<T>* my_value(T* val, std::string name="ARG") {
            My_typed_value<T>* r = new My_typed_value<T>(val, name);
            return r;        
        }

    template<class T>
        My_typed_value<T>* my_value(std::string name="ARG") {
            My_typed_value<T>* r = new My_typed_value<T>(name);
            return r;        
        }


} //namespace detail
} //namespace Permory

#endif


