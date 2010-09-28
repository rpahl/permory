// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_tokenizer_hpp
#define permory_detail_tokenizer_hpp

#include "config.hpp"
#include <string.h> //strlen, strtok
#include <stdlib.h> //realloc
#include <vector>

namespace Permory { namespace detail {
    // Tokenizer with the primary aim to be efficient
    class Tokenizer {
        public:
            // Ctor
            Tokenizer(const char* source, const char* delims = " \t\n\r")
                : curr_(0), delims_(delims) { this->assign(source); }
            ~Tokenizer(){}

            // Inspector
            bool empty() const { return curr_ == 0; }

            // Modifier
            void assign(const char* source)
            { 
                size_t len = strlen(source);
                if (len > buf_.size()) {             
                    buf_.resize(strlen(source) + 1); 
                }

                if (len > 0) {
                    strcpy(&buf_[0], source);
                    curr_ = strtok(&*buf_.begin(), delims_); //point to first token
                }
                else {
                    curr_ = 0;
                }
            }

            Tokenizer& operator++() //"++a" prefix
            {
                curr_ = strtok(0, delims_); 
                return *this;
            }

            const char* operator*() const { return curr_; }
            Tokenizer& next(const char* delims) //alternative delimiters
            {
                curr_ = strtok(0, delims); 
                return *this;
            }

        private:
            char* curr_;            //pointing at current token
            std::vector<char> buf_; //local buffer
            const char* delims_;
    };
} // namespace detail
} // namespace Permory

#endif

