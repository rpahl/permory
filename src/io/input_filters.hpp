// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_input_filters_hpp
#define permory_io_input_filters_hpp

#include <string.h>

#include <boost/iostreams/char_traits.hpp> //EOF, WOULD_BLOCK
#include <boost/iostreams/concepts.hpp>    //multichar_input_filter
#include <boost/iostreams/operations.hpp>  //get

#include "config.hpp"

namespace Permory { namespace io {
    namespace bio = boost::iostreams;

    //
    // Skips for three special delimiters, namely "empty char", tabulator 
    // and '\r'. In addition, the first n chunks ("defined" by the delimiters) 
    // are also skipped.
    //
    template<int n> class Skip_input_filter {
        public:
            Skip_input_filter() 
                : cnt_(0), skip_(false) 
            {}
            bool operator()(char c) {
                if (c == ' ' || c == '\t' || c == '\r') { //chars to skip
                    if (!skip_) { //catch more than one "skip-char" in a row
                        cnt_++;
                    }
                    skip_ = true;
                }
                else {
                    skip_ = false;
                }
                return (skip_ || cnt_ < n);
            }
        private:
            int cnt_;
            bool skip_;
    };

    //
    // Filters out each occurcence of some prespecified delimiter
    //
    class Delim_input_filter : public bio::multichar_input_filter {
        public:
            explicit Delim_input_filter(char delim = ' ')
                : delim_(delim) {}

            template<class Source> std::streamsize read(
                    Source& src, char* s, std::streamsize n)
            {
                int c;
                for (std::streamsize i=0; i<n; ++i) {
                    do {
                        c = bio::get(src);
                        if (c == EOF)
                            return i != 0 ? i : -1;
                        else if (c == bio::WOULD_BLOCK)
                            return i; //read i chars
                    }
                    while (c == delim_);
                    s[i] = c;
                }
                return n; //read n chars
            }
        private:
            char delim_;
    };
    

    //
    // Allows arbitrary string of characters as delimiters, and therefore is
    // slower than the Delim_input_filter
    //
    class String_delim_input_filter : public bio::multichar_input_filter {
        public:
            explicit String_delim_input_filter(const char* delims=" \t\r")
                : delims_(delims) {}

            template<class Source> std::streamsize read(
                    Source& src, char* s, std::streamsize n)
            {
                int c;
                for (std::streamsize i=0; i<n; ++i) {
                    do {
                        c = bio::get(src);
                        if (c == EOF)
                            return i != 0 ? i : -1;
                        else if (c == bio::WOULD_BLOCK)
                            return i;
                    }
                    while (strchr(delims_, c));
                    s[i] = c;
                }
                return n;
            }
        private:
            const char* delims_;
    };


    //
    // Same as the Delim_input_filter plus skip the first 'nskip' chunks, where
    // a chunk is one (or more) character(s) surrounded by delimiters
    //
    class Delim_and_skip_input_filter : public bio::multichar_input_filter {
        public:
            explicit Delim_and_skip_input_filter(char delim = ' ', int nskip=0)
                : delim_(delim), nskip_(nskip), chunkCount_(0), isChunkComplete_(true)
            { }

            template<class Source> std::streamsize read(
                    Source& src, char* s, std::streamsize n)
            {
                char c;
                for (std::streamsize i=0; i<n; ++i) {
                    do {
                        c = bio::get(src);
                        if (c == EOF)
                            return i != 0 ? i : -1;
                        else if (c == bio::WOULD_BLOCK)
                            return i;

                        if (c == '\n') {
                            chunkCount_ = 0;//reset
                            isChunkComplete_ = true;
                            break; //always break new lines
                        }
                        if (c == delim_) {
                            if(!isChunkComplete_) { //catch series of delimiters 
                                isChunkComplete_ = true;
                                chunkCount_++; 
                            }
                        }
                        else { 
                            isChunkComplete_ = false; 
                        }
                    }
                    while (c == delim_ || nskip_ > chunkCount_);
                    s[i] = c;
                }
                return n;
            }

            // Make filter re-usable in case stream is closed
            template<typename Source> void close(Source&) { 
                chunkCount_ = 0;
                isChunkComplete_ = true;
            }
        private:
            char delim_;
            int nskip_;
            int chunkCount_;
            bool isChunkComplete_;
    };
} // namespace io
} // namespace Permory

#endif

/* Boost iostreams filter example
 * const char* fn;
 * char line[4096];
 * in.push(bio::counter()); //count characters
 * in.push(bio::file_source(fn));
 * while (in.getline(line, 4096)) { }
 * int lines = in.component<0, bio::counter>()->lines();
 * int characters = in.component<0, counter>()->characters();
 * bool isOpen = in.component<bio::file_source>(0)->is_open();
 */
