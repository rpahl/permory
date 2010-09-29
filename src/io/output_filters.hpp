// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_output_filters_hpp
#define permory_io_outut_filters_hpp

#include <string.h>

#include <boost/iostreams/char_traits.hpp> //EOF, WOULD_BLOCK
#include <boost/iostreams/concepts.hpp>    //multichar_input_filter
#include <boost/iostreams/operations.hpp>  //get

#include "detail/config.hpp"

namespace Permory { namespace io {
    namespace bio = boost::iostreams;

    template<class charT> class Print_op {
        public:
            Print_op(std::ostreambuf_iterator<charT>& out, charT l, charT r) 
                : out_(out), left_(l), right_(r), cnt_(0) {}
            void operator()(charT x) {
                //if (cnt_++ > 4) {
                *out_++ = x;
                //}
                //else { *out_++ = '*'; }
                *out_++ = right_;
            }
        private:
            int cnt_;
            charT left_;
            charT right_;
            std::ostreambuf_iterator<charT>& out_; 
    };

    template<class T> class Outstream_op {
        public:
            Outstream_op(std::ostream_iterator<T>& out, T l, T r) 
                : out_(out), left_(l), right_(r), cnt_(0) {}
            void operator()(T x) {
                if (x != 2) {
                    *out_++ = x;
                    *out_++ = right_;
                }
            }
            //template <class T> void operator()(T const& x) { 
            //out_ << left_ << x << right_; }
        private:
            int cnt_;
            T left_;
            T right_;
            std::ostream_iterator<T>& out_; 
    };

    class Delim_output_filter : public bio::multichar_output_filter {
        public:
            explicit Delim_output_filter(char delim = ' ')
                : delim_(delim) {}

            template<typename Sink> std::streamsize write(
                    Sink& dest, const char* s, std::streamsize n)
            {
                std::streamsize i;
                for (i = 0; i < n; ++i) {
                    int c = s[i];
                    if (bio::put(dest, c)) {
                        if (c != '\n') { //no delimiter after last char
                            if (!bio::put(dest, delim_))
                                break;
                        }
                    }
                    else
                        break;
                }
                return i; 
            }

            template<typename Source> void close(Source&) { 
            } 
        private:
            char delim_;
    };
} // namespace io
} // namespace Permory

#endif

