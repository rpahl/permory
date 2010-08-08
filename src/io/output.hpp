// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_output_hpp
#define permory_io_output_hpp

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <zlib.h>

#include <boost/iostreams/copy.hpp> 
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/tee.hpp>

#include "config.hpp"
#include "detail/exception.hpp"
#include "io/file.hpp"
#include "io/output_filters.hpp"

namespace Permory { namespace io {
    namespace bio = boost::iostreams;
    namespace bfs = boost::filesystem;

    template<class T, class U=std::vector<T> > class Line_writer {
        public:
            typedef T char_t;
            typedef typename U::const_iterator iter;
            
            //  Ctor and Dtor
            Line_writer(const File_handle& f, bool a=false); 
            ~Line_writer() { out_.reset(); } //close all devices

            // Modification
            void insert(const char* s) { out_ << s; }   //no line feed
            void next() { out_ << std::endl; }          //empty line
            // straight output of s to out_ 
            void next(const char* s) { out_ << s << std::endl; }
            // write to out using delimiter after each word
            void next(iter start, iter end, const char* delim);

            // Inspection
            bool file_exists() const { return bfs::exists(*file_); }
            size_t line_count() const { return lineCount_; }
            size_t word_count() const { return wordCount_; }
        private:
            File_handle file_;
            size_t lineCount_;          
            size_t wordCount_;          
            bool useCout_;
            std::string ext_;   //the file name extension
            bio::filtering_ostream out_; //chain of filters and output device
    };

    template<class T, class U> inline Line_writer<T, U>::Line_writer(
            const File_handle& f, bool a)
        : file_(f), useCout_(a), lineCount_(0), wordCount_(0) 
    {
        bio::file_sink fs((*file_).filename());
        bool isOpen = fs.is_open();
        if (!isOpen)
            throw File_exception("failed to open file.");

        if (useCout_) { 
            // fork stream to both file stream and std::cout
            out_.push(bio::tee(fs));
            out_.push(std::cout);
        }
        else {
            if ((*file_).extension() ==  ".gz") {
                out_.push(bio::gzip_compressor()); 
            }
            out_.push(fs);
        }
        assert (out_.is_complete());
    }

    template<class T, class U> inline void Line_writer<T, U>::next(
            iter start, iter end, const char* delim)
    {
        std::ostream_iterator<T> to_out(out_, delim);
        std::copy(start, end, to_out);
        //std::for_each(start, end, Outstream_op<int>(out_, ' ', ' '));
        out_ << std::endl << std::flush;
        wordCount_ += distance(start, end);
        lineCount_++;
    }


    // Specialization to char for increased performance
    // ========================================================================
    template<> class Line_writer<char> {
        public:
            typedef char char_t;
            typedef std::vector<char>::iterator iter;
            
            //  Ctor and Dtor
            Line_writer(const File_handle& f, bool a=false); 
            ~Line_writer() { outbuf_.reset(); } //close all devices

            // Modification
            void next();          //empty line
            void next(iter start, iter end);

            // Inspection
            bool file_exists() const { return bfs::exists(*file_); }
            size_t char_count() const { return charCount_; }
            size_t line_count() const { return lineCount_; }
        private:
            File_handle file_;
            size_t charCount_;          
            size_t lineCount_;          
            bool useCout_;
            std::string ext_;   //the file name extension
            bio::filtering_streambuf<bio::output> outbuf_; //output device
    };


    inline Line_writer<char>::Line_writer(const File_handle& f, bool a)
        : file_(f), useCout_(a), charCount_(0), lineCount_(0)
    {
        bio::file_sink fs((*file_).filename());
        bool isOpen = fs.is_open();
        if (!isOpen)
            throw File_exception("failed to open file.");

        if ((*file_).extension() ==  ".gz") {
            outbuf_.push(bio::gzip_compressor()); 
        }
        outbuf_.push(fs);
        assert (outbuf_.is_complete());
    }

    inline void Line_writer<char>::next()
    {
        std::vector<char> v; //dummy string
        this->next(v.begin(), v.end()); //will just print an empty line
    }
    inline void Line_writer<char>::next(iter start, iter end)
    {
        if (useCout_) {
            std::ostreambuf_iterator<char> to_cout(std::cout);
            std::copy(start, end, to_cout);
            std::cout << std::endl << std::flush;
        }
        std::ostreambuf_iterator<char> out(&outbuf_);
        std::copy(start, end, out);
        //std::for_each(start, end, Print_op<char>(out, ' ', ' '));
#if defined(_WINDOWS) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)
        *out++ = '\r'; 
#endif
        *out++ = '\n';
        charCount_ += distance(start, end);
        lineCount_++;
        bio::flush(outbuf_); //mandatory to prevent Dtor from discarding chars
    }

} // namespace io
} // namespace Permory

#endif
