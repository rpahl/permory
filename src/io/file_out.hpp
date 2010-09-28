// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_file_out_hpp
#define permory_io_file_out_hpp

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

    class File_out {
        public:
            //  Ctor and Dtor
            File_out(const std::string& fn);
            ~File_out() { out_.reset(); } //close all devices

            // Output
            template<class T> File_out& operator<<(const T& x);
            File_out& endl(); 
            File_out& flush(); 
            // write to File_out using delimiter after each word
            template<class T> File_out& operator()(
                    typename std::vector<T>::const_iterator start,
                    typename std::vector<T>::const_iterator end,
                    const std::string& delim);

            // Inspection
            bool file_exists() const { return bfs::exists(*file_); }
            size_t line_count() const { return lineCount_; }
            size_t word_count() const { return wordCount_; }

        private:
            File_handle file_;
            size_t lineCount_;          
            size_t wordCount_;          
            std::string ext_;   //the file name extension
            bio::filtering_ostream out_; //chain of filters and output device
    };

    // File_out implementation
    // ========================================================================
    inline File_out::File_out(const std::string& fn)
        : file_(fn), lineCount_(0), wordCount_(0) 
    {
        bio::file_sink fs((*file_).filename());
        bool isOpen = fs.is_open();
        if (!isOpen) {
            throw File_exception("failed to open file.");
        }
        /*
        if (useCout_) { 
            // fork stream to both file stream and std::cout
            out_.push(bio::tee(fs));
            out_.push(std::cout);
        }
        */

        if ((*file_).extension() ==  ".gz") {
            out_.push(bio::gzip_compressor()); 
        }
        out_.push(fs);
        assert (out_.is_complete());
    }

    template<class T> inline File_out& File_out::operator<<(const T& x) 
    {
        out_ << x; 
        wordCount_++;
        return *this; 
    }

    inline File_out& File_out::endl()
    {
        out_ << std::endl;
        lineCount_++;
        return *this;
    }

    inline File_out& File_out::flush()
    {
        out_ << std::flush;
        return *this;
    }

    template<class T> inline File_out& File_out::operator()(
            typename std::vector<T>::const_iterator start,
            typename std::vector<T>::const_iterator end,
            const std::string& delim)
    {
        std::ostream_iterator<T> to_out(out_, delim);
        std::copy(start, end, to_out);
        wordCount_ += distance(start, end);
        return *this;
    }

    // Version using char buffer with increased performance
    // ========================================================================
    class File_out_char {
        public:
            typedef std::vector<char>::const_iterator iter;

            //  Ctor and Dtor
            File_out_char(const std::string& fn); 
            ~File_out_char() { outbuf_.reset(); } //close all devices

            // Modification
            File_out_char& endl(); 
            File_out_char& operator()(iter start, iter end);

            // Inspection
            size_t char_count() const { return charCount_; }
            size_t line_count() const { return lineCount_; }

        private:
            File_handle file_;
            size_t charCount_;          
            size_t lineCount_;          
            std::string ext_;   //the file name extension
            bio::filtering_streambuf<bio::output> outbuf_; //output device
    };


    inline File_out_char::File_out_char(const std::string& fn)
        : file_(fn), charCount_(0), lineCount_(0)
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

    inline File_out_char& File_out_char::operator()(iter start, iter end)
    {
        std::ostreambuf_iterator<char> out(&outbuf_);
        std::copy(start, end, out);
        charCount_ += distance(start, end);
        bio::flush(outbuf_); //mandatory to prevent Destructor discarding chars
        /*
           if (useCout_) {
           std::ostreambuf_iterator<char> to_cout(std::cout);
           std::copy(start, end, to_cout);
           }
           */
        return *this;
    }

    inline File_out_char& File_out_char::endl()
    {
        std::ostreambuf_iterator<char> out(&outbuf_);
#if defined(_WINDOWS) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)
        *out++ = '\r'; 
#endif
        *out++ = '\n';
        bio::flush(outbuf_); //mandatory 
        lineCount_++;
        return *this;
    }
} // namespace io
} // namespace Permory

#endif
