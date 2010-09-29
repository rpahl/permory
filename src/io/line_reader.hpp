// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_line_reader_hpp
#define permory_io_line_reader_hpp

#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>

#include "detail/config.hpp"
#include "detail/exception.hpp"
#include "io/file.hpp"
#include "io/input_filters.hpp"

namespace Permory { namespace io {
    namespace bio = boost::iostreams;
    namespace bfs = boost::filesystem;
    using namespace Permory::detail; //exceptions

    template<class T> class Line_reader {
        public:
            typedef T char_t;

            // Iterator pass through
            typedef typename std::vector<T>::const_iterator const_iterator;
            const_iterator begin() const { return buf_.begin(); }
            const_iterator end() const { return buf_.end(); }

            //  Ctor and Dtor
            Line_reader(const std::string& fn); 
            ~Line_reader() { in_.reset(); } //close all devices

            // Modification
            void next();        //read next line (formatted) into this->buf_
            void next(std::string&); //as above plus get the unformatted string
            void next_str(std::string& s); //gets just the unformatted string
            void skip();        //skip one line

            // Inspection
            bool empty() const { return buf_.empty(); }
            bool eof() const { return in_.eof(); }
            size_t size() const { return buf_.size(); }
            size_t char_count() const { return charCount_; }
            size_t line_count() const { return lineCount_; }
            size_t word_count() const { return wordCount_; }
            const File_handle get_file() const { return file_; }

        private:
            File_handle file_;
            std::vector<T> buf_;     //line buffer
            size_t charCount_;
            size_t lineCount_;          
            size_t wordCount_;          
            std::string ext_;   //the file name extension
            bio::filtering_istream in_;  //chain of filters and input device
    };

    // Line_reader implementation
    // ========================================================================
    template<class T> inline Line_reader<T>::Line_reader(const std::string& fn) 
        : file_(fn), charCount_(0), lineCount_(0), wordCount_(0) 
    {
        if (not bfs::is_regular(*file_)) {
            throw File_exception("not a regular file.");
        }
        if (not bfs::exists(*file_)) { 
            throw File_not_found("not found.");
        }

        ext_ =  (*file_).extension(); 
        if (ext_ ==  ".gz") {
            in_.push(bio::gzip_decompressor()); 
            ext_ = file_.extension(1); //discard compression extension
            if (ext_.empty())
                throw File_exception("not a regular file name.");
        }

        // the input streambuffer
        in_.push(bio::file_source((*file_).string()));
        assert (in_.is_complete());
        bool isOpen = 
            in_.component<bio::file_source>(in_.size()-1)->is_open();
        if (!isOpen) 
            throw File_exception("unable to open file.");
    }

    template<class T> inline void Line_reader<T>::next() 
    {
        std::string s;
        this->next(s);
    }

    template<class T> inline void Line_reader<T>::next(std::string& s) 
    {
        getline (in_, s);
        charCount_ += s.size();

        std::stringstream ss(s);
        std::istream_iterator<T> start(ss);
        std::istream_iterator<T> end;
        buf_.assign(start, end);

        wordCount_ += buf_.size();
        lineCount_ += (!buf_.empty());
    }

    template<class T> inline void Line_reader<T>::next_str(std::string& s)
    {
        getline (in_, s);
        charCount_ += s.size();
        lineCount_ += (!s.empty());
    }

    template<class T> inline void Line_reader<T>::skip()
    {
        std::string s;
        getline (in_, s);
    }


    // Specialization to char for increased performance
    // ========================================================================
    const size_t BUFFSIZE = 65536; 
    template<> class Line_reader<char> {
        public:
            typedef char char_t;

            // Iterator pass through
            typedef std::vector<char>::const_iterator const_iterator;
            const_iterator begin() const { return buf_.begin(); }
            const_iterator end() const { return buf_.end(); }

            //  Ctor and Dtor
            Line_reader(const std::string& fn); 
            ~Line_reader() { in_.reset(); } //close all devices

            // Modification
            void next();    //read next line 
            void skip();    //skip one line

            // Inspection
            bool empty() const { return buf_.empty(); }
            bool eof() const { return it_ == eos_; }
            size_t size() const { return buf_.size(); }
            size_t char_count() const { return charCount_; }
            size_t line_count() const { return lineCount_; }
            const File_handle get_file() const { return file_; }

        private:
            File_handle file_;
            std::vector<char> buf_; //line buffer
            size_t lineCount_;          
            size_t charCount_;          

            std::string ext_;   //the file name extension
            bio::filtering_istreambuf in_; //chain of filters and input device
            std::istreambuf_iterator<char> it_;
            std::istreambuf_iterator<char> eos_;
    };


    inline Line_reader<char>::Line_reader(const std::string& fn) 
        : file_(fn), charCount_(0), lineCount_(0)
    {
        buf_.resize(BUFFSIZE);
        ext_ =  (*file_).extension(); // file extension
        if (ext_ ==  ".gz") {
            in_.push(bio::gzip_decompressor()); 
            ext_ = file_.extension(1); //discard compression extension
            if (ext_.empty())
                throw File_exception("not a regular file name.");
        }

        // the input streambuffer
        in_.push(bio::file_source((*file_).string()));
        assert (in_.is_complete());
        bool isOpen = 
            in_.component<bio::file_source>(in_.size()-1)->is_open();
        if (!isOpen) 
            throw File_exception("unable to open file.");

        it_ = &in_;
        std::istreambuf_iterator<char> eos;
        eos_ = eos;
    }

    inline void Line_reader<char>::next() 
    {
        size_t nchar_ = 0;
        std::vector<char>::iterator i = buf_.begin();
        while (it_ != eos_ && *it_ != '\n') {
            if (i == buf_.end()) {          //buffer full?
                size_t sz = buf_.size();    //remember where we are
                buf_.resize(buf_.size() + BUFFSIZE); //enlarge buffer
                i = buf_.begin() + sz;      //read on from where we were
            }
            *i++ = *it_++;
            nchar_++;
        }
        buf_.resize(nchar_);    //chop off (possible) characters from last line
        charCount_+= nchar_;
        lineCount_ += (nchar_ > 0);
#if defined(_WINDOWS) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)
        if (!buf_.empty()) {
            if (buf_.back() == '\r') //remove possible windows line ending artefact
                buf_.pop_back();    
        }
#endif
        if (it_ != eos_) { 
            it_++; //point to next line
        }
    }

    inline void Line_reader<char>::skip() 
    {
        char ch;
        while (it_ != eos_ && *it_ != '\n') 
            ch = *it_++;

        if (it_ != eos_) { 
            it_++; //point to next line
        }
    }

} // namespace io
} // namespace Permory

#endif
