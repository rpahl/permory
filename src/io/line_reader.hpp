// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_line_reader_hpp
#define permory_io_line_reader_hpp

//#include <cstdlib>
#include <string>
#include <vector>
#include <zlib.h>

#include <boost/filesystem.hpp>
//#include <boost/filesystem/operations.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "config.hpp"
#include "exception.hpp"
#include "io/file_info.hpp"
#include "io/input_filters.hpp"

namespace Permory 
{
    //const size_t BUFFSIZE = 4096; 
    const size_t BUFFSIZE = 65536; 
    namespace bio = boost::iostreams;
    namespace bfs = boost::filesystem;

    template<class T> class Liner {
        public:
            // Iterator pass through
            typedef typename std::vector<T>::const_iterator const_iterator;
            const_iterator begin() const { return buf_.begin(); }
            const_iterator end() const { return buf_.end(); }
            
            //  Ctor and Dtor
            Liner(const char* fn, bool x = false);
            ~Liner() { in_.reset(); } //close all devices

        private:
            bfs::path p_;
            std::vector<T> buf_;     //line buffer
            size_t lineCount_;          
            size_t nchar_;
            size_t ncharTotal_;          

            const char* fn_;            //source file name
            bio::filtering_istream in_;  //chain of filters and input device
    };

    template<> class Liner<char> {
        public:
            // Iterator pass through
            typedef std::vector<char>::const_iterator const_iterator;
            const_iterator begin() const { return buf_.begin(); }
            const_iterator end() const { return buf_.end(); }

            //  Ctor and Dtor
            Liner(File_handler&); 
            ~Liner() { in_.reset(); } //close all devices

            // Modifier
            void next();            //jump to next line 
            void nextn(size_t n);   //jump n lines forward

            // Inspector
            bool empty() const { return nchar_ == 0; }
            bool eof() const { return it_ == eos_; }
            size_t char_count() const { return ncharTotal_; }
            size_t line_count() const { return lineCount_; }

        private:
            File_handler f_;
            std::vector<char> buf_; //line buffer
            size_t lineCount_;          
            size_t nchar_;
            size_t ncharTotal_;          

            std::string ext_;   //file extension
            bio::filtering_istreambuf in_; //chain of filters and input device
            std::istreambuf_iterator<char> it_;
            std::istreambuf_iterator<char> eos_;
    };

    /*
    //TODO move to reading function
    if (ext_ =  "slide") { //SLIDE format
    }
    if (ext_ =  "012") { //Permory simple genotype format
//no filter needed
}
*/

inline Liner<char>::Liner(File_handler& f) : f_(f) 
{
    buf_.resize(BUFFSIZE);
    ext_ =  (*f_).extension(); // file extension
    if (ext_ ==  "z") 
        in_.push(bio::zlib_decompressor()); 
    if (ext_ ==  "gz")
        in_.push(bio::gzip_decompressor()); 
    bool isCompressed = (in_.size() > 0);
    if (isCompressed) {
        ext_ = f_.extension(1); //ignore compression extension
        if (ext_.empty())
            throw File_exception("not a regular file name.");
    }

    // the input streambuffer
    in_.push(bio::file_source((*f_).string()));
    assert (in_.is_complete());
    bool isOpen = 
        in_.component<bio::file_source>(in_.size()-1)->is_open();
    if (!isOpen) 
        throw File_exception("unable to open file.");

    it_ = &in_;
    std::istreambuf_iterator<char> eos;
    eos_ = eos;
}
/*
   Move to where a Liner object is initialized, i.e. 
   try {
   Liner<char>("my_file")
   }

   catch (File_exception& e) {
   std::cerr << fn << ": " << e.what() << std::endl; 
   exit(-1);
   }
   catch (std::exception& e) {
   std::cerr << e.what() << std::endl;
   exit(-1);
   }
   */

inline void Liner<char>::next() 
{
    nchar_ = 0;
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
    ncharTotal_+= nchar_;
    lineCount_ += (nchar_ > 0);
    if (buf_.back() == '\r')    //windows line ending artefact
        buf_.pop_back();    
    if (it_ != eos_) { 
        it_++; //point to next line
    }
}

inline void Liner<char>::nextn(size_t n) 
{
    while (n-- > 0)
        this->next();
}

class Line_reader {
    public:
        // Iterator pass through
        typedef std::vector<char>::const_iterator const_iterator;
        const_iterator begin() const { return buf_.begin(); }
        const_iterator end() const { return buf_.begin() + nchar_; }

        //  Ctor and Dtor
        Line_reader(const char* fn, bool x = false);
        ~Line_reader() { in_.reset(); } //close all devices

        // Modifier
        Line_reader& operator++();      //next line (prefix operator)
        Line_reader& operator+=(int n); //skip n lines 

        // Inspector
        const char* operator*() const { return &buf_[0]; }
        const char* file_name() const { return fn_; }
        bool empty() const { return nchar_ == 0; }
        bool eof() const { return in_.eof(); }
        bool good() const { return in_.good(); }
        size_t line_count() const { return lineCount_; }
        size_t char_count() const { return ncharTotal_; }
        size_t get_file_size() const { return fileSize_; }

    private:
        void init() throw (File_exception); 
        std::vector<char> buf_;     //line buffer
        size_t lineCount_;          
        size_t nchar_;
        size_t ncharTotal_;          
        size_t fileSize_;

        const char* fn_;            //source file name
        bio::filtering_istream in_;  //chain of filters and input device
        bool extract_;
};


// ========================================================================
// Alternative line reader (not used earlier implementation)
// ========================================================================
class Line_reader_alt {
    public:
        // Iterator pass through
        typedef std::vector<char>::const_iterator const_iterator;
        const_iterator begin() const { return buf_.begin(); }
        const_iterator end() const { return buf_.begin() + nchar_; }

        //  Ctor and Dtor
        Line_reader_alt(const char* fn, bool x = false);
        ~Line_reader_alt() { in_.reset(); } //close all devices

        // Modifier
        Line_reader_alt& operator++();      //next line (prefix operator)
        Line_reader_alt& operator+=(int n); //skip n lines 

        // Inspector
        //std::vector<char> operator*() const;
        const char* operator*() const { return &buf_[0]; }
        const char* file_name() const { return fn_; }
        bool empty() const { return nchar_ == 0; }
        bool eof() const { return it_ == eos_; }
        size_t char_count() const { return ncharTotal_; }
        size_t file_size() const { return fileSize_; }
        size_t line_count() const { return lineCount_; }

    private:
        void init() throw (File_exception); 
        std::vector<char> buf_;             //line buffer
        size_t lineCount_;          
        size_t nchar_;
        size_t ncharTotal_;          
        size_t fileSize_;

        const char* fn_;                    //source file name
        bio::filtering_istreambuf in_;      //chain of filters and input device
        std::istreambuf_iterator<char> it_;
        std::istreambuf_iterator<char> eos_;
        bool extract_;
};

} // namespace Permory

#endif
