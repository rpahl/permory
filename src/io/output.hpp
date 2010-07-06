/**
 * @author Roman Pahl
 */

#ifndef permory_io_output_hpp
#define permory_io_output_hpp

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/tee.hpp>

#include "config.hpp"
#include "exception.hpp"
#include "io/file_info.hpp"
#include "io/output_filters.hpp"

namespace Permory 
{
    namespace bio = boost::iostreams;

    class Output {
        public:
            typedef std::vector<int>::const_iterator iter;
            // Ctor and Dtor
            Output(const char* fn="", 
                    bool useCout=false,
                    bool z=false, 
                    bool i=false);
            ~Output() { out_.reset(); } //close all devices
            void write_form(const char*);
            void write_stream(iter, iter);
            void write_unform(iter, iter);
        private:
            void init() throw (File_exception); 
            size_t lineCount_;
            size_t charCount_;

            const char* fn_;                //sink file name
            bio::filtering_ostream out_;    //chain of filters and output device
            bio::filtering_streambuf<bio::output> outbuf_;    //chain of filters and output device
            bool useCout_;
            bool zip_;                      //use compression
            bool interactive_;              //if true, prompt before overwrite
    };

    class file_out : public Output {
        const char* fn_;            //source file name
    };


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
} // namespace Permory

//TODO two classes: (1) formatted output using ostream_iterator<T> and (2)
//unformatted faster output using ostreambuf_iterator<char> or both in one or
//via specialization for char
//
/*
   filebuf fb;
   fb.open ("test.txt",ios::in);
   istream is(&fb);
   istream_iterator<std::string> eos;         
   istream_iterator<std::string> iit (is);   
   std::vector<std::string> v;
   while (iit != eos) {
   v.push_back(*iit++);
   }
   PRINT(v.size());
   print_vec(v, " ");
//ostream os(&fb);
//ostream_iterator<double> out_it (os,", ");
//ostreambuf_iterator<char> out_it (cout);
//copy ( v.begin(), v.end(), out_it );
//return 0;
*/

#endif
