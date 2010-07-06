/**
 * @author Roman Pahl
 */

#include "output.hpp"
using namespace std;
using namespace Permory;

Output::Output(const char* fn, bool useCout, bool z, bool i)
    : fn_(fn), useCout_(useCout), zip_(z), interactive_(i)
{
    try {
        init();
    }
    catch (File_exception& e) {
        cerr << fn_ << ": " << e.what() << std::endl; 
        exit(-1);
    }
    catch (exception& e) {
        cerr << e.what() << std::endl;
        exit(-1);
    }
}

void Output::init() throw (File_exception)
{
    // Create filter-device chain for the input stream
    //if (zip_) in_.push(bio::zlib_compressor()); //TODO
    //out_.push(Delim_output_filter(','));
    //outbuf_.push(Delim_output_filter(','));
    
    if (strlen(fn_) > 0) { //no file name means to use std::cout 
        File_info fi(fn_);
        if (interactive_ && fi.file_exists()) {
            cerr << "Overwrite existing file '" << fn_ << "'? (y/n): ";
            char c; cin >> c;
            if(c != 'y') 
                exit(-1);
        }
        bio::file_sink fs(fn_);
        bool isOpen = fs.is_open();
        if (!isOpen) 
            throw File_exception("failed to open file.");
        if (!useCout_) {
            out_.push(fs); //only file stream
            //outbuf_.push(fs); //only file stream FIXME
        }
        else { //both file stream and cout 
            out_.push(bio::tee(fs));
            out_.push(std::cout);
        }
    }
    else {
        // use only std::cout
        out_.push(std::cout);
    }
    //assert (out_.is_complete());
}

void Output::write_form(const char* s)
{
    out_ << s << std::endl;
    //out_.flush();
    //std::ostream_iterator<int> to_cout(std::cout, " ");
}

void Output::write_stream(iter start, iter end) 
{
    std::ostream_iterator<int> out(out_, " ");
    //std::copy(start, end, out);
    std::for_each(start, end, Outstream_op<int>(out, 0, 0));
    // if want to use out_ directly, i.e. //FIXME 
    //std::for_each(start, end, Outstream_op<int>(out_, ' ', ' '));
    //maybe just respecify Outstream_op
    
    out_ << std::endl;
}

void Output::write_unform(iter start, iter end) 
{
    //std::ostream_iterator<class T::value_type> out_it(out_, " ");
    //std::ostream_iterator<int> to_cout(std::cout, " ");
    //std::copy(start, end, to_cout);

    std::ostreambuf_iterator<char> out = &outbuf_;
    //std::copy(start, end, out);
    std::for_each(start, end, Print_op<char>(out, ' ', ' '));
#ifdef __WIN32__
    *out++ = '\r'; 
#endif
    *out++ = '\n';
    //outbuf_ << std::endl;
}

