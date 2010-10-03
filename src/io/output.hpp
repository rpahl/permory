// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_output_hpp
#define permory_io_output_hpp

#include <boost/scoped_ptr.hpp>

#include "detail/parameter.hpp"
#include "detail/enums.hpp"
#include "file_out.hpp"

namespace Permory { namespace io {

    // 
    // Define own manipulators
    //
    std::ostream& stdpre(std::ostream& os) { //standard prefix
        return os << "..."; 
    }
    std::ostream& errpre(std::ostream& os) { //error prefix
        return os << "!!!"; 
    }
    std::ostream& warnpre(std::ostream& os) { //warning prefix
        return os << "!"; 
    }
    class indent { 
        public: 
            indent(size_t n, char c=' ') : indent_string_(n, c) {}
            std::ostream& operator()(std::ostream& os) const {
                return os << indent_string_;
            }
            io::File_out& operator()(io::File_out& fo) const {
                return fo << indent_string_;
            }
        private:
            const std::string indent_string_;
    };

    //
    // Output class to allow/suppress output to console (and/or file) 
    // for different verbose levels
    //
    class Myout : std::ostream {
        public:
            // Ctor
            Myout(detail::Parameter* par);

            // support standard I/O manipulators like std::endl, std::flush, ...
            Myout& operator<<(std::ostream& (*f)(std::ostream&)); 
            template<class T> Myout& operator<<(const T& x);
            template<class T> void operator()(const T& x);
            Myout& all() { verb_ = detail::all; return *this; }
            Myout& verbose() { verb_ = detail::verbose; return *this; }
            Myout& normal() { verb_ = detail::normal; return *this; }
            Myout& mute() { verb_ = detail::muted; return *this; }

            // Modification
            void set_verbosity(detail::Verbosity v) { verb_ref_ = v; }
            void set_logfile(std::string, bool); 

        private:
            detail::Verbosity verb_ref_;
            detail::Verbosity verb_;
            boost::scoped_ptr<File_out> out_;
    };

    // Myout implementation
    // ========================================================================
    Myout::Myout(detail::Parameter* par)
        : verb_ref_(par->verbose_level), verb_(verb_ref_)
    {
        bool useLogfile = not par->log_file.empty();
        if (useLogfile) {
            this->set_logfile(par->log_file, par->interactive);
        }
    }

    inline void Myout::set_logfile(std::string fn, bool interactive)
    {
        File_handle file(fn);
        if (bfs::exists(*file) && interactive) {
            char c;
            std::cerr << "Logfile: overwrite `" << fn << "'? (y/n) ";
            std::cin >> c;
            if (c != 'y') {
                return;
            }
        }
        out_.reset(new File_out(fn));
    }

    template<> inline Myout& Myout::operator<<(const detail::Verbosity& v)  
    {
        this->verb_ = v;
        return *this;
    } 

    template<> inline Myout& Myout::operator<<(const indent& x)  
    {
        using namespace detail;
        if (verb_ >= verb_ref_ && verb_ref_ != muted) {
            if (out_.get() != 0) {
                x(*out_);
            }
            x(std::cout);
        }
        return *this;
    } 

    inline Myout& Myout::operator<<(std::ostream& (*f)(std::ostream&))
    {
        using namespace detail;
        if (verb_ >= verb_ref_ && verb_ref_ != muted) {
            if (out_.get() != 0) {
                *out_ << f;
            }
            std::cout << f;
        }
        return *this;
    }

    template<class T> inline Myout& Myout::operator<<(const T& x)  
    {
        using namespace detail;
        if (verb_ >= verb_ref_ && verb_ref_ != muted) {
            if (out_.get() != 0) {
                (*out_) << x;
            }
            std::cout << x;
        }
        return *this;
    } 

    template<class T> inline void Myout::operator()(const T& x) 
    {
        using namespace detail;
        if (verb_ >= verb_ref_ && verb_ref_ != muted) {
            if (out_.get() != 0) {
                (*out_) << x;
            }
            std::cout << x;
        }
    } 

} // namespace io
} // namespace Permory

#endif

