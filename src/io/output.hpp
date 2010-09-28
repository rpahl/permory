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
    std::ostream& myendl(std::ostream& os) {
        return os << std::endl; 
    }
    std::ostream& myflush(std::ostream& os) {
        return os << std::flush; 
    }

    class Out_log {
        public:
            // Ctor
            Out_log(detail::Parameter* par);

            // Output
            template<class T> Out_log& operator<<(const T& x);
            template<class T> void operator()(const T& x);
            Out_log& all() { verb_ = detail::all; return *this; }
            Out_log& verbose() { verb_ = detail::verbose; return *this; }
            Out_log& normal() { verb_ = detail::normal; return *this; }
            Out_log& mute() { verb_ = detail::muted; return *this; }

            // Modification
            void set_verbosity(Verbosity v) { verb_ref_ = v; }

        private:
            Verbosity verb_ref_;
            Verbosity verb_;
            boost::scoped_ptr<File_out> out_;
    };

    // Out_log implementation
    // ========================================================================
    Out_log::Out_log(detail::Parameter* par)
        : verb_ref_(par->verbose_level), verb_(verb_ref_)
    {
        bool useLogfile = not par->log_file.empty();
        if (useLogfile) {
            bool ok = true;
            File_handle file(par->log_file);
            if (bfs::exists(*file) && par->interactive) {
                char c;
                std::cerr << "Logfile: overwrite `" << par->log_file << "'? (y/n)";
                std::cin >> c;
                if (c != 'y') {
                    ok = false;
                }
            }
            if (ok) {
                out_.reset(new File_out(par->log_file));
            }
        }
    }

    template<> inline Out_log& Out_log::operator<<(const Verbosity& v)  
    {
        this->verb_ = v;
        return *this;
    } 

    template<class T> inline Out_log& Out_log::operator<<(const T& x)  
    {
        if (verb_ >= verb_ref_ && verb_ref_ != muted) {
            if (out_.get() != 0) {
                (*out_) << x;
            }
            std::cout << x;
        }
        return *this;
    } 

    template<class T> inline void Out_log::operator()(const T& x) 
    {
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

