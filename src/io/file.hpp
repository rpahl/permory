// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_file_handle_hpp
#define permory_io_file_handle_hpp

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>

#include "detail/config.hpp"
#include "detail/exception.hpp"


namespace Permory { namespace io {
    namespace bfs = boost::filesystem;

    class File_handle {
        public:
            // Ctor
            File_handle(const std::string& fn) 
                : p_(fn) 
            {}

            // Inspection
            const bfs::path& operator*() const { return p_; }
            //bool file_exists() const { return bfs::exists(p_); }
            //bool isRegular() const { return bfs::is_regular(p_); } 
            //size_t file_size() const { return bfs::file_size(p_); }

            // Conversion
            std::string extension(int a=0) const;
        private:
            bfs::path p_;
    };

    // File_handle implementation
    // ========================================================================
    inline std::string File_handle::extension(int a) const
    {
        // example: hello.txt.gz
        // file_extension(0) returns "gz"
        // file_extension(1) returns "txt"
        // file_extension(2 or greater) returns ""
        std::string s = p_.filename().string();
        int x = s.find_last_of(".");
        for (int i=0; i<a; i++) 
        {
            if (x > 0)
                s.erase(x); //remove file extension
            x = s.find_last_of(".");
        }
        if (x > 0)
            return s.substr(x + 1);
        else 
            return "";
    }


    // Deprecated
    class File_info {
        public:
            // Ctor
            File_info(const char* fn) : p_(fn) {}
            bool file_exists() const { return bfs::exists(p_); }
            bool isRegular() const { return bfs::is_regular(p_); } 
            size_t file_size() const { return bfs::file_size(p_); }
            std::string file_extension(int a=0) const;

        private:
            bfs::path p_;
    };

} // namespace io
} // namespace Permory
#endif
