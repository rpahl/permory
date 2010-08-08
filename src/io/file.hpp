// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_file_handle_hpp
#define permory_io_file_handle_hpp

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>

#include "config.hpp"
#include "detail/exception.hpp"


namespace Permory { namespace io {
    namespace bfs = boost::filesystem;

    class File_handle {
        public:
            // Ctor
            File_handle(const char* fn) 
                : p_(fn, bfs::native) 
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


    // Deprecated
    class File_info {
        public:
            // Ctor
            File_info(const char* fn) : p_(fn, bfs::native) {}
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
