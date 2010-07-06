// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#include "file_info.hpp"
using namespace std;
using namespace Permory;

std::string File_info::file_extension(int a) const
{
    // example: hello.txt.gz
    // file_extension(0) returns "gz"
    // file_extension(1) returns "txt"
    // file_extension(2 or greater) returns ""
    std::string s = (std::string) p_.filename(); 
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

File_handler::File_handler(const char* fn) : p_(fn, bfs::native) 
{
    if (!bfs::exists(p_)) { 
        throw File_exception("not found.");
    }
    if ( !bfs::is_regular(p_)) {
        throw File_exception("not a regular file.");
    }
}

std::string File_handler::extension(int a) const
{
    // example: hello.txt.gz
    // file_extension(0) returns "gz"
    // file_extension(1) returns "txt"
    // file_extension(2 or greater) returns ""
    std::string s = (std::string) p_.filename(); 
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

