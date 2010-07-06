/**
 * @author Roman Pahl (c) 2010
 */

#ifndef permory_config_hpp
#define permory_config_hpp

//#define BOOST_IOSTREAMS_NO_LIB

#define USE_STAT64 0

// Enables large files (> 2GB) at 32-bit architectures
#define _FILE_OFFSET_BITS  64

// A simple print macro 
#define PRINT(X) std::cerr << (#X) << " = "<< (X) << std::endl

//#include <dlib/all/source.cpp>

#include <iostream>

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include "helper/print.hpp" //mainly debug printing
namespace Permory 
{
    typedef std::string string;
}
#endif

