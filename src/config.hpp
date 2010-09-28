// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_config_hpp
#define permory_config_hpp

//#include <cstdlib> 
#include <stdlib.h> //assert
#include <math.h>
#include <iostream>
#include <limits>

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "detail/print.hpp" //mainly debug printing

// Enables large files (> 2GB) at 32-bit architectures
#define _FILE_OFFSET_BITS  64

// A simple print macro 
#define PRINT(X) std::cerr << (#X) << " = "<< (X) << std::endl

namespace Permory 
{
    //typedef std::string string;
    typedef unsigned int uint;
}
#endif

