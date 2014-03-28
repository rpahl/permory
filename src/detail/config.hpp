// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_config_hpp
#define permory_detail_config_hpp

#include <stdlib.h> //assert
#include <math.h>
#include <iostream>
#include <limits>

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/utility.hpp>

#include "detail/print.hpp" //mainly debug printing
#include "detail/hooks.hpp"

#ifdef USE_MPI
#include "mpi/mpi.hpp"
#endif  // USE_MPI

// A simple print macro 
#define PRINT(X) std::cerr << (#X) << " = "<< (X) << std::endl

// Enables large files (> 2GB) at 32-bit architectures
#define _FILE_OFFSET_BITS  64

// Ensure locale independence using boost::lexical_cast
#define BOOST_LEXICAL_CAST_ASSUME_C_LOCALE 1

namespace Permory 
{
    typedef unsigned int uint;

    namespace hook {
        //
        // Activation of global static hooks for namespace "Permory".
        //

        //
        // Argument_hook
        //   Default: None
        //
        // Activation is done in the USE_MPI cases below.
        //struct Argument_hook : public Argument_hook_impl<None> { };
    } // namespace hook


    #ifdef USE_MPI
        namespace hook {
            struct Argument_hook : public Argument_hook_impl<MPI> { };
        }

        namespace gwas { class Mpi_analyzer_factory; }
        typedef gwas::Mpi_analyzer_factory analyzer_factory_t;
    #else
        namespace hook {
            struct Argument_hook : public Argument_hook_impl<None> { };
        }

        namespace gwas { class Default_analyzer_factory; }
        typedef gwas::Default_analyzer_factory analyzer_factory_t;
    #endif  // USE_MPI
}
#endif

