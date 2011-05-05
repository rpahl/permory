// Copyright (c) 2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_hooks_hpp
#define permory_detail_hooks_hpp

#include "detail/config.hpp"


//
// Definition of static hooks.
// Hooks are points of entry for custumization of control flow. Static hooks are
// compile time hooks. Only one implementation of a hook is called. To call
// multiple implementations composition of the different implementations is needed.
//
// In this file hooks are defined only. Concrete implementations should be done
// in "detail/config.hpp".
//
// An example for a static hook implementation is support for MPI.
//
namespace Permory {

    namespace hook { // Global definitions for hooks

        // NIL hook.
        class None { };
    }

    namespace hook { // Static hooks for namespace "Permory"

        //
        // Argument hook is called at the start of the application. This hook
        // can be used to alter the arguments passed from commandline.
        //
        template<class T> struct Argument_hook_impl
        {
            public:
                void operator()(int *argc, char ***argv);
        };
        template<> void Argument_hook_impl<None>::operator()
            (int *argc, char ***argv) { }


    } // namespace hook

} // namespace Permory

#endif  // include guard
