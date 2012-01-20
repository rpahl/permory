// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_permutation_permutation_hpp
#define permory_permutation_permutation_hpp

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "detail/exception.hpp"
#include "detail/matrix.hpp"
#include "detail/vector.hpp"

namespace Permory { namespace permutation {

    //
    //  Providing permutating/shuffling function. Here we prefer the GSL over 
    //  the boost random generator interface as (at the time of implementation) 
    //  it was considerably more efficient
    //
    class Permutation {
        public:
            // Ctor
            explicit Permutation(size_t seed=17061979) : seed_(seed) {
                rg = gsl_rng_alloc(gsl_rng_mt19937); //Mersenne Twister
                gsl_rng_set(rg, seed); //init random generator with seed
            }
            // Dtor
            ~Permutation() { gsl_rng_free (rg); }

            // Inspector
            size_t get_seed() const { return seed_; }
            template<typename T> void shuffle(T* pv, size_t sz) const {
                gsl_ran_shuffle(rg, pv, sz, sizeof(T));
            }

            // Modifier
            void set_seed(size_t x) { seed_ = x; gsl_rng_set(rg, seed_); }
            void reset_seed() { gsl_rng_set(rg, seed_); }

        private:
            size_t seed_;   //random seed
            gsl_rng* rg;    //random generator
    };
} // namespace permutation
} // namespace Permory

#endif // include guard


