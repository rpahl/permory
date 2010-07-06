/**
 * @author Roman Pahl
 */

#ifndef permory_permutation_perm_hpp
#define permory_permutation_perm_hpp

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "exception.hpp"
#include "helper/matrix.hpp"
#include "helper/vector.hpp"

namespace Permory 
{
    class Permutation
    {
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
} // namespace Permory

#endif // include guard


