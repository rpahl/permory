/**
 * @author Roman Pahl
 */

#ifndef permory_parameter_hpp
#define permory_parameter_hpp
#include "config.hpp"

namespace Permory 
{
    class Parameter {
        public:
            // statistical tests
            enum Var_type{pooled=0, controls=1, cc=2};
            static Var_type var;                //variance estimator
            static double w[3];                 //weights of CA trend test

            // permutation
            static size_t nperm;                //number of permutations
            
            // permutation optimization
            static size_t tail_sz;            //size of tail (PAM method)
            static bool useBar;                      //bit arithmetics flag
    };

    // Declare static variables
    // ========================
    
    // statistical tests
    Parameter::Var_type Parameter::var = Parameter::pooled;
    double Parameter::w[3] = {0, 1, 2};

    // permutation
    size_t Parameter::nperm = 0;
    size_t Parameter::tail_sz = 100;

    // permutation optimization
    bool Parameter::useBar = true;
}

#endif

