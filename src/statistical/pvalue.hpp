/**
 * @author Roman Pahl
 */

#ifndef permory_pvalue_hpp
#define permory_pvalue_hpp

#include <iostream>
#include <valarray>
#include <vector>
#include <gsl/gsl_cdf.h>
#include "config.hpp"
#include "helper/vector.hpp"

namespace Permory 
{
    typedef std::vector<double> vector_double;
    typedef std::vector<std::vector<double> > matrix_type;

    // asymptotic
    double quisq(double, double); // FIXME just use the gsl function
    
    // Computes single step p-values (or counts)
    std::vector<double> single_step_pvalues(
            const std::vector<double>&, //test statistics sorted _decreasingly_ 
            const std::vector<double>&, //max test statistic per permutation
            bool cc=true);              //if true, raw counts are returned

    // Computes step down p-values (or counts)
    std::vector<double> step_down_pvalues(
            const std::vector<double>&,//top X test stats sorted _decreasingly_
            const matrix_type&,//max teststats per permutation for top X markers
            const std::vector<double>&,//condensed max test stats for the 
            //non-top markers
            bool cc=true);  //if true, raw counts are returned
} // namespace Permory

#endif // include guard
