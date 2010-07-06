/**
 * @author Roman Pahl
 */

#ifndef permory_helper_functors_hpp
#define permory_helper_functors_hpp

#include "config.hpp"
#include <functional>

namespace Permory 
{
    // Predicates
    // ========================================================================
    // Compares two pairs only based on their second argument
    template<class T> class comp_second : 
        public std::binary_function<T, T, bool>
    {
        public:
            bool operator()(const T& x, const T& y) const {
                return x.second < y.second;
            }
    };

    // True if the second argument of a pair is greater than some value x
    template<class T1, class T2> class greater_than_second : 
        public std::unary_function<T2, bool> 
    {
        public:
            explicit greater_than_second(const T2& x) : arg2(x) {}
            bool operator()(const T1& x) const { return x.second > arg2; }
        private:
            T2 arg2;
    };

    // Operators
    // =======================================================================
} // namespace Permory

#endif // include guard
