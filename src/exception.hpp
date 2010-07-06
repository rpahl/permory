

/**
 * @author Roman Pahl
 */

#ifndef permory_exception_file_hpp
#define permory_exception_file_hpp

#include <exception>
#include <stdexcept>
#include "config.hpp"

namespace Permory {
    class File_exception : public std::exception {
        public:
            explicit File_exception(const char* s) : s_(s) {}
            virtual const char* what() const throw() { return s_; }
        private:
            const char* s_;
    };

    class Dimension_error : public std::exception { 
        public:
            explicit Dimension_error(const char* s) : s_(s) {}
            virtual const char* what() const throw() { return s_; }
        private: 
            const char* s_;
    };

    class Math_error {
    };

} //namespace Permory

#endif
