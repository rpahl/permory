/**
 * @author Roman Pahl
 */

#ifndef permory_io_output_filters_hpp
#define permory_io_outut_filters_hpp

#include <string.h>

#include <boost/iostreams/char_traits.hpp> //EOF, WOULD_BLOCK
#include <boost/iostreams/concepts.hpp>    //multichar_input_filter
#include <boost/iostreams/operations.hpp>  //get

#include "config.hpp"

namespace Permory 
{
    namespace bio = boost::iostreams;

    class Delim_output_filter : public bio::multichar_output_filter {
        public:
            explicit Delim_output_filter(char delim = ' ')
                : delim_(delim) {}

            template<typename Sink> std::streamsize write(
                    Sink& dest, const char* s, std::streamsize n)
            {
                std::streamsize i;
                for (i = 0; i < n; ++i) {
                    int c = s[i];
                    if (bio::put(dest, c)) {
                        if (c != '\n') { //no delimiter after last char
                            if (!bio::put(dest, delim_))
                                break;
                        }
                    }
                    else
                        break;
                }
                return i; 
            }

            template<typename Source> void close(Source&) { 
            } 
        private:
            char delim_;
    };
} // namespace Permory

#endif

