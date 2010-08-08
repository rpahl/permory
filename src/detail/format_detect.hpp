// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_format_detect_hpp
#define permory_detail_format_detect_hpp

#include "config.hpp"
#include "detail/exception.hpp"
#include "detail/parameter.hpp" //Data_format
#include "io/file.hpp" 
#include "io/line_reader.hpp" 

namespace Permory { namespace detail {
    using io::Line_reader;

    bool isInt(char c) {
        int i = int(c);
        return (i >= 48 && i<=57);
    }

    // The format detection works very rudimentary and can be easily fooled to
    // given false positives. On the other hand, it should do the job unless
    // the user attentionally tries to abuse and/or uses non-supported formats.
    enum Data_format{unknown=0, simple=1, slide=2, beagle=3, plink=4};
    Data_format detect_data_format(
            const io::File_handle& f, 
            char mc='?') //the character for the missing value
    {
        Line_reader<char> lr(f);
        while (!lr.eof()) {
            lr.next();
            Line_reader<char>::const_iterator il = lr.begin();
            char c1 = *il++;
            char c2 = *il++;
            char c3 = *il++;

            if (c1 == '#') //ignore lines starting with '#' 
                continue;

            if (c1 == mc || isInt(c1)) {
                if (c2 == mc || isInt(c2)) {
                    // either simple or plink
                    if (c3 != ' ') 
                        return simple;
                    else 
                        return plink;
                }
                if (c2 == ' ') {
                    // either slide or plink
                    if (isInt(c3)) 
                        // assume SNP name in plink file starts with non-int
                        return slide;
                    else
                        return plink;
                }
            }
            else {
                // beagle
                if (c1 == 'M' || c1 == 'A' || c1 == 'S' && c2 == ' ') {
                    return beagle;
                }
                // check for plink non-autosomal chromosomes
                if ((c1 == 'X' && (c2 == 'Y' || c2 == ' ')) ||  //X or XY
                    (c1 == 'Y' && c2 == ' ') ||                 //Y
                    (c1 == 'M' && c2 == 'T' && c3 == ' ')       //MT
                   )
                    return plink;
            }

        }
        return unknown;
    }

} //namespace detail
} //namespace Permory

#endif
