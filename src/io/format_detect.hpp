// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_format_detect_hpp
#define permory_io_format_detect_hpp

#include "config.hpp"
#include "detail/enums.hpp" //Datafile_format
#include "detail/exception.hpp"
#include "detail/parameter.hpp" 
#include "io/file.hpp" 
#include "io/line_reader.hpp" 

namespace Permory { namespace detail {
    using io::Line_reader;

    bool isInt(char c) {
        int i = int(c);
        return (i >= 48 && i<=57);
    }

    //
    // The format detection works very rudimentary and can be easily fooled to
    // given false positives. On the other hand, it should do the job unless
    // the user attentionally tries to abuse and/or uses non-supported formats.
    //
    Datafile_format detect_marker_data_format(
            const std::string& fn,  //file name
            char mc='?')            //the character for the missing value
    {
        Line_reader<char> lr(fn);
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
                    // either permory or plink
                    if (c3 != ' ') {
                        return permory_data;
                    }
                    else {
                        return plink_tped;
                    }
                }
                if (c2 == ' ') {
                    // either slide or plink
                    if (isInt(c3)) {
                        return slide;
                    }
                    else {
                        return plink_tped; //given that SNP name in plink file starts with non-int
                    }
                }
            }
            else {
                // presto
                if (c1 == 'M' || c1 == 'A' || c1 == 'S' && c2 == ' ') {
                    return presto;
                }
                // check for plink non-autosomal chromosomes
                if ((c1 == 'X' && (c2 == 'Y' || c2 == ' ')) ||  //X or XY
                        (c1 == 'Y' && c2 == ' ') ||                 //Y
                        (c1 == 'M' && c2 == 'T' && c3 == ' ')       //MT
                   ) {
                    return plink_tped;   
                }
            }
        }
        return unknown;
    }

    // 
    // Again, the detection works very rudimentary and can be easily fooled to
    // given false positives but should do the job unless if not abused.
    //
    Datafile_format detect_phenotype_data_format(
            const std::string& fn)  //file name
    {
        try {
            Line_reader<char> lr(fn);
            while (!lr.eof()) {
                lr.next();
                Line_reader<char>::const_iterator il = lr.begin();
                char c1 = *il++;
                char c2 = *il++;
                char c3 = *il++;

                if (c1 == '#') //ignore lines starting with '#' 
                    continue;

                if (isInt(c1)) {
                    if (c2 == ' ') {
                        return plink_tfam;   
                    }
                    else {
                        return permory_data;
                    }
                }
                else {
                    // presto
                    if (c1 == 'M' || c1 == 'A' || c1 == 'S' && c2 == ' ') {
                        return presto;
                    }
                }
            }
            return unknown;
        }
        catch (const std::exception& e) {
            std::cerr << "Error during data format detection: " << e.what() << std::endl;
        }
    }

    Datafile_format detect_meta_data_format(
            const std::string& fn)  //file name
    {
        try {
            Line_reader<char> lr(fn);
            while (!lr.eof()) {
                lr.next();
                Line_reader<char>::const_iterator il = lr.begin();
                char c1 = *il++;
                char c2 = *il++;
                char c3 = *il++;

                if (c1 == '#') {    //ignore lines with '#' 
                    continue;
                }
                else if (isInt(c1)) {
                    if (c2 == ' ') {
                        return plink_tped;   
                    }
                }
                else {
                    if (c1 == 'M' || c1 == 'A' || c1 == 'S' && c2 == ' ') {
                        return presto;
                    }
                    if (c1 == '/' && c2 =='/' && c3 == '!') {
                        return permory_meta;
                    }
                }
            }
            return unknown;
        }
        catch (const std::exception& e) {
            std::cerr << "Error during data format detection: " << e.what() << std::endl;
        }
    }
} //namespace detail
} //namespace Permory

#endif
