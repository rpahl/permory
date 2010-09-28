// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_read_locus_data_hpp
#define permory_io_read_locus_data_hpp

#include <sstream>
#include <string>
#include <vector>
#include <deque>

#include "boost/lexical_cast.hpp"

#include "config.hpp"
#include "detail/exception.hpp"
#include "io/format_detect.hpp"
#include "io/line_reader.hpp"
#include "io/input_filters.hpp"

namespace Permory { namespace io {

    //
    // Read loci informatio from file. Supports PERMORY, PRESTO, PLINK, and 
    // SLIDE where in case of SLIDE, the loci names are simply formed by the 
    // ID. All newly read loci are appended to the loci deque.
    //
    void read_loci(const Datafile_format& format, const std::string& fn, 
            //std::back_insert_iterator<std::deque<Locus> > loci_back)
            std::deque<Locus>* loci)
    {
        Line_reader<std::string> lr(fn);
        std::string line;
        size_t id = 0;
        if (not loci->empty()) {
            id = loci->back().id() + 1;
        }
        while (not lr.eof()) {
            lr.next_str(line);
            if (line.empty()) {
                continue;
            }
            switch(format) {
                //case permory: break; //TODO
                case presto: //read from *.bgl file as used by presto
                    if (line[0] == 'M') {
                        std::istringstream iss(line);
                        std::string rs;
                        iss >> rs; //now rs == 'M'
                        iss >> rs; //now rs contains the markers name
                        Locus loc(id++, rs);
                        loci->push_back(loc);
                    }
                    break;
                case plink: //read from trans.tped file
                    if (line[0] != '#') {   //skip comments
                        std::istringstream iss(line);
                        std::string chr;    //chromosome
                        std::string rs;    //rs# or snp identifier
                        double cm;          //cM map position
                        size_t bp;          //base pair position in bp units
                        iss >> chr >> rs >> cm >> bp;
                        Locus loc(id++, rs, "", string2chr(chr), bp, cm);
                        loci->push_back(loc);
                    }
                    break;
                case slide:
                    if (line[0] != '#') {   //skip comments
                        std::string rs = "SNP";
                        rs.append(boost::lexical_cast<std::string>(id));
                        Locus loc(id++, rs);
                        loci->push_back(loc);
                    }
                    break;
                default:
                    throw std::invalid_argument("Data format not supported.");
            }
        }
    }

    template<class T> class Locus_data_reader {
        public:
            // Ctor
            Locus_data_reader(
                    const std::string&, //file name
                    char mc='?');       //the character for the missing value

            // Inspection
            bool hasData() const { return !lr_.eof(); }
            Datafile_format get_format() const { return format_; }

            // Conversion
            std::vector<T> get_next();

        private:
            Datafile_format format_;
            Line_reader<T> lr_;
    };

    template<class T> inline Locus_data_reader<T>::Locus_data_reader(
            const std::string& fn, char mc)
        : lr_(fn)
    {
        format_ = detect_marker_data_format(fn, mc);

        if (format_ == unknown) {
            throw std::runtime_error("Unknown data format.");
        }
    }

    template<> inline std::vector<char> Locus_data_reader<char>::get_next()
    {
        std::vector<char> v;

        while (hasData()) {   
            lr_.next();
            if (*lr_.begin() == '#')    
                continue;   //skip comments
            v.reserve(lr_.size());

            switch(format_) {
                case permory: //straight copy
                    std::copy(lr_.begin(), lr_.end(), std::back_inserter(v)); 
                    break;
                case slide: //white space delimiters
                    std::remove_copy(lr_.begin(), lr_.end(), 
                            std::back_inserter(v), ' ');
                    break; 
                case presto: //skip first two chunks and then white space delims
                    if (*(lr_.begin()) != 'M')
                        continue;
                    std::remove_copy_if(lr_.begin(), lr_.end(), 
                            std::back_inserter(v), Skip_input_filter<2>());
                    break;
                case plink: //skip first four chunks and then white space delims
                    std::remove_copy_if(lr_.begin(), lr_.end(), 
                            std::back_inserter(v), Skip_input_filter<4>());
                    break;
            }
            if (not v.empty()) {
                break;
            }
        }
        return v;
    }

} // namespace io
} // namespace Permory

#endif
