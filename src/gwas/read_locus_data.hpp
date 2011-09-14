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

#include "detail/config.hpp"
#include "detail/exception.hpp"
#include "io/format_detect.hpp"
#include "io/line_reader.hpp"
#include "io/input_filters.hpp"

namespace Permory { namespace gwas {
    template<class T> class Locus_data_reader {
        public:
            // Ctor
            Locus_data_reader(
                    const std::string&, //file name
                    char mc='?');       //the character for the missing value

            // Inspection
            bool hasData() const { return !lr_.eof(); }
            detail::Datafile_format get_format() const { return format_; }

            // Conversion
            size_t get_next(std::vector<T>&);

        private:
            detail::Datafile_format format_;
            io::Line_reader<T> lr_;
    };

    // Locus_data_reader<T> implementation
    // ========================================================================
    template<class T> inline Locus_data_reader<T>::Locus_data_reader(
            const std::string& fn, char mc)
        : lr_(fn)
    {
        this->format_ = io::detect_marker_data_format(fn, mc);

        if (format_ == detail::unknown) {
            throw std::runtime_error("Unknown data format.");
        }
    }

    template<> inline size_t Locus_data_reader<char>::get_next(std::vector<char>& v)
    {
        using namespace Permory::detail;
        size_t nskipped = 0;
        v.clear();
        v.reserve(lr_.size());

        while (hasData()) {   
            lr_.next();
            if (*lr_.begin() == '#') {
                nskipped++;
                continue;   //skip comments
            }

            switch(format_) {
                case compact: //straight copy
                    std::copy(lr_.begin(), lr_.end(), std::back_inserter(v)); 
                    break;
                case slide: //white space delimiters
                    std::remove_copy(lr_.begin(), lr_.end(), 
                            std::back_inserter(v), ' ');
                    break; 
                case presto: //skip first two chunks and white space delims thereafter
                    if (*(lr_.begin()) != 'M') {
                        nskipped++;
                        continue;
                    }
                    std::remove_copy_if(lr_.begin(), lr_.end(), 
                            std::back_inserter(v), io::Skip_input_filter<2>());
                    break;
                case plink_tped: //skip first four chunks and white space delims thereafter
                    std::remove_copy_if(lr_.begin(), lr_.end(), 
                            std::back_inserter(v), io::Skip_input_filter<4>());
                    break;
                default:
                    throw std::invalid_argument("Data format not supported.\n");
            }
            if (not v.empty()) {
                break;
            }
        }
        return nskipped;
    }

    // ========================================================================
    //
    // Read loci information from file. Supports PERMORY, PRESTO, PLINK, and 
    // SLIDE where in case of SLIDE, the loci names are simply formed by the 
    // ID. All newly read loci are appended to the loci deque.
    //
    void read_loci(const detail::Datafile_format& format, const std::string& fn, 
            //std::back_insert_iterator<std::deque<Locus> > loci_back)
            std::deque<Locus>* loci)
    {
        using namespace std;
        using namespace Permory::io;
        using namespace Permory::detail;

        io::Line_reader<string> lr(fn);
        string line;
        size_t id = 1;
        if (not loci->empty()) {
            id = loci->back().id() + 1;
        }
        while (not lr.eof()) {
            lr.next_str(line);
            if (line.empty()) {
                continue;
            }
            switch(format) {
                case compact:
                    if (line[0] != '#') {   //skip comments
                        loci->push_back(Locus(id++));
                    }
                    break;
                case presto: //read from *.bgl file as used by presto
                    if (line[0] == 'M') {
                        istringstream iss(line);
                        string rs;
                        iss >> rs; //now rs == 'M'
                        iss >> rs; //now rs contains the markers name
                        loci->push_back(Locus(id++, rs));
                    }
                    break;
                case plink_tped: //read from trans.tped file
                    if (line[0] != '#') {   //skip comments
                        istringstream iss(line);
                        string chr;     //chromosome
                        string rs;      //rs# or snp identifier
                        double cm;      //cM map position
                        size_t bp;      //base pair position in bp units
                        iss >> chr >> rs >> cm >> bp;
                        loci->push_back(Locus(id++, rs, "", string2chr(chr), bp, cm));
                    }
                    break;
                case slide: //create SNP names with serial number
                    if (line[0] != '#') {   //skip comments
                        loci->push_back(Locus(id++));
                    }
                    break;
                default:
                    throw std::invalid_argument("Data format not supported.");
            }
        }
    }
} // namespace gwas
} // namespace Permory

#endif
