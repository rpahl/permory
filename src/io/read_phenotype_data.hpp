// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_io_read_phenotype_data_hpp
#define permory_io_read_phenotype_data_hpp

#include <sstream>
#include <string>
#include <vector>
#include <deque>

#include "boost/lexical_cast.hpp"

#include "config.hpp"
#include "detail/exception.hpp"
#include "detail/parameter.hpp"
#include "individual.hpp"
#include "io/format_detect.hpp"
#include "io/line_reader.hpp"
#include "io/input_filters.hpp"

namespace Permory { namespace io {

    //
    // Read individuals from PLINK's *.tfam file. This function is called by
    // the function 'read_individuals' (see below).
    //
    void read_individuals_from_tfam(
            const Parameter& par,   //data format, missing code
            const std::string& fn,  //file name 
            std::vector<Individual>* individuals)
    {
        size_t id = 0;
        if (not individuals->empty()) 
            id = individuals->back().id() + 1;

        Line_reader<std::string> lr(fn);   
        bool has2 = false;//for auto correction of affection status coding

        while (not lr.eof()) {
            lr.next();
            std::vector<std::string> v(lr.begin(), lr.end());
            if (v.size() < 6) {
                throw std::runtime_error("Less than six entries in PLINK *.tfam file.");
            }
            // Determine sex
            Individual::Sex sex = Individual::nosex;
            int i = boost::lexical_cast<int>(v[4]);
            if (i == 1) {
                sex = Individual::male;
            }
            else if (i == 2) {
                sex = Individual::female;
            }
            // Get the phenotype
            i = boost::lexical_cast<int>(v[5]);
            has2 = (i == 2);
            Record r(boost::lexical_cast<double>(v[5]), par.val_type);
            if (v[5] == par.undef_phenotype_code) { 
                r.theType = Record::undefined; 
            }
            Individual ind(id++, "", sex);
            ind.add_measurement(r);
            individuals->push_back(ind);
        }
        bool requiresAutocorrection = (par.val_type == Record::dichotomous && has2);
        if (requiresAutocorrection) {
            // Auto correct affection status as follows:
            //  unaffected: 1 -> 0
            //  affected:   2 -> 1
            BOOST_FOREACH(Individual i, *individuals) {
                Individual::iterator it = i.end()-1; //iterator to last record
                if (it->val == 2) {
                    it->val = 1;
                }
                else if (it->val == 1) {
                    it->val = 0;
                }
            }
        }
    }

    //
    // Read individuals from file. Supports PERMORY, PRESTO, and PLINK. If 
    // vector contains individuals, the newly read individuals are appended.
    //
    void read_individuals(const Parameter& par, const std::string& fn,      
            std::vector<Individual>* individuals)
    {
        Line_reader<std::string> lr(fn);   
        size_t id = 0;
        if (not individuals->empty()) {
            id = individuals->back().id() + 1;
        }
        switch (par.phenotype_data_format) {
            //case permory: break; //TODO
            case plink:
                read_individuals_from_tfam(par, fn, individuals);
                break;
            case presto:
                while (not lr.eof()) {
                    lr.next();  //read next line
                    if (*lr.begin() == "A") { //"A" indicates affection status data
                        // Ignore the "A" as well as the next entry 
                        std::vector<std::string> vs(lr.begin()+2, lr.end());  

                        for (int i=0; i<vs.size(); i++) {
                            Individual ind(id++);
                            double val = boost::lexical_cast<double>(vs[i]);
                            // Presto uses 1 (unaffected) and 2 (affected), but
                            // we use 0 (unaffected) and 1 (affected), thus -1
                            val--;  
                            Record r(val, par.val_type);
                            ind.add_measurement(r);
                            individuals->push_back(ind);
                            if (par.gen_type == haplotype) {
                                i++;
                            }
                        }
                        break;
                    }
                }
                break;
            default:
                throw std::invalid_argument("Data format not supported.");
        }
    }

} // namespace io
} // namespace Permory

#endif
