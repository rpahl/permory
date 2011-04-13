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

#include "detail/config.hpp"
#include "detail/exception.hpp"
#include "detail/parameter.hpp"
#include "individual.hpp"
#include "io/format_detect.hpp"
#include "io/line_reader.hpp"
#include "io/input_filters.hpp"

namespace Permory { namespace gwas {
    //
    // Read individuals from PLINK's *.tfam file. This function is called by
    // the function 'read_individuals' (see below).
    //
    void read_individuals_from_tfam(
            const detail::Parameter& par,   //data format, missing code
            const std::string& fn,  //file name 
            std::vector<Individual>* individuals)
    {
        using namespace std;
        using namespace Permory::io;
        using namespace Permory::detail;
        size_t id = 0;
        if (not individuals->empty()) 
            id = individuals->back().id() + 1;

        Line_reader<string> lr(fn);   
        bool has2 = false;//for auto correction of affection status coding

        while (not lr.eof()) {
            lr.next();
            if (lr.size() < 6) {
                continue;
            }
            vector<string> v(lr.begin(), lr.begin()+6);

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
            if (par.val_type == Record::dichotomous) {
                i = boost::lexical_cast<int>(v[5]);
                has2 = has2 || (i == 2);
            }
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
            for (vector<Individual>::iterator itInd = individuals->begin(); 
                    itInd != individuals->end(); ++itInd) {
                // get iterator to last added record
                Individual::iterator itRecord = itInd->end()-1; 
                if (itRecord->val == 2.0) {
                    itRecord->val = 1.0;
                }
                else if (itRecord->val == 1.0) {
                    itRecord->val = 0.0;
                }
            }
        }
    }

    //
    // Read individuals from file. Supports PERMORY, PRESTO, and PLINK. If 
    // vector contains individuals, the newly read individuals are appended.
    //
    void read_individuals(const detail::Parameter& par, const std::string& fn,      
            std::vector<Individual>* individuals)
    {
        using namespace std;
        using namespace Permory::io;
        Line_reader<string> lr(fn);   
        size_t id = 0;
        if (not individuals->empty()) {
            id = individuals->back().id() + 1;
        }
        switch (par.phenotype_data_format) {
            case compact:
                while (not lr.eof()) {
                    lr.next();  //read next line
                    if (lr.empty()) {
                        continue;
                    }
                    if (*lr.begin() != "#") { //skip comments
                        vector<string> vs(lr.begin(), lr.end());  

                        for (size_t i=0; i<vs.size(); i++) {
                            Individual ind(id++);
                            double val = //0 means unaffected and 1 affected
                                boost::lexical_cast<double>(vs[i]);
                            Record r(val, par.val_type);
                            ind.add_measurement(r);
                            individuals->push_back(ind);
                        }
                        break;
                    }
                }
                break;
            case plink_tfam:
                read_individuals_from_tfam(par, fn, individuals);
                break;

            case presto:
                while (not lr.eof()) {
                    lr.next();  //read next line
                    if (lr.empty()) {
                        continue;
                    }
                    if (*lr.begin() == "A") { //"A" indicates affection status data
                        // Ignore the "A" as well as the next entry 
                        vector<string> vs(lr.begin()+2, lr.end());  

                        for (size_t i=0; i<vs.size(); i+=2) {
                            Individual ind(id++);
                            double val = boost::lexical_cast<double>(vs[i]);
                            // Presto uses 1 (unaffected) and 2 (affected), but
                            // we use 0 (unaffected) and 1 (affected), thus
                            // subract 1
                            val--;  
                            Record r(val, par.val_type);
                            ind.add_measurement(r);
                            individuals->push_back(ind);
                        }
                        break;
                    }
                }
                break;
            default:
                throw std::invalid_argument("Phenotype data format not supported.");
        }
    }
} // namespace gwas
} // namespace Permory

#endif
