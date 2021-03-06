// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
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
            if (par.phenotype_domain == Record::dichotomous) {
                i = boost::lexical_cast<int>(v[5]);
                has2 = has2 || (i == 2);
            }
            Record r(boost::lexical_cast<double>(v[5]), par.phenotype_domain);
            if (v[5] == par.undef_phenotype_code) { 
                r.theType = Record::undefined; 
            }
            Individual ind(id++, "", sex);
            ind.add_measurement(r);
            individuals->push_back(ind);
        }
        bool requiresAutocorrection = (par.phenotype_domain == Record::dichotomous && has2);
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

    Record::Value_type determine_phenotype_domain(
            const detail::Parameter& par,
            const std::string& fn)
    {
        using namespace std;
        using namespace Permory::io;
        using namespace Permory::detail;

        Record::Value_type result = Record::undefined;
        Line_reader<string> lr(fn);
        set<double> values; // different values of phenotypes for plink_tfam

        switch (par.phenotype_data_format) {
            case plink_tfam:
                while (not lr.eof()) {
                    lr.next();
                    if (lr.size() < 6) {
                        continue;
                    }
                    vector<string> v(lr.begin(), lr.begin()+6);

                    if (v[5] == par.undef_phenotype_code) {
                        throw std::domain_error("Undefined phenotypes not allowed.");
                    }
                    values.insert(boost::lexical_cast<double>(v[5]));
                    if (values.size() > 2) {
                        result = Record::continuous;
                        break;
                    }
                }
                if (values.size() == 2) {
                    result = Record::dichotomous;
                }
                break;

            case presto:
               while (not lr.eof()) {
                    lr.next();  //read next line
                    if (not lr.empty()) {
                        if (*lr.begin() == "A") {
                            result = Record::dichotomous;
                        }
                        if (*lr.begin() == "T") {
                            result = Record::continuous;
                        }
                    }
                }
                break;

            default:
                throw std::invalid_argument("Phenotype data format not supported.");
        }
        return result;
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
        using namespace Permory::detail;
        Line_reader<string> lr(fn);   
        size_t id = 0;
        if (not individuals->empty()) {
            id = individuals->back().id() + 1;
        }
        switch (par.phenotype_data_format) {
            case plink_tfam:
                read_individuals_from_tfam(par, fn, individuals);
                break;

            case presto:
                while (not lr.eof()) {
                    lr.next();  //read next line
                    if (lr.empty()) {
                        continue;
                    }
                    if (*lr.begin() == "A" || *lr.begin() == "T") {
                        //"A" indicates affection status data
                        //"T" indicates quantitative phenotypes
                        if (*lr.begin() == "A" &&
                            par.phenotype_domain != Record::dichotomous) {
                            throw std::invalid_argument(
                                "Trait indicator needs to be \"A\" for dichotomous phenotypes.");
                        }
                        if (*lr.begin() == "T" &&
                            par.phenotype_domain != Record::continuous) {
                            throw std::invalid_argument(
                                "Trait indicator needs to be \"T\" for quantitative phenotypes.");
                        }
                        // Ignore the "A" or "T" as well as the next entry
                        vector<string> vs(lr.begin()+2, lr.end());  

                        for (size_t i=0; i<vs.size(); i+=2) {
                            Individual ind(id++);
                            double val = boost::lexical_cast<double>(vs[i]);
                            // Presto uses 1 (unaffected) and 2 (affected), but
                            // we use 0 (unaffected) and 1 (affected), thus
                            // subtract 1 if dichotomous phenotypes are used.
                            if (par.phenotype_domain == Record::dichotomous) {
                                val--;
                            }
                            Record r(val, par.phenotype_domain);
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
