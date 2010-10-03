// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_parameter_hpp
#define permory_detail_parameter_hpp

#include <set>
#include <string>

#include "detail/config.hpp"
#include "detail/enums.hpp" 
#include "individual.hpp" //Record::Value_type

namespace Permory { namespace detail {

    class Parameter {
        public:
            //
            // General
            //
            static double version;                  //program version

            //
            // Analysis
            //
            static double min_maf;                  //minor allele freq threshold
            static double max_maf;                  //minor allele freq threshold
            static Genetic_type gen_type;           //haplotype or genotype analysis
            static size_t m;                        //number of markers
            static size_t ncase, ncontrol;          //number of cases/controls

            //
            // Data
            //
            // Supported data file formats: PERMORY, PRESTO, PLINK, and SLIDE
            static Datafile_format phenotype_data_format;     
            static Record::Value_type val_type;     //dichotomous or continous 
            static char undef_allele_code;          //code of missing allele
            static std::string undef_phenotype_code;//code of missing phenotype
            
            //
            // Input
            //
            static std::set<std::string> fn_marker_data;//data file names
            static std::set<std::string> fn_bad_files;  //files of unknown format
            static std::string fn_trait;        //trait/phenotype file name
            //static std::string fn_meta;         //meta information file name

            //
            // Output
            //
            static bool interactive;        //ask before overwriting files 
            static bool pval_counts;        //output raw "p-value counts" yes/no
            static bool quiet;              //no output to console
            static bool verbose;            //verbose output 
            static Verbosity verbose_level; //verbose, normal, or muted
            static size_t ntop;             //show the top n results
            static std::string out_prefix;  //to derive output file names 
            static std::string log_file;    //log console output to file

            //
            // Statistical testing
            //
            static Var_estimate ve;             //variance estimator
            static double w[3];                 //weights of CA trend test
            static std::set<Test_type> tests;   //statistical tests

            //
            // Permutation
            //
            int seed;                   //random seed;
            static size_t nperm_total;  //total number of permutations
            static size_t nperm_block;  //block-wise number of permutations

            // speed optimization
            static size_t tail_size;    //size of tail (REM method)
            static bool useBar;         //use bit arithmetics yes/no

    };

    // Declare static variables
    // ========================
    
    //
    // General
    //
    double Parameter::version = 1.0;

    //
    // Analysis
    //
    double Parameter::min_maf = 0.0;
    double Parameter::max_maf = 1.0;
    Genetic_type Parameter::gen_type = genotype;
    size_t Parameter::m = 0;
    size_t Parameter::ncase = 0;
    size_t Parameter::ncontrol = 0;

    //
    // Data
    //
    Datafile_format Parameter::phenotype_data_format = unknown;
    char Parameter::undef_allele_code = '?';
    std::string Parameter::undef_phenotype_code = "?"; 
    Record::Value_type Parameter::val_type = Record::dichotomous;     

    //
    // Input
    //
    std::set<std::string> Parameter::fn_marker_data; 
    std::set<std::string> Parameter::fn_bad_files; 
    std::string Parameter::fn_trait = "";                
    //std::string Parameter::fn_meta = "";                

    //
    // Output
    //
    bool Parameter::interactive = true;
    bool Parameter::pval_counts = false;    
    bool Parameter::quiet = false;
    bool Parameter::verbose = false;
    Verbosity Parameter::verbose_level = detail::normal;
    size_t Parameter::ntop = 100;
    std::string Parameter::out_prefix = "out";  
    std::string Parameter::log_file = "";  

    //
    // Statistical testing
    //
    Var_estimate Parameter::ve = pooled;
    double Parameter::w[3] = {0, 1, 2};
    std::set<Test_type> Parameter::tests;    

    //
    // Permutation
    //
    size_t Parameter::nperm_total = 10000;
    size_t Parameter::nperm_block = 10000;
    size_t Parameter::tail_size = 100;
    bool Parameter::useBar = true;

} //namespace detail
} //namespace Permory

#endif

