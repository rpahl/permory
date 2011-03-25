// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#include <iostream>
#include <iomanip>
#include <time.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/progress.hpp>   //timer

#include "detail/config.hpp"
#include "detail/parameter.hpp"
#include "gwas/analysis.hpp"
#include "io/file.hpp"
#include "io/line_reader.hpp"
#include "io/output.hpp"

int main(int ac, char* av[])
{
    namespace cls = boost::program_options::command_line_style;
    using namespace std;
    using namespace boost;
    using namespace boost::program_options;
    using namespace Permory;
    using namespace Permory::detail;
    using namespace Permory::gwas;
    using namespace Permory::io;

    timer t;    //t.restart(); //start clock
    time_t rawtime; 
    struct tm * timeinfo;   //time and date
    Parameter par;
    Myout myout(&par); 
    string config_file;
    std::vector<string> marker_data_files;

    analyzer_factory_t factory(ac, av);

    myout << endl;
    myout << "+-----------------+-----------------+------------------+" << endl;
    myout << "|    PERMORY      |      v"<< fixed << setprecision(2) << 
        par.version <<                     "      |    6/Oct/2010    |" << endl;
    myout << "+-----------------+-----------------+------------------+" << endl;
    myout << "|           Copyright (c) 2010 Roman Pahl              |" << endl;
    myout << "|  Distributed under the Boost Software License, v1.0  |" << endl;
    myout << "+------------------------------------------------------+" << endl;
    myout << "|                   www.permory.org                    |" << endl;
    myout << "+------------------------------------------------------+" << endl;

    try {
        //
        // General
        //
        options_description general("General");
        general.add_options()
            ("help,h", "show important options")
            ("help-adv", "show advanced options")
            ("help-all,H", "show all options")
            ("interactive,i", "prompt before overwriting files")
            ("quiet", "suppress console output")
            ("verbose,v", "detailed output")
            ("debug,d", "more detailed output")
            ;
        //
        // Analysis
        //
        options_description analysis("Analysis");
        analysis.add_options()
            ("nco", value<size_t>(&par.ncontrol), 
             "number of controls - if specified, shadows option "
             "'--trait-file' and requires option '--nca'\n" 
             "Assumes: data format [...controls...|...cases...]")
            ("nca", value<size_t>(&par.ncase), 
             "number of cases (requires option '--nco')")
            ("min-maf", 
             value<double>(&par.min_maf)->default_value(0.0),
             "lower minor allele frequency threshold") 
            ("max-maf", 
             value<double>(&par.max_maf)->default_value(0.5),
             "upper minor allele frequency threshold")
            ;
        //
        // Data
        //
        options_description data("Data");
        data.add_options()
            ("allelic", "allelic-based permutation (default: genotype)") 
            ("missing",  
             value<char>(&par.undef_allele_code)->default_value('?'),
             "code of missing marker data (single character)")
            ;
        //
        // I/O
        //
        options_description io("Input/Output");
        io.add_options()
            ("trait-file,f", value<string>(&par.fn_trait), 
             "file with binary trait data")
            ("out-prefix,o", 
             value<string>(&par.out_prefix)->default_value("out"), 
             "prefix for all output files") 
            ("config,c", value<string>(&config_file), "configuration file")
            ("log", value<string>(&par.log_file), "log file")
            //("stat,s", value<string>(&tests), "statistical tests to be used")
            ;
        //
        // Permutation
        //
        options_description perm("Permutation");
        perm.add_options()
            ("nperm,n",
             value<size_t>(&par.nperm_total)->default_value(10000),
             "Number of permutations")
            ("seed", 
             value<int>(&par.seed)->default_value(12345678), 
             "random seed")
            ;
        //
        // Advanced
        //
        // Statistical testing
        // speed optimization
        options_description advanced("Advanced options");
        advanced.add_options()
            ("counts", "output counts in addition to p-values")
            ("block",value<size_t>(&par.nperm_block)->default_value(10000),  
             "permutation block size")
            ("ntop", value<size_t>(&par.ntop)->default_value(100), 
             "number of top markers shown in the *.top output file")
            ("tail", value<size_t>(&par.tail_size)->default_value(100), 
             "size of sliding tail (REM method)")

            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        options_description hidden("Hidden");
        hidden.add_options()
            ("data-file", value<vector<string> >(&marker_data_files)->composing(), 
             "input data file")
            ;

        // Create option "profiles"
        options_description cmd_line_options;
        cmd_line_options.add(general).add(analysis).add(data).add(io).add(perm).
            add(advanced).add(hidden);

        options_description config_file_options;
        config_file_options.add(general).add(analysis).add(data).add(io).add(perm).
            add(advanced).add(hidden);

        options_description visible("Important options");
        visible.add(general).add(analysis).add(data).add(io).add(perm);

        // data-file is specified without option flag
        positional_options_description pos;
        pos.add("data-file", -1);

        // parse command line
        variables_map vm;
        store(command_line_parser(ac, av).
                options(cmd_line_options).
                positional(pos).
                // turn off guessing, otherwise --help-all shadows --help
                style(command_line_style::default_style ^ command_line_style::allow_guessing).
                run(), vm);


        if (vm.count("help")) {
            cout << "Usage:\n\tpermory [options] <data_file1> [data_file2 ...]\n\n"; 
            cout << visible << endl;
            return 0;
        }
        if (vm.count("help-all")) {
            cout << "Usage:\n\tpermory [options] <data_file1> [data_file2 ...]\n\n"; 
            cout << visible << endl;
            cout << advanced << endl;
            return 0;
        }
        if (vm.count("help-adv")) {
            cout << advanced << endl;
            return 0;
        }
        notify(vm);

        bool hasConfigFile = vm.count("config") > 0;
        if (hasConfigFile) {
            ifstream ifs(config_file.c_str());
            if (!ifs) {
                cerr << errpre << "Unable to open config file: " << config_file << endl;
                return 1;
            }
            else {
                // in any case, config file is shadowed by cmd_line inputs
                store(parse_config_file(ifs, config_file_options), vm);
                notify(vm);
            }
        }
        par.quiet = vm.count("quiet") > 0;
        par.interactive = vm.count("interactive") > 0;
        par.verbose = vm.count("verbose") > 0;
        par.debug = vm.count("debug") > 0;
        if (par.quiet) {
            myout.set_verbosity(muted);
        }
        else if (par.debug) {
            myout.set_verbosity(all);
        }
        else if (par.verbose) {
            myout.set_verbosity(verbose);
        }
        if (not par.log_file.empty()) {
            myout.set_logfile(par.log_file, par.interactive);
            myout << normal << stdpre << "Writing this to log file `" << 
                par.log_file<<"'." << endl;
        }

        par.fn_marker_data.insert(marker_data_files.begin(), marker_data_files.end());
        bool hasData = vm.count("data-file") > 0;
        if (!hasData) {
            cerr << errpre << "Please specify a data-file." << endl;
            cerr << errpre << "For more information try ./permory --help" << endl;
            return 1;
        }
        else {
            set<string> failed;
            BOOST_FOREACH(string fn, par.fn_marker_data) {
                try {
                    Line_reader<char> lr(fn);
                }
                catch (const std::exception& e) {
                    cerr << errpre << fn << ": Could not open file - will be ignored" << endl;
                    failed.insert(fn);
                }
            }
            // now discard the bad files (could not be done during loop)
            BOOST_FOREACH(string s, failed) {
                par.fn_marker_data.erase(s);
            }
        }
        if (par.fn_marker_data.size() < 1) {
            cerr << errpre << "No valid data file." << endl;
            return 1;
        }

        bool hasNco = vm.count("nco") > 0;
        bool hasNca = vm.count("nca") > 0;
        bool hasTraitFile = vm.count("trait-file") > 0;
        bool createTrait = (hasNco || hasNca);
        if (createTrait) {
            bool ok = true;
            if (not hasNco) {
                cerr << errpre << "Option '--nco' is missing." << endl;
                ok = false;
            }
            if (not hasNca) {
                cerr << errpre << "Option '--nca' is missing." << endl;
                ok = false;
            }
            if (par.ncontrol == 0 || par.ncase == 0) {
                cerr << errpre << "Number of controls or cases cannot be 0." << endl;
                ok = false;
            }
            if (ok) {
                par.fn_trait = ""; //do not read trait file anymore
            }
            else if (hasTraitFile) {    //last chance
                cerr << errpre << "Will try to read trait from " << par.fn_trait << "." << endl;
            }
            else {
                return 1;   //failure
            }
        }
        else {
            if (hasTraitFile) {
                // check file 
                try {
                    Line_reader<char> lr(par.fn_trait);
                }
                catch (const std::exception& e) {
                    cerr << e.what() << endl;
                    return 1;
                }
            }
            else {
                cerr << errpre << "Please specify either trait file (option " <<
                    "'--trait-file') or both number of controls ('--nco') " <<
                    "AND cases ('--nca')." << endl; 
                return 1;
            }
        }

        if (vm.count("counts")) {
            par.pval_counts = true;
        }
        if (vm.count("allelic")) {
            par.marker_type = allelic;
            par.tests.insert(chisq);
        }
        else {
            par.tests.insert(trend);
        }

        if (par.min_maf >= par.max_maf) {
            cerr << errpre << "Bad maf filter specification will left no " <<
                "markers to analyze." << endl;
            return 1;
        }

        time(&rawtime);
        timeinfo = localtime(&rawtime);
        myout << normal << stdpre << "Started at " << asctime(timeinfo) << endl;

        myout << normal << stdpre << "Data file(s):" << endl;
        BOOST_FOREACH(string fn, par.fn_marker_data) {
            myout << indent(4) << fn << endl;
        }
        myout << stdpre<< "Output to: " << par.out_prefix << ".*" << endl;
        if (createTrait) {
            myout.verbose();
            myout << stdpre << "User specified case/control numbers:" << endl;
            myout << stdpre << "Number of controls: " << par.ncontrol << endl;
            myout << stdpre << "Number of cases: " << par.ncase << endl;
        }
        myout << normal << stdpre << "Number of permutations: " << 
            par.nperm_total << endl << endl;
    }
    catch(std::exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }    

    try {
        gwas_analysis(&par, myout, factory);
    }
    catch(std::exception& e)
    {
        cerr << errpre << "Error: " << e.what() << endl;
        return 1;
    }
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    myout << stdpre << "Finished at " << asctime(timeinfo);
    myout << stdpre << "Runtime: " << t.elapsed() << " s" << endl;
    return 0;
}

