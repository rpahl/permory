// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#include <iostream>
#include <iomanip>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/progress.hpp>

#include "detail/config.hpp"
#include "detail/parameter.hpp"
#include "gwas.hpp"
#include "io/file.hpp"
#include "io/line_reader.hpp"
#include "io/output.hpp"

using namespace std;

// A helper function to simplify the main part.
template<class T> ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(cout, " ")); 
    return os;
}

namespace cls = boost::program_options::command_line_style;
int main(int ac, char* av[])
{
    using namespace std;
    using namespace boost;
    using namespace boost::program_options;
    using namespace Permory;
    using namespace Permory::detail;
    using namespace Permory::io;

    timer t; //t.restart();
    Parameter par;
    Out_log myout(&par); 
    string config_file;
    std::vector<string> marker_data_files;

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
            ("verbose,v", "detailed status output")
            ("version", "output the version number")
            ;
        //
        // Analysis
        //
        options_description analysis("Analysis");
        analysis.add_options()
            ("nco", value<size_t>(&par.ncontrol), 
             "number of controls - if specified, shadows option "
             "'--trait-file' and requires option '--nca'\n" 
             "    Assumes: data format [...controls...|...cases...]")
            ("nca", value<size_t>(&par.ncase), 
             "number of cases (requires option '--nco')")
            ("min-maf", 
             value<double>(&par.min_maf)->default_value(0.0),
             "only consider markers with minor allele frequency greater "
             "than min-maf")
            ;
        //
        // Data
        //
        options_description data("Data");
        data.add_options()
            ("haplo", "assume haplotype data format (default: genotype)") 
            ("missing,m",  
             value<char>(&par.undef_allele_code)->default_value('?'),
             "code for missing marker data")
            ;
        //
        // I/O
        //
        options_description io("Input/Output");
        io.add_options()
            ("trait-file,f", value<string>(&par.fn_trait), 
             "file containing binary trait data")
            ("out-prefix,o", 
             value<string>(&par.out_prefix)->default_value("out"), 
             "prefix used for all output files") 
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
            ("tail", 
             value<size_t>(&par.tail_size)->default_value(100), 
             "size of sliding tail (REM method)")
            //("ndisp", value<int>()->default_value(100), "number of top markers to be displayed in the *.top output file")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        options_description hidden("Hidden");
        hidden.add_options()
            //("data-file", value< set<string> >(&par.fn_marker_data)->composing(), 
            ("data-file", value< vector<string> >(&marker_data_files)->composing(), 
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
            cout << "Usage:\n\tpermory [options] data-file1 [data-file2 ...]\n\n"; 
            cout << visible;
            return 0;
        }
        if (vm.count("help-all")) {
            cout << "Usage:\n\tpermory [options] data-file1 [data-file2 ...]\n\n"; 
            cout << visible << endl;
            cout << advanced << endl;
            return 0;
        }
        if (vm.count("help-adv")) {
            cout << advanced << endl;
            return 0;
        }
        if (vm.count("version")) {
            cout.precision(2);
            cout << "PERMORY version: " << showpoint << par.version << endl;
            return 0;
        }
        notify(vm);

        bool hasConfigFile = vm.count("config") > 0;
        if (hasConfigFile) {
            ifstream ifs(config_file.c_str());
            if (!ifs) {
                cerr << "!!! Unable to open config file: " << config_file << endl;
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
        if (par.quiet) {
            myout.set_verbosity(muted);
        }
        else if (par.verbose) {
            myout.set_verbosity(verbose);
        }
        if (not par.log_file.empty()) {
            myout.set_logfile(par.log_file, par.interactive);
        }

        par.fn_marker_data.insert(marker_data_files.begin(), marker_data_files.end());
        bool hasData = vm.count("data-file") > 0;
        if (!hasData) {
            cerr << "!!! Please specify a data-file." << endl;
            cerr << "!!! For more information try ./permory --help" << endl;
            return 1;
        }
        else {
            set<string> failed;
            BOOST_FOREACH(string fn, par.fn_marker_data) {
                try {
                    Line_reader<char> lr(fn);
                }
                catch (const std::exception& e) {
                    cerr << "!!! " << fn << ": Could not open file - will be ignored" << endl;
                    failed.insert(fn);
                }
            }
            // now discard the bad files (could not be done during loop)
            BOOST_FOREACH(string s, failed) {
                par.fn_marker_data.erase(s);
            }
        }
        if (par.fn_marker_data.size() < 1) {
            cerr << "!!! No valid data file." << endl;
            return 1;
        }

        bool hasNco = vm.count("nco") > 0;
        bool hasNca = vm.count("nca") > 0;
        bool hasTraitFile = vm.count("trait-file") > 0;
        bool createTrait = (hasNco || hasNca);
        if (createTrait) {
            bool ok = true;
            if (not hasNco) {
                cerr << "!!! Option '--nco' is missing." << endl;
                ok = false;
            }
            if (not hasNca) {
                cerr << "!!! Option '--nca' is missing." << endl;
                ok = false;
            }
            if (par.ncontrol == 0 || par.ncase == 0) {
                cerr << "!!! Number of controls or cases cannot be 0." << endl;
                ok = false;
            }
            if (ok) {
                par.fn_trait = ""; //do not read trait file anymore
            }
            else if (hasTraitFile) {    //last chance
                cerr << "!!! Will try to read trait from " << par.fn_trait << "." << endl;
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
                cerr << "!!! Please specify either trait file (option '--trait-file')" <<
                    " or both number of controls ('--nco') AND cases ('--nca')." << endl; 
                return 1;
            }
        }

        if (vm.count("counts")) {
            par.pval_counts = true;
        }
        if (vm.count("haplo")) {
            par.gen_type = haplotype;
            par.tests.insert(chisq);
        }
        else {
            par.tests.insert(trend);
        }

        myout << normal <<"...User specified data file(s):" << myendl;
        BOOST_FOREACH(string fn, par.fn_marker_data) {
            myout << "\t" << fn << myendl;
        }
        myout << "...Output to: " << par.out_prefix << ".*" << myendl;
        if (createTrait) {
            myout.verbose();
            myout << "...User specified case/control numbers:" << myendl;
            myout << "...Number of controls: " << par.ncontrol << myendl;
            myout << "...Number of cases: " << par.ncase << myendl;
        }
        myout << normal << "...Number of permutations: " << par.nperm_total << myendl;
        myout << myendl;
    }
    catch(std::exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }    

    try {
        gwas_analysis(&par, myout);
    }
    catch(std::exception& e)
    {
        cout << "Error: " << e.what() << endl;
        return 1;
    }    
    myout << "...PERMORY finished successful." << myendl;
    myout << "...Runtime: " << t.elapsed() << " s" << myendl;
    return 0;
}

