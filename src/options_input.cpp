// Copyright (c) 2011 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#include <fstream>


#include "individual.hpp"
#include "detail/config.hpp"
#include "detail/enums.hpp"
#include "detail/exception.hpp"
#include "detail/parameter.hpp"
#include "detail/program_options.hpp"
#include "io/file.hpp"
#include "io/output.hpp"

namespace Permory { 

    boost::program_options::variables_map get_options(int ac, char* av[])
    {
        using namespace std;
        using namespace boost;
        using namespace boost::program_options;
        using namespace detail;

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
            ;
        //
        // Analysis
        //
        options_description analysis("Analysis");
        analysis.add_options()
            ("max-maf", my_value<double>("NUM")->my_default_value(0.5),
             "upper minor allele frequency threshold")
            ("min-maf", my_value<double>("NUM")->my_default_value(0.0),
             "lower minor allele frequency threshold") 
            ("nca", my_value<size_t>("NUM")->my_default_value(0,""), 
             "number of cases (requires option '--nco')\n" 
             "Assumes: data format [...controls...|...cases...]")
            ("nco", my_value<size_t>("NUM")->my_default_value(0,""), 
             "number of controls (requires option '--nca')")
            ("phenotype",
             //my_value<Record::Value_type>("NUM")->my_default_value(Record::dichotomous, " (=1)"),
             my_value<size_t>("NUM")->my_default_value(1),
             "1=dichotomous, 2=continuous")
            ;
        //
        // Data
        //
        options_description data("Data");
        data.add_options()
            ("allelic", "perform allelic-based permutation") 
            ("missing", my_value<string>("CHAR")->my_default_value("?"),
             "code of missing marker data")
            ;
        //
        // I/O
        //
        options_description io("Input/Output");
        io.add_options()
            ("config,c", my_value<string>("FILE"),
             "include program options from FILE")
            ("log", my_value<string>("FILE")->my_default_value("", ""), 
             "write console output into FILE")
            ("out-prefix,o", my_value<string>("STRING")->my_default_value("out"),
             "prefix for all output files") 
            ("trait-file,f", my_value<string>("FILE")->my_default_value("", ""), 
             "read trait data from FILE")
            //("stat,s", my_value<string>(&tests, "STRING"), "statistical tests to be used")
            ;
        //
        // Permutation
        //
        options_description perm("Permutation");
        perm.add_options()
            ("nperm,n", my_value<size_t>("NUM")->my_default_value(10000),
             "Number of permutations")
            ("seed", my_value<int>("NUM")->my_default_value(12345678), 
             "random seed")
            ;
        //
        // Advanced
        //
        // Statistical testing
        // speed optimization
        options_description advanced("Advanced options");
        advanced.add_options()
            ("alpha", my_value<double>("NUM")->my_default_value(0.05), 
             "genome-wide significance threshold determining the effective number of tests")
            ("block",my_value<size_t>("NUM")->my_default_value(10000),  
             "permutation block size")
            ("counts", "in addition to p-values, output #(T_perm > T_orig)")
            ("debug,d", "most detailed output")
            ("ntop", my_value<size_t>("NUM")->my_default_value(100), 
             "number of top markers listed in *.top output file")
            ("tail", my_value<size_t>("NUM")->my_default_value(100), 
             "size of sliding tail (REM method)")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        options_description hidden("Hidden");
        hidden.add_options()
            ("data-file", value<vector<string> >()->composing(), 
             "input data file")
            ;

        // Create option "profiles"
        options_description cmd_line_options;
        cmd_line_options.add(general).add(analysis).add(data).add(io).add(perm).
            add(advanced).add(hidden);

        options_description config_file_options;
        config_file_options.add(general).add(analysis).add(data).add(io).add(perm).
            add(advanced).add(hidden);

        options_description visible("Options");
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
            exit(0);
        }
        if (vm.count("help-all")) {
            cout << "Usage:\n\tpermory [options] <data_file1> [data_file2 ...]\n\n"; 
            cout << visible << endl;
            cout << advanced << endl;
            exit(0);
        }
        if (vm.count("help-adv")) {
            cout << advanced << endl;
            exit(0);
        }
        notify(vm);

        bool hasConfigFile = vm.count("config") > 0;
        if (hasConfigFile) {
            string config_file = vm["config"].as<string>();
            ifstream ifs(config_file.c_str());
            if (!ifs) {
                string msg("Unable to open config file: ");
                msg.append(config_file);
                throw runtime_error(msg);
            }
            else {
                // in any case, config file is shadowed by cmd_line inputs
                store(parse_config_file(ifs, config_file_options), vm);
                notify(vm);
            }
        }
        return vm;
    }

    void check_options(
            Permory::io::Myout& myout,
            const boost::program_options::variables_map& vm) 
    {
        using namespace std;
        using namespace Permory::io;
        using namespace Permory::detail;

        // Missing options 
        if (vm.count("data-file") == 0) {
            throw Missing_option("data file"); 
        }

        vector<string> file_names = vm["data-file"].as<vector<string> >();
        bool hasValidFile = false;
        BOOST_FOREACH(string fn, file_names) {
            ifstream ifs(fn.c_str());
            if (ifs) {
                hasValidFile = true;
            }
            else {
                myout << warnpre << fn << ": Could not open file - will be ignored" << endl;
            }
        }
        if (not hasValidFile) {
            throw runtime_error("No valid data file.");
        }

        bool hasNca = not vm["nca"].defaulted();
        bool hasNco = not vm["nco"].defaulted();
        bool hasTraitFile = not vm["trait-file"].defaulted();
        if (hasTraitFile) {
            string fn = vm["trait-file"].as<string>();
            ifstream ifs(fn.c_str());
            if (!ifs) {
                throw runtime_error("Unable to open trait file: " + fn);
            }
        }
        else {
            if (hasNca && !hasNco) {
                throw Missing_option("number of controls (option --nco)"); 
            }
            if (!hasNca && hasNco) {
                throw Missing_option("number of cases (option --nca)"); 
            }
            if (!hasNca && !hasNco) {
                throw Missing_option("trait file (option --trait-file) OR "
                        "number of cases and controls (--nca and --nco)"); 
            }
        }

        // Bad arguments
        if (!hasTraitFile && hasNca && hasNco) {
            if (vm["nca"].as<size_t>() == 0 || vm["nco"].as<size_t>() == 0) {
                throw invalid_argument("number of cases/controls must not be 0");
            }
        }
        if (vm["phenotype"].as<size_t>() < 1 && vm["phenotype"].as<size_t>() > 2) {
            throw invalid_argument("invalid phenotype.");
        }
        double minmaf = vm["min-maf"].as<double>();
        double maxmaf = vm["max-maf"].as<double>();
        if (minmaf < 0) {
            throw invalid_argument("--min-maf must not be < 0");
        }
        if (maxmaf > 1) {
            throw invalid_argument("--max-maf must not be > 1");
        }
        if (minmaf > maxmaf) {
            throw invalid_argument("--min-maf must be < than --max-maf");
        }
        if (vm["alpha"].as<double>() <= 0 || vm["alpha"].as<double>() > 1) {
            throw invalid_argument("significance threshold --alpha must be in [0,1].");
        }

        // Obsolete options
        if (hasTraitFile && hasNca) {
            myout << normal << warnpre << "Ignoring option --nca, because " 
                "option --trait-file is used." << endl; 
        }
        if (hasTraitFile && hasNco) {
            myout << normal << warnpre << "Ignoring option --nco, because " 
                "option --trait-file is used." << endl; 
        }
        string s = vm["missing"].as<string>();
        if (s.size() > 1) {
            myout << normal << warnpre << "Only using first character '" <<
            s.at(0) << "' from option --missing = " << s << endl; 
        }
    }

    void set_parameter(
            detail::Parameter& par, 
            const boost::program_options::variables_map& vm) 
    {
        using namespace std;
        using namespace boost;
        using namespace Permory::io;
        using namespace Permory::detail;

        // General
        par.interactive = vm.count("interactive") > 0;
        par.quiet = vm.count("quiet") > 0;
        par.verbose = vm.count("verbose") > 0;

        // Analysis
        par.max_maf = vm["max-maf"].as<double>();
        par.min_maf = vm["min-maf"].as<double>();
        par.ncase = vm["nca"].as<size_t>();
        par.ncontrol = vm["nco"].as<size_t>();
        switch (vm["phenotype"].as<size_t>()) {
            case 1:
                par.phenotype_domain = Record::dichotomous;
                break;
            case 2:
                par.phenotype_domain = Record::continuous;
                break;
            default:
                throw invalid_argument("invalid phenotype.");
        }
        par.useBar = par.phenotype_domain == Record::dichotomous;

        // Data
        if (vm.count("allelic")) {
            par.marker_type = allelic;
            par.tests.insert(chisq);
        }
        else {
            par.tests.insert(trend);
        }
        par.undef_allele_code = vm["missing"].as<string>().at(0);

        // I/O
        vector<string> file_names = vm["data-file"].as<vector<string> >();
        BOOST_FOREACH(string fn, file_names) {
            ifstream ifs(fn.c_str());
            if (ifs) {
                par.fn_marker_data.insert(fn);
            }
        }
        par.log_file = vm["log"].as<string>();
        par.out_prefix = vm["out-prefix"].as<string>();
        par.fn_trait = vm["trait-file"].as<string>();


        // Permutation
        par.nperm_total = vm["nperm"].as<size_t>();
        par.seed = vm["seed"].as<int>();

        // Advanced
        par.alpha = vm["alpha"].as<double>();
        par.nperm_block = vm["block"].as<size_t>();
        par.pval_counts = vm.count("counts") > 0;
        par.debug = vm.count("debug") > 0;
        par.ntop = vm["ntop"].as<size_t>();
        par.tail_size = vm["tail"].as<size_t>();
    }
}   //namespace Permory

