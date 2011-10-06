// Copyright (c) 2010-2011 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#include <iostream>
#include <iomanip>
#include <time.h>

#include <boost/progress.hpp>   //timer

#include "options_input.cpp"

#include "detail/config.hpp"
#include "detail/enums.hpp"
#include "detail/parameter.hpp"
#include "detail/program_options.hpp"
#include "gwas/analysis.hpp"
#include "io/file.hpp"
#include "io/line_reader.hpp"
#include "io/output.hpp"

//#define BOOST_RESULT_OF_USE_DECLTYPE
#include <boost/utility/result_of.hpp>

namespace Permory { 

    void print_head(io::Myout& myout, double version) 
    {
        using namespace std;
        myout << endl;
        myout << "+-----------------+-----------------+------------------+" << endl;
        myout << "|    PERMORY      |      v"<< fixed << setprecision(2) << 
            version << "      |    6/Oct/2011    |" << endl;
        myout << "+-----------------+-----------------+------------------+" << endl;
        myout << "|          Copyright (c) 2010-2011 Roman Pahl          |" << endl;
        myout << "|                             2011 Volker Steiß        |" << endl;
        myout << "|  Distributed under the Boost Software License, v1.0  |" << endl;
        myout << "+------------------------------------------------------+" << endl;
        myout << "|                   www.permory.org                    |" << endl;
        myout << "+------------------------------------------------------+" << endl;
    }

    void print_options_in_effect(io::Myout& myout, 
            const boost::program_options::variables_map& vm)
    {
        using namespace std;
        using namespace Permory::io;
        boost::program_options::variables_map::const_iterator itMap;

        myout << stdpre << "Options set by user:" << endl;
        for (itMap = vm.begin(); itMap != vm.end(); ++itMap) {
            if (not (itMap->second.defaulted())) {
                string option = itMap->first;
                if (option != "data-file") {//skip vector<string> of file names
                    myout << "   --" << option;

                    // Workaround: since we dont know the type, we play try and catch :)
                    try { myout  << " = " << itMap->second.as<double>(); }
                    catch (...) {
                        try { myout  << " = " << itMap->second.as<size_t>(); }
                        catch (...) {
                            try { 
                                string s = itMap->second.as<string>();
                                if (not s.empty()) {
                                    myout  << " = " << s;
                                }
                            }
                            catch (...) {
                                try {myout  << "4= " << itMap->second.as<int>(); }
                                catch (...) { }
                            }
                        }
                    }
                    myout << endl;
                }
            }
        }
    }

    void set_program_verbosity(detail::Parameter& par, io::Myout& myout)
    {
        using namespace std;
        using namespace Permory::detail;
        if (par.quiet) {
            myout.set_verbosity(muted);
        }
        else if (par.debug) {
            myout.set_verbosity(all);
        }
        else if (par.verbose) {
            myout.set_verbosity(verbose);
        }
    }

}   //namespace permory

int main(int ac, char* av[])
{
    namespace cls = boost::program_options::command_line_style;
    using namespace std;
    using namespace Permory;
    using namespace Permory::detail;
    using namespace Permory::gwas;
    using namespace Permory::io;

    Parameter par;
    Myout myout(&par); 
    Permory::print_head(myout, par.version);

    // Activation/usage of MPI is handled by a static hook
    Permory::hook::Argument_hook()(&ac, &av);    

    try {
        boost::program_options::variables_map vm = Permory::get_options(ac, av);

        string logfn = vm["log"].as<string>();
        if (not logfn.empty()) {
            myout.set_logfile(logfn, par.interactive);
            myout << normal << stdpre << "Writing this to log file `" << logfn <<"'." << endl;
        }

        Permory::check_options(myout, vm);
        Permory::set_parameter(par, vm);
        Permory::set_program_verbosity(par, myout);
        Permory::print_options_in_effect(myout, vm);
        myout << normal << indent(3) << "data file(s):" << endl;
        BOOST_FOREACH(string fn, par.fn_marker_data) {
            myout << indent(6) << fn << endl;
        }
        myout << stdpre<< "Output to: " << par.out_prefix << ".*" << endl;
        myout << normal << stdpre << "Number of permutations: " << 
            par.nperm_total << endl << endl;

        boost::timer t;                
        time_t rawtime; 
        struct tm * timeinfo;   //time and date

        time(&rawtime);
        timeinfo = localtime(&rawtime);
        myout << normal << stdpre << "Started at " << asctime(timeinfo) << endl;

        // Start main analysis 
        t.restart();    //start clock
        analyzer_factory_t factory;
        gwas_analysis(&par, myout, factory);    //see src/gwas/analysis.hpp

        time(&rawtime);
        timeinfo = localtime(&rawtime);
        myout << stdpre << "Finished at " << asctime(timeinfo);
        myout << stdpre << "User runtime: " << t.elapsed() << " s" << endl;
    }
    catch(const Ambigous_option& e) {
        myout << errpre << "Ambigous option: " << e.what() << endl;
    }    
    catch(const Missing_option& e) {
        myout << errpre << "Missing option - please specify " << e.what() << endl;
        myout << errpre << "For more information try ./omd --help" << endl;
    }    
    catch(const Data_length_mismatch_error& e) {
        myout << errpre << "Data length mismatch: " << e.what() << endl;
    }
    catch(const std::exception& e) {
        myout << errpre << "Error: " << e.what() << endl;
    }    
}

