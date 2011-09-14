// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#include <iostream>
#include <iomanip>
#include <time.h>

#include "boost/lexical_cast.hpp"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/progress.hpp>   //timer

#include "detail/config.hpp"
#include "gwas/locusdata.hpp"
#include "gwas/read_locus_data.hpp"
#include "io/file.hpp"
#include "io/file_out.hpp"
#include "io/line_reader.hpp"
#include "io/output.hpp"


//
// Compact data files of various formats by discarding space characters
void file_compact(const std::string& fn_in, const std::string& fn_out, 
        char undef, size_t nc=0)
{
    using namespace std;
    using namespace boost;
    using namespace Permory::detail;
    using namespace Permory::gwas;
    using namespace Permory::io;
    cout << stdpre << "Scanning file..." << endl;
    Line_reader<string> lr(fn_in);
    string line;
    size_t n = 0;
    while (not lr.eof()) {
        lr.next_str(line);
        n++;
    }
    n--;
    Datafile_format dff = detect_marker_data_format(fn_in, undef);
    if (dff == presto) {
        n--;
    }
    cout << stdpre << n << " lines found" << endl;
    cout << stdpre << "Writing `" << fn_out << "'..." << endl;
    progress_display progress(n);

    File_out_char fout(fn_out);
    Locus_data_reader<char> loc_reader(fn_in, undef);
    while (loc_reader.hasData()) {
        std::vector<char> v;
        size_t nskipped = loc_reader.get_next(v);
        if (v.empty()) {
            continue;
        }
        if (nc == 0) {
            fout(v.begin(), v.end());
        }
        else {
            Locus_data<char> data(v, undef);
            data = data.condense_alleles_to_genotypes(nc);
            fout(data.begin(), data.end());
        }
        fout.endl();
        for (size_t i=0; i<nskipped+1; ++i) {
            ++progress;
        }
    }
    cout << endl;
}


int main(int ac, char* av[])
{
    namespace cls = boost::program_options::command_line_style;
    using namespace std;
    using namespace boost;
    using namespace boost::program_options;
    using namespace Permory;
    using namespace Permory::detail;
    using namespace Permory::io;

    timer t;            //t.restart(); //start clock
    char undef;         //allele code for undefined data
    size_t ncon=0;      //number of alleles to condense
    string fn_in;       //data file name
    string fn_out;      //out file name

    cout << endl;
    cout << "+------------------------------------------------------+" << endl;
    cout << "|             PERMORY data compactor, v1.0             |" << endl;
    cout << "+------------------------------------------------------+" << endl;
    cout << "|            Copyright (c) 2010 Roman Pahl             |" << endl;
    cout << "|  Distributed under the Boost Software License, v1.0  |" << endl;
    cout << "+------------------------------------------------------+" << endl;
    cout << "|                   www.permory.org                    |" << endl;
    cout << "+------------------------------------------------------+" << endl;
    cout << endl;

    try {
        // Command line options
        options_description cmd_line_options("Option");
        cmd_line_options.add_options()
            ("help,h", "show this helps")
            ("file,f", value<string>(&fn_in), "file to convert")
            ("condense,c",  value<size_t>(&ncon),
             "condense 'arg' alleles to genotypes (maximal 9)")
            ("missing,m",  value<char>(&undef)->default_value('?'),
             "code of missing marker data (single character)")
            ("out,o", value<string>(&fn_out), "output file name (optional)")
            ;

        positional_options_description pos;
        pos.add("file", -1);

        // parse command line
        variables_map vm;
        store(command_line_parser(ac, av).
                options(cmd_line_options).
                positional(pos).
                run(), vm);

        if (vm.count("help")) {
            cout << "Converts marker data files into compact format (*.comp)." <<
            endl << "Supported input formats (automatic detection):" << endl <<
            "\tPLINK transposed fileset (*.tped)" << endl << 
            "\tPRESTO (*.bgl)" << endl <<
            "\tSLIDE (*.slide)" << endl <<
            "Usage:" << endl << "\tconverter <file_to_convert>" << endl << 
            endl << cmd_line_options << endl;
            return 0;
        }

        bool hasFile = vm.count("file") > 0;
        if (not hasFile) {
            cerr << errpre << "Please specify a file." << endl;
            cerr << errpre << "For more information try ./converter --help" << endl;
            return 1;
        }
        notify(vm);

        if (ncon > 9) {
            cerr << errpre << "Maximal 9 alleles allowed for condension." << endl;
            return 1;
        }
        File_handle fin(fn_in);
        if (not bfs::exists(*fin)) {
            cerr << errpre << "`" << fn_in << "': No such file." << endl;
            return 1;
        }
        if (fn_out.empty()) {
            // Derive fn_out file name from input file
            size_t found = fn_in.find_last_of(".");
            fn_out = fn_in.substr(0,found);
            std::string ext = fin.extension(0);
            if (ext == "gz") {
                found = fn_out.find_last_of(".");
                fn_out = fn_out.substr(0,found);
            }
            fn_out.append(".comp");
            if (ext == "gz") {
                fn_out.append(".gz");
            }
        }
        Datafile_format dff = detect_marker_data_format(fn_in, undef);
        if (dff == unknown) {
            cerr << errpre << "Unknown data format." << endl;
        }
        cout << stdpre << "Input file: " << "`" << fn_in << 
            "' -> assuming file format " << 
            Permory::detail::datafile_format_to_string(dff) << endl;
    }
    catch(std::exception& e)
    {
        cerr << e.what() << endl;
        return 1;
    }    

    try {
        file_compact(fn_in, fn_out, undef, ncon);
    }
    catch(std::exception& e)
    {
        cerr << errpre << "Error: " << e.what() << endl;
        return 1;
    }    
    cout << stdpre << "Compaction complete." << endl;
    cout << stdpre << "Runtime: " << t.elapsed() << " s" << endl;
    return 0;
}

