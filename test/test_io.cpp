// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#include <cassert>
#include <exception>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <valarray>

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/line.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/stream.hpp>

#include "config.hpp"
#include "detail/matrix.hpp"
#include "detail/tokenizer.hpp"
#include "locus.hpp"
#include "locusdata.hpp"
#include "io/format_detect.hpp"
#include "io/line_reader.hpp"
#include "io/file_out.hpp"
#include "io/read_phenotype_data.hpp"

int main(int argc, char** argv) 
{
    using namespace std;
    using namespace Permory;
    using namespace Permory::io;
    using namespace Permory::detail;
    using namespace boost::iostreams;

    const char* source = argv[1];
    const char* sink = argv[2];

    //Locus_data_reader ldr(source);

    //Out<char, std::vector<char> > lw(outf, false);
    double dd[] = {1.7, 0.6, 1.9, 7.9};
    vector<double> vv(dd, dd+4);

    std::string ss("hello world!\nhello people!");
    vector<char> vc(ss.begin(), ss.end());
    //lw.next(vv.begin(), vv.end(), " hello ");
    //lw.insert("hello world!\n");

    File_out_char lw(sink);

    //const char* source = argv[1];
    int n = atoi(argv[3]);

    PRINT(detect_marker_data_format(source));
    return 0;

    boost::shared_ptr<Locus> loc(new Locus(17, "dummy locus", "", Locus::none));

    for (int i=0; i<n; i++) {
        Line_reader<int> lr(source);
        Line_reader<int>::const_iterator it = lr.begin();
        std::string s;
        //Skip_input_filter<0> sk();
        while (!lr.eof()) {
            lr.next(s);
            PRINT(lr.char_count());
            PRINT(lr.word_count());
            PRINT(lr.line_count());
            //lr.next_unf(s);
            vector<int> v;
            v.reserve(lr.size());
            //SLIDE
            //remove_copy(lr.begin(), lr.end(), std::back_inserter(v), ' ');
            // Beagle
            remove_copy_if(lr.begin(), lr.end(), std::back_inserter(v), 
                    Skip_input_filter<2>());
            print_seq(v, "", "");
            cout << endl;
            //Locus_data<int, genotype> ld(v.begin(), v.end(), '?');
            /*
               it = lr.begin();
               while (it != lr.end()) {
               cout << *it++;
               }
            //cout << endl;
            */
        }
    }
    //PRINT(lr.word_count());
    //PRINT(lr.line_count());

    return 0;
    const char* dest = argv[2];
    double maf_thresh = atof(argv[3]);

    typedef int data_t;

    data_t myints[] = {0, 1, 2};
    set<data_t> genoset (myints, myints+3);

    //Test::perm_test(argc, argv);

}


