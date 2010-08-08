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
#include "detail/format_detect.hpp"
#include "detail/matrix.hpp"
#include "detail/tokenizer.hpp"
#include "locus.hpp"
#include "locusdata.hpp"
#include "io/line_reader.hpp"
#include "io/output.hpp"

int main(int argc, char** argv) 
{
    using namespace std;
    using namespace Permory;
    using namespace Permory::io;
    using namespace Permory::detail;
    using namespace boost::iostreams;

    const char* source = argv[1];
    const char* sink = argv[2];

    //Line_writer<char, std::vector<char> > lw(outf, false);
    double dd[] = {1.7, 0.6, 1.9, 7.9};
    vector<double> vv(dd, dd+4);

    std::string ss("hello world!\nhello people!");
    vector<char> vc(ss.begin(), ss.end());
    //lw.next(vv.begin(), vv.end(), " hello ");
    //lw.insert("hello world!\n");

    File_handle outf(sink);
    Line_writer<char> lw(outf, true);

    //const char* source = argv[1];
    int n = atoi(argv[3]);
    File_handle fh(source);

    PRINT(detect_data_format(fh));
    return 0;

    boost::shared_ptr<Locus> loc(new Locus(17, "dummy locus", Locus::none));

    for (int i=0; i<n; i++) {
        Line_reader<int> lr(fh);
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
            //Locus_data<char> ld(v.begin(), v.end(), 9, loc, Locus_data<char>::genotype);
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

    //boost::shared_ptr<Locus> loc(new Locus(17, "dummy locus", Locus::none));
    typedef int data_t;

    data_t myints[] = {0, 1, 2};
    set<data_t> genoset (myints, myints+3);

    /*
       Tokenizer tok("");
       Output out(dest);
       for (Line_reader lr(source); lr.good(); ++lr) {
       tok.assign(*lr);
       Locus_data<int> ld(tok, 9, loc, Locus_data<int>::genotype);
       ld.add_to_domain(genoset);
    //Locus_data<char> ld(lr.begin(), lr.end(), 9, loc);
    //ld.print();

    //if (ld.maf() > maf_thresh) out.write_stream(ld.begin(), ld.end());
    }
    return 0;

    /*
    io::filtering_ostream teeOut;
    teeOut.push( io::tee(io::file_sink("out1.txt")) );
    //teeOut.push( io::tee(io::file_sink("out1.txt")) );
    teeOut.push( io::tee(io::file_sink("out2.txt")) );
    //teeOut.push( std::cout );

    teeOut << "Hello World!\n"
    << 23 << '*' << 45 << '=' << 23*45
    << std::endl;
    //teeOut.pop();

    teeOut << "Another chunk of output!" << std::endl;      
    */
    return 0;
    //Test::perm_test(argc, argv);

    /*
       std::map<int, int> m(10);
       std::copy(m.begin(), m.end(); std:back_inserter(unic));
       */
    //int a = atoi(argv[1]);
    //int b = atoi(argv[2]);

    return 0;
}


