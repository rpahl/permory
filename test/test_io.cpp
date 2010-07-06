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
#include "helper/matrix.hpp"
#include "helper/tokenizer.hpp"
#include "locus.hpp"
#include "locusdata.hpp"
#include "io/line_reader.hpp"
#include "io/output.hpp"

int main(int argc, char** argv) 
{
    using namespace std;
    using namespace Permory;
    using namespace boost::iostreams;

    //boost::shared_ptr<Permutation> pp(new Permutation());
    
    /*
    Tokenizer toki("Hello World! -      how are   you?   ");
    while (!toki.empty()) {
        std::cerr << *toki << std::endl;
        ++toki;
    }
    return 0;
    */

    const char* source = argv[1];
    File_handler fh(source);
    Liner<char> lin(fh);
    Liner<char>::const_iterator it = lin.begin();
    while (!lin.eof()) {
        lin.next();
        while (it != lin.end()) {
            cout << *it++;
        }
        cout << endl;
    }

    return 0;
    const char* dest = argv[2];
    double maf_thresh = atof(argv[3]);

    boost::shared_ptr<Locus> loc(new Locus(17, "dummy locus", Locus::none));
    typedef int data_t;

    data_t myints[] = {0, 1, 2};
    set<data_t> genoset (myints, myints+3);

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


