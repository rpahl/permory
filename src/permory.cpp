//
// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)
//

#include <cassert>
#include <exception>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <valarray>

#include <boost/circular_buffer.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/line.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/stream.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "config.hpp"
#include "helper/matrix.hpp"
#include "helper/tokenizer.hpp"
#include "locus.hpp"
#include "locusdata.hpp"
//#include "io/line_reader.hpp"
//#include "io/output.hpp"
//#include "test/test1.hpp"

#include <fstream>
#include <iostream>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>


int main(int argc, char** argv) 
{
    //TODO: use boost::program_options for getting parameters from command line
    //TODO: use boost::serialization for example for MRU list of last x calls
    //(see Example 9 in boost::MultiIndex)
    //
    using namespace std;
    using namespace Permory;
    using namespace boost::iostreams;
    using namespace boost::program_options;

    const char* fn = argv[1];
    ifstream file(fn, ios_base::in | ios_base::binary);
    filtering_streambuf<input> in;
    in.push(gzip_decompressor());
    in.push(file);
    boost::iostreams::copy(in, cout);

    //boost::shared_ptr<Permutation> pp(new Permutation());
    
    /*
    Tokenizer toki("Hello World! -      how are   you?   ");
    while (!toki.empty()) {
        std::cerr << *toki << std::endl;
        ++toki;
    }
    return 0;
    */

    /*
    const char* source = argv[1];
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

        if (ld.maf() > maf_thresh) 
            out.write_stream(ld.begin(), ld.end());
    }
    return 0;
    */

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
