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

    const char* dest = argv[1];
    int x = atoi(argv[2]);
    File_info fi(dest);
    PRINT(fi.file_extension(x));
    return 0;
    
    istreambuf_iterator<char> eos;
    ifstream in(argv[1]);
    istreambuf_iterator<char> iit(in);

    Output out(dest);
    int ii[] = {1, 2, 3};
    vector<int> cc(ii, ii+3); 
    //out.write_unform(cc.begin(), cc.end());
    out.write_stream(cc.begin(), cc.end());
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


