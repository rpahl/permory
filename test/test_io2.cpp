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
#include "io/line_reader.hpp"
#include "io/output.hpp"

int main(int argc, char** argv) 
{
    using namespace std;
    using namespace Permory;
    using namespace Permory::io;
    using namespace boost::iostreams;

}


