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
#include "detail/matrix.hpp"
#include "detail/tokenizer.hpp"
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

}
