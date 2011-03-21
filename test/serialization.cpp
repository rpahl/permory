// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#define PERMORY_TEST serialization_test

#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "individual.hpp"
#include "gwas/locus.hpp"
#include "test.hpp"

using namespace std;
using namespace boost;
using namespace unit_test;
using namespace Permory;
using namespace Permory::gwas;

template<class T> void serialize_deserialize(const T& orig, T& loaded)
{
    const string filename = "serialization.test";

    {
        ofstream ofile(filename.c_str());
        boost::archive::text_oarchive oa(ofile);
        oa << orig;
    }

    {
        ifstream ifile(filename.c_str());
        boost::archive::text_iarchive ia(ifile);
        ia >> loaded;
    }
}

void record_test()
{
    Record r_orig(0.42);
    Record r_loaded;

    serialize_deserialize(r_orig, r_loaded);

    BOOST_CHECK_EQUAL(r_orig.val, r_loaded.val);
    BOOST_CHECK_EQUAL(r_orig.theType, r_loaded.theType);
}

void individual_test()
{
    const size_t id = 42;
    const string name = "Freeman";
    const Individual::Sex sex = Individual::female;
    const double c = 3.14;

    Individual ind_orig(id, name, sex, c);
    Individual ind_loaded(0, "");

    serialize_deserialize(ind_orig, ind_loaded);
    BOOST_CHECK_EQUAL(ind_orig.id(), ind_loaded.id());
    BOOST_CHECK_EQUAL(ind_orig.name(), ind_loaded.name());
    BOOST_CHECK_EQUAL(ind_orig.sex(), ind_loaded.sex());
    BOOST_CHECK_EQUAL(ind_orig.costs(), ind_loaded.costs());
}

void locus_test()
{
    Locus orig(1, "rs", "gene", Locus::chr22, 23, 42, false);
    Locus loaded;

    serialize_deserialize(orig, loaded);
    BOOST_CHECK_EQUAL(orig.id(), loaded.id());
    BOOST_CHECK_EQUAL(orig.rs(), loaded.rs());
    BOOST_CHECK_EQUAL(orig.gene(), loaded.gene());
    BOOST_CHECK_EQUAL(orig.chr(), loaded.chr());
    BOOST_CHECK_EQUAL(orig.bp(), loaded.bp());
    BOOST_CHECK_EQUAL(orig.cm(), loaded.cm());
}

test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("Functions and classes from src/detail");

    test->add(BOOST_TEST_CASE(&record_test));
    test->add(BOOST_TEST_CASE(&individual_test));
    test->add(BOOST_TEST_CASE(&locus_test));

    return test;
}

