// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#define PERMORY_TEST statictic_test
#include "statistical/pvalue.hpp"
#include "test.hpp"

using namespace std;
using namespace boost;
using namespace unit_test;
using namespace Permory::statistic;


void single_step_counts_test()
{
    deque<double> l1;
    l1.push_back(0.1);
    l1.push_back(0.5);
    l1.push_back(0.9);
    l1.push_back(1.0);

    deque<double> l2;
    l2.push_back(0.1);
    l2.push_back(0.2);
    l2.push_back(0.3);
    l2.push_back(0.4);
    l2.push_back(0.5);
    l2.push_back(0.6);
    l2.push_back(0.7);
    l2.push_back(0.8);
    l2.push_back(0.9);
    l2.push_back(1.0);

    deque<size_t> result;
 
    result = single_step_counts(l1, &l2);

    BOOST_CHECK_EQUAL( result.at(0), size_t(10) );
    BOOST_CHECK_EQUAL( result.at(1), size_t(6) );
    BOOST_CHECK_EQUAL( result.at(2), size_t(2) );
    BOOST_CHECK_EQUAL( result.at(3), size_t(1) );

    deque<double> l3;
    l3.push_back(0.3);
    l3.push_back(0.2);
    l3.push_back(0.1);
    l3.push_back(0.5);
    l3.push_back(0.4);
    l3.push_back(0.6);
    l3.push_back(0.7);
    l3.push_back(0.8);
    l3.push_back(0.9);

    result = single_step_counts(l1, &l3);

    BOOST_CHECK_EQUAL( result.at(0), size_t(9) );
    BOOST_CHECK_EQUAL( result.at(1), size_t(5) );
    BOOST_CHECK_EQUAL( result.at(2), size_t(1) );
    BOOST_CHECK_EQUAL( result.at(3), size_t(0) );
}

test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("Functions and classes from src/statistical");

    test->add(BOOST_TEST_CASE(&single_step_counts_test));

    return test;
}

