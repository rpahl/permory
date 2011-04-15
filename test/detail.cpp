// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#define PERMORY_TEST detail_test
#include "detail/config.hpp"
#include "detail/vector.hpp"
#include "detail/matrix.hpp"
#include "detail/functors.hpp"
#include "detail/pair.hpp"
#include "test.hpp"

#include <utility>

using namespace std;
using namespace boost;
using namespace unit_test;
using namespace Permory::detail;


void vector_functions_test()
{
    int a1[4] =  {1,2,3,4};
    int a2[4] =  {1,7,0,6};

    vector<int> v(a1, a1+4);
    BOOST_CHECK_EQUAL( sequence_is_sorted(v), true );
    BOOST_CHECK_EQUAL( sequence_is_sorted(v, true), false );

    v.assign(a2, a2+4);
    BOOST_CHECK_EQUAL( sequence_is_sorted(v), false );
    BOOST_CHECK_EQUAL( sequence_is_sorted(v, true), false );
}


void matrix_class_test()
{
    // Constructor test
    Matrix<int> m1(5, 20);
    BOOST_CHECK_EQUAL( m1.empty(), false );
    BOOST_CHECK_EQUAL( m1.nrow(), size_t(5) );
    BOOST_CHECK_EQUAL( m1.ncol(), size_t(20) );

    // Data access
    // TODO
    // Modification
    // TODO

    // Exception
    Matrix<int> m2(0, 0);
    BOOST_CHECK_THROW (m2.row_apply(&std::valarray<int>::sum), std::out_of_range );

}


void deque_concat_test()
{
    struct deque_concat<double> concatenator;
    deque<double> l1(4, 1.);
    deque<double> l2(4, 2.);

    BOOST_CHECK_EQUAL((concatenator(l1, l2)).size(), size_t(8));
    BOOST_CHECK_EQUAL((concatenator(l2, l1)).size(), size_t(8));

    const deque<double>& q1 = concatenator(l1, l2);
    BOOST_CHECK_EQUAL(q1.at(0), 1.);
    BOOST_CHECK_EQUAL(q1.at(1), 1.);
    BOOST_CHECK_EQUAL(q1.at(2), 1.);
    BOOST_CHECK_EQUAL(q1.at(3), 1.);
    BOOST_CHECK_EQUAL(q1.at(4), 2.);
    BOOST_CHECK_EQUAL(q1.at(5), 2.);
    BOOST_CHECK_EQUAL(q1.at(6), 2.);
    BOOST_CHECK_EQUAL(q1.at(7), 2.);


    deque<double> l3;
    l3.push_back(1.);
    l3.push_back(2.);
    l3.push_back(3.);

    const deque<double>& q2 = concatenator(l3, l3);
    BOOST_CHECK_EQUAL(q2.at(0), 1.);
    BOOST_CHECK_EQUAL(q2.at(1), 1.);
    BOOST_CHECK_EQUAL(q2.at(2), 2.);
    BOOST_CHECK_EQUAL(q2.at(3), 2.);
    BOOST_CHECK_EQUAL(q2.at(4), 3.);
    BOOST_CHECK_EQUAL(q2.at(5), 3.);
}


void pair_helper_test() {
    typedef int T;
    typedef pair<T, T> P;

    P a(1,2);
    P b(1,2);
    P c;

    c = a + b;
    BOOST_CHECK_EQUAL( c.first,  2 );
    BOOST_CHECK_EQUAL( c.second, 4 );

    c += a;
    BOOST_CHECK_EQUAL( c.first,  3 );
    BOOST_CHECK_EQUAL( c.second, 6 );
    BOOST_CHECK_EQUAL( a.first,  1 );
    BOOST_CHECK_EQUAL( a.second, 2 );
}


test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("Functions and classes from src/detail");

    test->add(BOOST_TEST_CASE(&vector_functions_test));
    test->add(BOOST_TEST_CASE(&matrix_class_test));

    test->add(BOOST_TEST_CASE(&deque_concat_test));

    test->add(BOOST_TEST_CASE(&pair_helper_test));

    return test;
}

