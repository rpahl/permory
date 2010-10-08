// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#define PERMORY_TEST detail_test
#include "detail/config.hpp"
#include "detail/vector.hpp"
#include "detail/matrix.hpp"
#include "test.hpp"

using namespace std;
using namespace boost::unit_test;
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
    BOOST_CHECK_EQUAL( m1.nrow(), 5 );
    BOOST_CHECK_EQUAL( m1.ncol(), 20 );

    // Data access
    // TODO
    // Modification
    // TODO
    
    // Exception
    Matrix<int> m2(0, 0);
    BOOST_CHECK_THROW (m2.row_apply(&std::valarray<int>::sum), std::out_of_range );

}

test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("Test some functions and classes from src/detail");

    test->add(BOOST_TEST_CASE(&vector_functions_test));
    test->add(BOOST_TEST_CASE(&matrix_class_test));

    return test;
}

/*
int main(int argc, char** argv)
{
    //init_unit_test_suite(argc, argv);
    //return unit_test_main( &init_unit_test_suite, argc, argv );
}
*/
