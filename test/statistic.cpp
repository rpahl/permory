// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#define PERMORY_TEST statictic_test
#include "statistical/pvalue.hpp"
#include "statistical/quantitative.hpp"
#include "test.hpp"

#include "detail/parameter.hpp"
#include "permutation/perm.hpp"
#include "gwas/locusdata.hpp"

#include <vector>
#include <utility>

using namespace std;
using namespace boost;
using namespace unit_test;
using namespace Permory::statistic;
using namespace Permory::detail;
using namespace Permory::gwas;


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


template<class T, uint L>
Locus_data<T> create_locus_data(const T (&list)[L], const Parameter& par) {
    vector<T> v(&list[0], &list[0]+L);
    return Locus_data<T>(v, par.undef_allele_code);
}

void quantitative_test() {
    typedef double T;
    const double tolerance = 0.0001;

    Parameter par;
    vector<T> trait;
    Permory::permutation::Permutation perm;

    par.useBar = true;
    typedef Quantitative<3, T> check_throw_t;
    BOOST_CHECK_THROW( check_throw_t(par, trait, &perm), invalid_argument );

    par.useBar = false;

    trait.resize(5);
    trait[0] = 0.8;
    trait[1] = 0.7;
    trait[2] = 0.5;
    trait[3] = 0.1;
    trait[4] = 0.2;

    uint marker1[5] = {0,0,1,0,2};
    uint marker2[5] = {0,0,1,0,2};
    Locus_data<uint> locus_data_1 = create_locus_data<uint>(marker1, par);
    Locus_data<uint> locus_data_2 = create_locus_data<uint>(marker2, par);

    Quantitative<3, T> q(par, trait, &perm);

    {
        vector<pair<T, T> > r = q.make_table(locus_data_1);
        BOOST_CHECK_CLOSE( r[0].first,   0.22   , tolerance );
        BOOST_CHECK_CLOSE( r[0].second,  0.3028 , tolerance );
        BOOST_CHECK_CLOSE( r[1].first,   0.04   , tolerance );
        BOOST_CHECK_CLOSE( r[1].second,  0.0016 , tolerance );
        BOOST_CHECK_CLOSE( r[2].first,  -0.26   , tolerance );
        BOOST_CHECK_CLOSE( r[2].second,  0.0676 , tolerance );
    }

    {
        vector<double> r = q.test(locus_data_1);
        BOOST_REQUIRE_EQUAL( r.size(), size_t(1) );
        BOOST_CHECK_CLOSE( r[0], 0.9530113, tolerance );
    }
}

test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("Functions and classes from src/statistical");

    test->add(BOOST_TEST_CASE(&single_step_counts_test));
    test->add(BOOST_TEST_CASE(&quantitative_test));

    return test;
}

