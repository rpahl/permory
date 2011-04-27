// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#define PERMORY_TEST statictic_test
#include "statistical/pvalue.hpp"
#include "statistical/quantitative.hpp"
#include "test.hpp"

#include "individual.hpp"
#include "detail/parameter.hpp"
#include "permutation/perm.hpp"
#include "gwas/locusdata.hpp"

#include <vector>

using namespace std;
using namespace boost;
using namespace unit_test;
using namespace Permory::statistic;
using namespace Permory::detail;
using namespace Permory::gwas;

using Permory::Individual;
using Permory::Record;


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


Individual make_individual(double phenotype) {
    Individual individual(0);
    Record record(phenotype);
    individual.add_measurement(record);
    return individual;
}

template<class T, size_t L>
Locus_data<T> create_locus_data(const T (&list)[L], const Parameter& par) {
    vector<T> v(&list[0], &list[0]+L);
    return Locus_data<T>(v, par.undef_allele_code);
}

template<class T, class Q>
deque<double> calculate_pvalue(const deque<double>& orig, Quantitative<3, T>& q,
        const Parameter &par) {

    deque<double> r(q.tmax_begin(), q.tmax_end());
    BOOST_REQUIRE_EQUAL( r.size(), par.nperm_block );

    deque<size_t> counts = single_step_counts(orig, &r);

    deque<double> pvalues = single_step_pvalues(counts, par.nperm_block);
    BOOST_REQUIRE_EQUAL( pvalues.size(), orig.size() );

    return pvalues;
}

template<class D, size_t S>
void do_quantitative_test(const D (&marker1)[S], const D (&marker2)[S]) {
    typedef double T;
    const double tolerance = 0.0001;

    Parameter par;
    par.nperm_block = 100000;
    vector<Individual> trait;
    Permory::permutation::Permutation perm;

    trait.reserve(5);
    trait.push_back(make_individual(0.8));
    trait.push_back(make_individual(0.7));
    trait.push_back(make_individual(0.5));
    trait.push_back(make_individual(0.1));
    trait.push_back(make_individual(0.2));

    par.useBar = true;
    typedef Quantitative<3, T> check_throw_t;
    BOOST_CHECK_THROW( check_throw_t(par, trait.begin(), trait.end(), &perm), invalid_argument );

    par.useBar = false;

    Locus_data<D> locus_data_1 = create_locus_data(marker1, par);
    Locus_data<D> locus_data_2 = create_locus_data(marker2, par);

    {
        Quantitative<3, T> q(par, trait.begin(), trait.end(), &perm);
        vector<Pair<T> > r = q.make_table(locus_data_1);
        BOOST_CHECK_CLOSE( r[0].first,   0.22   , tolerance );
        BOOST_CHECK_CLOSE( r[0].second,  0.3028 , tolerance );
        BOOST_CHECK_CLOSE( r[1].first,   0.04   , tolerance );
        BOOST_CHECK_CLOSE( r[1].second,  0.0016 , tolerance );
        BOOST_CHECK_CLOSE( r[2].first,  -0.26   , tolerance );
        BOOST_CHECK_CLOSE( r[2].second,  0.0676 , tolerance );
    }

    {
        Quantitative<3, T> q(par, trait.begin(), trait.end(), &perm);
        vector<double> r = q.test(locus_data_1);
        BOOST_REQUIRE_EQUAL( r.size(), size_t(1) );
        BOOST_CHECK_CLOSE( r[0], 0.9530113, tolerance );
    }
    {
        Quantitative<3, T> q(par, trait.begin(), trait.end(), &perm);
        vector<double> r = q.test(locus_data_2);
        BOOST_REQUIRE_EQUAL( r.size(), size_t(1) );
        BOOST_CHECK_CLOSE( r[0], 3.887689, tolerance );
    }

    double tolerance_permutation = 5;
    {
        Quantitative<3, T> q(par, trait.begin(), trait.end(), &perm);
        deque<double> orig;
        orig.push_back(q.test(locus_data_1)[0]);
        q.permutation_test(locus_data_1);
        deque<double> pvalues = calculate_pvalue<T, D>(orig, q, par);
        BOOST_CHECK_CLOSE( pvalues[0], 0.45, tolerance_permutation );
    }
    {
        Quantitative<3, T> q(par, trait.begin(), trait.end(), &perm);
        deque<double> orig;
        orig.push_back(q.test(locus_data_2)[0]);
        q.permutation_test(locus_data_2);
        deque<double> pvalues = calculate_pvalue<T, D>(orig, q, par);
        BOOST_CHECK_CLOSE( pvalues[0], 0.065, tolerance_permutation );
    }
    {
        Quantitative<3, T> q(par, trait.begin(), trait.end(), &perm);
        deque<double> orig;
        orig.push_back(q.test(locus_data_1)[0]);
        orig.push_back(q.test(locus_data_2)[0]);
        q.permutation_test(locus_data_1);
        q.permutation_test(locus_data_2);
        deque<double> pvalues = calculate_pvalue<T, D>(orig, q, par);
        BOOST_CHECK_CLOSE( pvalues[0], 0.63, tolerance_permutation );
        BOOST_CHECK_CLOSE( pvalues[1], 0.066, tolerance_permutation );
    }
}

void quantitative_test() {
    uint uint_marker1[5] = {0,0,1,0,2};
    uint uint_marker2[5] = {2,2,1,0,0};
    do_quantitative_test(uint_marker1, uint_marker2);

    char char_marker1[5] = {'0','0','1','0','2'};
    char char_marker2[5] = {'2','2','1','0','0'};
    do_quantitative_test(char_marker1, char_marker2);
}


test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("Functions and classes from src/statistical");

    test->add(BOOST_TEST_CASE(&single_step_counts_test));
    test->add(BOOST_TEST_CASE(&quantitative_test));

    return test;
}

