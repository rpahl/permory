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

template<class T> set<T> create_domain();
template<> set<uint> create_domain() {
    set<uint> result;
    result.insert(0);
    result.insert(1);
    result.insert(2);
    return result;
}
template<> set<char> create_domain() {
    set<char> result;
    result.insert('0');
    result.insert('1');
    result.insert('2');
    return result;
}

template<class T, size_t L>
Locus_data<T> create_locus_data(const T (&list)[L], const Parameter& par) {
    vector<T> v(&list[0], &list[0]+L);
    Locus_data<T> locus_data(v, par.undef_allele_code);
    locus_data.add_to_domain(create_domain<T>());
    return locus_data;
}

template<class T, size_t S, size_t L>
vector<Locus_data<T> > create_locus_datas(const T (&list)[S][L], const Parameter& par) {
    vector<Locus_data<T> > result;
    for (size_t i = 0; i < S; ++i) {
        result.push_back(create_locus_data(list[i], par));
    }
    return result;
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

template<class T> vector<T> create_trait();
template<> vector<Individual> create_trait() {
    vector<Individual> result;
    result.reserve(5);
    result.push_back(make_individual(0.8));
    result.push_back(make_individual(0.7));
    result.push_back(make_individual(0.5));
    result.push_back(make_individual(0.1));
    result.push_back(make_individual(0.2));
    return result;
}
template<> vector<double> create_trait() {
    vector<double> result;
    result.reserve(5);
    result.push_back(0.8);
    result.push_back(0.7);
    result.push_back(0.5);
    result.push_back(0.1);
    result.push_back(0.2);
    return result;
}
template<> vector<Pair<double> > create_trait() {
    vector<Pair<double> > result;
    result.reserve(5);
    result.push_back(make_pair( 0.34, 0.1156));
    result.push_back(make_pair( 0.24, 0.0576));
    result.push_back(make_pair( 0.04, 0.0016));
    result.push_back(make_pair(-0.36, 0.1296));
    result.push_back(make_pair(-0.26, 0.0676));
    return result;
}

template<class T, class D, uint K>
deque<T> create_orig(const vector<Locus_data<D> > data, Quantitative<K, T> &q) {
    deque<T> result;
    BOOST_FOREACH(Locus_data<D> locus_data, data) {
        result.push_back(q.test(locus_data)[0]);
    }
    return result;
}

template<class D, size_t S>
void do_quantitative_test(const D (&marker1)[S], const D (&marker2)[S]) {
    typedef double T;
    const double tolerance = 0.0001;

    Parameter par;
    par.nperm_block = 100000;
    vector<Individual> trait(create_trait<Individual>());
    vector<Pair<T> > trait_pair = create_trait<Pair<T> >();
    Permory::permutation::Permutation perm;

    par.useBar = true;
    typedef Quantitative<3, T> check_throw_t;
    BOOST_CHECK_THROW( check_throw_t(par, trait.begin(), trait.end(), &perm), invalid_argument );

    par.useBar = false;

    Locus_data<D> locus_data_1 = create_locus_data(marker1, par);
    Locus_data<D> locus_data_2 = create_locus_data(marker2, par);

    {
        Quantitative<3, T> q(par, trait.begin(), trait.end(), &perm);
        vector<Pair<T> > r;
        q.make_table<D>(trait_pair.begin(), trait_pair.end(),
                        locus_data_1.begin(), locus_data_1.end(),
                        r);
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

    { // Test with boosters.
        typedef double T;
        typedef char D;

        Parameter par;
        par.useBar = false;
        par.nperm_block = 100000;
        const double tolerance = 0.0001;
        const double tolerance_permutation = 5;

        char char_markers[6][5] = {
                {'2','2','1','0','0'},
                {'2','2','1','0','1'},
                {'2','1','0','0','0'},
                {'0','0','1','0','2'},
                {'0','0','0','0','2'},
                {'0','0','1','1','2'}
            };

        vector<Individual> trait(create_trait<Individual>());
        vector<Locus_data<D> > data = create_locus_datas(char_markers, par);
        Permory::permutation::Permutation perm(134687313);
        Quantitative<3, T> q(par, trait.begin(), trait.end(), &perm);

        deque<T> orig(create_orig(data, q));
        BOOST_CHECK_CLOSE( orig[0], 3.8876890, tolerance );
        BOOST_CHECK_CLOSE( orig[1], 2.9429790, tolerance );
        BOOST_CHECK_CLOSE( orig[2], 2.7537741, tolerance );
        BOOST_CHECK_CLOSE( orig[3], 0.9530113, tolerance );
        BOOST_CHECK_CLOSE( orig[4], 1.2193362, tolerance );
        BOOST_CHECK_CLOSE( orig[5], 3.3058471, tolerance );

        BOOST_FOREACH(Locus_data<D> l, data) {
            q.permutation_test(l);
        }
        deque<double> pvalues = calculate_pvalue<T, D>(orig, q, par);

        BOOST_CHECK_CLOSE( pvalues[0], 0.066, tolerance_permutation );
        BOOST_CHECK_CLOSE( pvalues[1], 0.332, tolerance_permutation );
        BOOST_CHECK_CLOSE( pvalues[2], 0.398, tolerance_permutation );
        BOOST_CHECK_CLOSE( pvalues[3], 0.899, tolerance_permutation );
        BOOST_CHECK_CLOSE( pvalues[4], 0.833, tolerance_permutation );
        BOOST_CHECK_CLOSE( pvalues[5], 0.233, tolerance_permutation );
    }
}


void teststat_test() {
    //
    // Dichotom
    {
        Con_tab<2,3> contab;
        contab[0][0] = 0;
        contab[0][1] = 0;
        contab[0][2] = 0;
        contab[1][0] = 0;
        contab[1][1] = 0;
        contab[1][2] = 0;

        Trend tt;
        BOOST_CHECK_EQUAL( tt(contab), 0 );
    }
    {
        Con_tab<2,3> contab;
        contab[0][0] = 0;
        contab[0][1] = 0;
        contab[0][2] = 100;
        contab[1][0] = 100;
        contab[1][1] = 0;
        contab[1][2] = 0;

        Trend tt;
        BOOST_CHECK_EQUAL( tt(contab), 200 );
    }

    //
    // Quantitative
    {
        const double tolerance = 0.0001;
        typedef double T;
        typedef Pair<T> P;

        Parameter par;
        par.val_type = Record::continuous;
        par.useBar = false;

        vector<double> trait(create_trait<double>());

        Trend_continuous<T> test(trait);

        BOOST_CHECK_CLOSE( test.get_mu_y(), 0.46, tolerance );

        const vector<P> buffer(test.get_buffer());
        BOOST_CHECK_CLOSE( buffer[0].first,  0.34,   tolerance );
        BOOST_CHECK_CLOSE( buffer[0].second, 0.1156, tolerance );
        BOOST_CHECK_CLOSE( buffer[1].first,  0.24,   tolerance );
        BOOST_CHECK_CLOSE( buffer[1].second, 0.0576, tolerance );
        BOOST_CHECK_CLOSE( buffer[2].first,  0.04,   tolerance );
        BOOST_CHECK_CLOSE( buffer[2].second, 0.0016, tolerance );
        BOOST_CHECK_CLOSE( buffer[3].first, -0.36,   tolerance );
        BOOST_CHECK_CLOSE( buffer[3].second, 0.1296, tolerance );
        BOOST_CHECK_CLOSE( buffer[4].first, -0.26,   tolerance );
        BOOST_CHECK_CLOSE( buffer[4].second, 0.0676, tolerance );

        char marker[5] = {'0','0','1','0','2'};
        Locus_data<char> locus_data = create_locus_data(marker, par);

        test.update(locus_data);

        BOOST_CHECK_CLOSE( test.get_mu_j(), 0.6, tolerance );
        BOOST_CHECK_CLOSE( test.get_denom_invariant(), 0.13392, tolerance );

        vector<P> sums(3);
        sums[0] = buffer[0] + buffer[1] + buffer[3];
        sums[1] = buffer[2];
        sums[2] = buffer[4];
        T result = test(sums);
        BOOST_CHECK_CLOSE( result, 0.9530113, tolerance );
    }
}


test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("Functions and classes from src/statistical");

    test->add(BOOST_TEST_CASE(&single_step_counts_test));
    test->add(BOOST_TEST_CASE(&quantitative_test));
    test->add(BOOST_TEST_CASE(&teststat_test));

    return test;
}

