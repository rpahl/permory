// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#define PERMORY_TEST permutation_test
#include "detail/config.hpp"
#include "permutation/perm.hpp"
#include "permutation/booster.hpp"
#include "test.hpp"

#include "detail/parameter.hpp"
#include "detail/pair.hpp"

#include <vector>
#include <valarray>

using namespace std;
using namespace boost;
using namespace unit_test;
using namespace Permory::permutation;
using namespace Permory::detail;
using namespace Permory::gwas;


void git_test()
{
    Parameter par;
    par.useBar = false;
    par.nperm_block = 10000;

    {
        double d[] = {1, 1, 1, 1, 1};
        vector<double> phenotypes(&d[0], &d[0]+5);
        valarray<double> result(0., par.nperm_block);
        vector<int> index;
        index.push_back(1);
        index.push_back(3);

        Permutation pp;
        Perm_matrix<double> pmat(par.nperm_block, pp, phenotypes, par.useBar);
        git(pmat, index, result);

        BOOST_CHECK_CLOSE( result[0], 2., 0.0001);
    }

    {
        typedef Pair<double> P;
        P d[] = {
                make_pair(1., 2.), make_pair(1., 2.), make_pair(1., 2.),
                make_pair(1., 2.), make_pair(1., 2.)
            };
        vector<P> phenotypes(&d[0], &d[0]+5);
        valarray<P> result(make_pair(0., 0.), par.nperm_block);
        vector<int> index;
        index.push_back(1);
        index.push_back(3);

        Permutation pp;
        Perm_matrix<P> pmat(par.nperm_block, pp, phenotypes, par.useBar);
        git(pmat, index, result);

        BOOST_CHECK_CLOSE( result[0].first, 2., 0.0001);
        BOOST_CHECK_CLOSE( result[0].second, 4., 0.0001);
    }
}

test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("Functions and classes from src/permutation");

    test->add(BOOST_TEST_CASE(&git_test));

    return test;
}

