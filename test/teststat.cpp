// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#include <valarray>
#include <vector>

#include "config.hpp"
#include "parameter.hpp"
#include "statistical/algorithms.hpp"
#include "statistical/contab.hpp"
#include "statistical/teststat.hpp"
#include "statistical/testpool.hpp"

using namespace std;
using namespace Permory;

int main(int argc, char** argv) 
{
    //void Permory::Test::teststat(int a, int n) { 
    int a = atoi(argv[1]);
    int n = atoi(argv[2]);

    Con_tab<2, 3> ct;
    vector<Con_tab<2, 3> > v23(n); 
    //for (size_t i=0; i<v23.size(); ++i) v23[i].print(); cout << endl;

    ct[0][0] = 1; ct[0][1] = 2; ct[0][2] = 3;
    ct[1][0] = 4; ct[1][1] = 5; ct[1][2] = 6;
    ct.print();

    vector<valarray<int> > vv(3);
    for (size_t i=0; i<vv.size(); ++i) {
        vv[i].resize(n);
        vv[i] = i+3;
    }

    //fill_2xL_tabs(v23, vv, ct);
    //for (size_t i=0; i<v23.size(); ++i) v23[i].print(); cout << endl;
    Parameter par;
    Test_pool<2,3> pool("t", par);
    
    double d;
    vector<double> v(v23.size()); 
    if (a == 0) {
        for (int i=0; i<n; ++i) {
            d = (*pool.begin())(v23[i]);
            //d = (*(pool.begin()+1))(v23[i]);
        }
    }
    else
    {
        for_each_test_and_tab<2, 3>(v23, pool, v.begin());
        d = v.back();
    }
    PRINT(d);

}


