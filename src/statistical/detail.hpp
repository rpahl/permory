// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_statistical_detail_hpp
#define permory_statistical_detail_hpp

#include <algorithm>
#include <vector>

#include "config.hpp"
#include "contab.hpp"
#include "statistical/testpool.hpp"
#include "statistical/teststat.hpp"

namespace Permory { namespace statistic {
    template<int K, int L> inline typename std::vector<double>::iterator 
        for_each_tab(
            typename std::vector<Con_tab<K, L> >::const_iterator start,
            typename std::vector<Con_tab<K, L> >::const_iterator end,
            const typename Test_pool<K, L>::const_iterator t,
            typename std::vector<double>::iterator itResult)
    {
        while (start != end) {
            *itResult++ = (*t)(*start++);
        }
        return itResult;
    }
    template<int K, int L> inline typename std::vector<double>::iterator 
        for_each_test(
            const Con_tab<K, L>& tab,
            typename Test_pool<K, L>::const_iterator start, 
            typename Test_pool<K, L>::const_iterator end, 
            typename std::vector<double>::iterator itResult)
    {
        while (start != end) {
            *itResult++ = (*start++)(tab);
        }
        return itResult;
    }
    template<int K, int L> inline typename std::vector<double>::iterator 
        for_each_test_and_tab(
            const typename std::vector<Con_tab<K, L> >& tab,
            const Test_pool<K, L>& pool,
            typename std::vector<double>::iterator itResult)
        {
            //vector<double> v(tab.size());
            typename std::vector<double>::iterator result_begin = itResult;
            for (typename Test_pool<K, L>::const_iterator 
                    itTest = pool.begin(); itTest!=pool.end(); itTest++)   
            {
                itResult = result_begin; //for each test start from beginning
                for (typename std::vector<Con_tab<K, L> >::const_iterator 
                        itTab=tab.begin(); itTab!=tab.end(); itTab++)   
                { 
                    double d = max(*itResult, (*itTest)(*itTab));
                    *itResult++ = d;    //Tmax
                }
                // "STL-version" requires temporary vector v, and hence is slower:
                //for_each_tab<K, L>(tab.begin(), tab.end(), i, v.begin());
                //transform(v.begin(), v.end(), itResult, itResult, op_max<double>());
            }
            return itResult;
        }

} // namespace statistic
} // namespace Permory

#endif // include guard

