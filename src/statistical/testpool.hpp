// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_testpool_hpp
#define permory_testpool_hpp

#include <boost/ptr_container/ptr_vector.hpp>

#include <algorithm>
#include <set>
#include <vector>

#include "detail/config.hpp"
#include "detail/parameter.hpp"
#include "contab.hpp"
#include "teststat.hpp"

namespace Permory { namespace statistic {
    template<class T> class Test_pool {
        public:
            // Iterator pass through
            typedef typename boost::ptr_vector<Test_stat<T> >::const_iterator
                const_iterator;
            const_iterator begin() const { return ts_.begin(); }
            const_iterator end() const { return ts_.end(); }

            Test_pool(){}
            Test_pool(const detail::Parameter&);
            // Inspector
            size_t size() const { return ts_.size(); }

            // Modifier
            void add(const detail::Parameter&);
            void add(Test_stat<T>*);
            void remove(const std::set<detail::Test_type>&);
        private:
            boost::ptr_vector<Test_stat<T> > ts_;
    };

    // Test_pool implementation
    // ========================================================================
    template<class T> inline Test_pool<T>::Test_pool(const detail::Parameter& par)
    {
        add(par);
    }
    template<> inline void Test_pool<Con_tab<2, 3> >::add(const detail::Parameter& par)
    {
        using namespace detail;
        BOOST_FOREACH(Test_type t, par.tests) {
            switch (t) {
                case trend: // standard trend test
                    ts_.push_back(new Trend);
                    break;
                case trend_extended: // extended trend test
                    ts_.push_back(new Trend_ext(par));
                    break;
                default:
                    //ignore unsupported tests
                    break;
            }
        }
    }

    template<> inline void Test_pool<Con_tab<2, 2> >::add(const detail::Parameter& par)
    {
        using namespace detail;
        BOOST_FOREACH(Test_type t, par.tests) {
            switch (t) {
                case chisq: // Qui square test
                    ts_.push_back(new Chi_squ);
                    break;
                default:
                    //ignore unsupported tests
                    break;
            }
        }
    }

    template<class T> inline void Test_pool<T>::add(Test_stat<T>* ts)
    {
        ts_.push_back(ts);
    }

    template<class T> inline void Test_pool<T>::remove(
            const std::set<detail::Test_type>& tests)
    {
        BOOST_FOREACH(detail::Test_type t, tests) {
            ts_.erase(remove_if(ts_.begin(), ts_.end(), t), ts_.end());
        }
    }

    // ========================================================================
    // Non-member functions
    template<class T> inline typename std::vector<double>::iterator
        for_each_tab(
            typename std::vector<T>::const_iterator start,
            typename std::vector<T>::const_iterator end,
            const typename Test_pool<T>::const_iterator t,
            typename std::vector<double>::iterator itResult)
    {
        while (start != end) {
            *itResult++ = (*t)(*start++);
        }
        return itResult;
    }
    template<class T> inline typename std::vector<double>::iterator
        for_each_test(
            const T& tab,
            typename Test_pool<T>::const_iterator start,
            typename Test_pool<T>::const_iterator end,
            typename std::vector<double>::iterator itResult)
    {
        while (start != end) {
            *itResult++ = (*start++)(tab);
        }
        return itResult;
    }
    template<class T> inline typename std::vector<double>::iterator
        for_each_test_and_tab(
            const typename std::vector<T>& tab,
            const Test_pool<T>& pool,
            typename std::vector<double>::iterator itResult)
        {
            //vector<double> v(tab.size());
            typename std::vector<double>::iterator result_begin = itResult;
            for (typename Test_pool<T>::const_iterator
                    itTest = pool.begin(); itTest!=pool.end(); itTest++)
            {
                itResult = result_begin; //for each test start from beginning
                for (typename std::vector<T>::const_iterator
                        itTab=tab.begin(); itTab!=tab.end(); itTab++)
                {
                    double d = std::max(*itResult, (*itTest)(*itTab));
                    *itResult++ = d;    //Tmax
                }
                // "STL-version" requires temporary vector v, and hence is slower:
                //for_each_tab<T>(tab.begin(), tab.end(), i, v.begin());
                //transform(v.begin(), v.end(), itResult, itResult, op_max<double>());
            }
            return itResult;
        }
} // namespace statistic
} // namespace Permory

#endif // include guard
