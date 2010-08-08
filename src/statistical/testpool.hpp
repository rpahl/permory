// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_testpool_hpp
#define permory_testpool_hpp

#include <boost/ptr_container/ptr_vector.hpp>

#include "config.hpp"
#include "detail/parameter.hpp"
#include "teststat.hpp"

namespace Permory { namespace stat {
    template<int K, int L> class Test_pool {
        public:
            // Iterator pass through
            typedef typename boost::ptr_vector<Test_stat<K, L> >::const_iterator
                const_iterator;
            const_iterator begin() const { return ts_.begin(); }
            const_iterator end() const { return ts_.end(); }

            Test_pool(){}
            Test_pool(const char*, const Parameter&); 
            // Inspector
            size_t size() const { return ts_.size(); }

            // Modifier
            void add(const char*, const Parameter&);
            void remove(const char*);
        private:
            boost::ptr_vector<Test_stat<K, L> > ts_; //smart pointers 
    };

    // ========================================================================
    // Test_pool implementation
    template<int K, int L> inline Test_pool<K, L>::Test_pool(
            const char* tests, const Parameter& par) {
        add(tests, par);
    }
    template<> inline void Test_pool<2, 3>::add(
            const char* tests, const Parameter& par)
    {
        BOOST_FOREACH(char c, tests) {
            switch (c) {
                case 't': // standard trend test
                    ts_.push_back(new Trend);
                    break;
                case 'e': // extended trend test
                    ts_.push_back(new Trend_ext(par));
                    break;
            }
        }
    }

    template<> inline void Test_pool<2, 2>::add(
            const char* tests, const Parameter& par)
    {
        par.nperm = par.nperm; //suppress warning of unused par
        BOOST_FOREACH(char c, tests) {
            switch (c) {
                case 'c': // Qui square test
                    ts_.push_back(new Chi_squ);
                    break;
            }
        }
    }
    template<int K, int L> inline void Test_pool<K, L>::remove(
            const char* tests) {
        BOOST_FOREACH(char c, tests) {
            ts_.erase(remove_if(ts_.begin(), ts_.end(), c), ts_.end());
        }
    }

} // namespace stat
} // namespace Permory

#endif // include guard
