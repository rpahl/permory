// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_testpool_hpp
#define permory_testpool_hpp

#include <boost/ptr_container/ptr_vector.hpp>

#include <set>

#include "config.hpp"
#include "detail/parameter.hpp"
#include "teststat.hpp"

namespace Permory { namespace statistic {
    template<uint K, uint L> class Test_pool {
        public:
            // Iterator pass through
            typedef typename boost::ptr_vector<Test_stat<K, L> >::const_iterator
                const_iterator;
            const_iterator begin() const { return ts_.begin(); }
            const_iterator end() const { return ts_.end(); }

            Test_pool(){}
            Test_pool(const Parameter&); 
            // Inspector
            size_t size() const { return ts_.size(); }

            // Modifier
            void add(const Parameter&);
            void remove(const std::set<Test_type>&);
        private:
            boost::ptr_vector<Test_stat<K, L> > ts_; 
    };

    // ========================================================================
    // Test_pool implementation
    template<uint K, uint L> inline Test_pool<K, L>::Test_pool(const Parameter& par) 
    {
        add(par);
    }
    template<> inline void Test_pool<2, 3>::add(const Parameter& par)
    {
        BOOST_FOREACH(detail::Test_type t, par.tests) {
            switch (t) {
                case trend: // standard trend test
                    ts_.push_back(new Trend);
                    break;
                case trend_extended: // extended trend test
                    ts_.push_back(new Trend_ext(par));
                    break;
            }
        }
    }

    template<> inline void Test_pool<2, 2>::add(const Parameter& par)
    {
        BOOST_FOREACH(detail::Test_type t, par.tests) {
            switch (t) {
                case chisq: // Qui square test
                    ts_.push_back(new Chi_squ);
                    break;
            }
        }
    }
    template<uint K, uint L> inline void Test_pool<K, L>::remove(
            const std::set<Test_type>& tests)
    {
        BOOST_FOREACH(detail::Test_type t, tests) {
            ts_.erase(remove_if(ts_.begin(), ts_.end(), t), ts_.end());
        }
    }

} // namespace statistic
} // namespace Permory

#endif // include guard
