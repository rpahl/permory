// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_dichotom_hpp
#define permory_dichotom_hpp

#include <algorithm>
#include <numeric>
#include <set>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "config.hpp"
#include "contab.hpp"
#include "detail/recode.hpp" //dummy_code(), index_code()
#include "detail/functors.hpp" //pair_comp_2nd
#include "locus.hpp"
#include "locusdata.hpp"
#include "parameter.hpp"
#include "permutation/booster.hpp"
#include "permutation/perm.hpp"
#include "statistical/detail.hpp"
#include "statistical/testpool.hpp"

namespace Permory { namespace stat {

    // Analyze genotype data with binary/dichotomous trait
    template<int K, int L, class T = unsigned short int> class Dichotom {
        public:
            // Iterator pass through
            typedef typename std::vector<double>::const_iterator 
                const_iterator;
            const_iterator tmax_begin() const { return tMax_.begin(); }
            const_iterator tmax_end() const { return tMax_.end(); }

            Dichotom(
                    const char*,                    //char coded tests
                    const Parameter&, 
                    const std::vector<bool>& trait, //dichotomous trait
                    const Permutation* pp=0);

            // Inspection
            size_t size() const { return testPool_.size(); }

            // Test locus that is associated with the data
            template<class D> boost::shared_ptr<Locus> locus_test(
                    const Locus_data<D>&);
            // Compute permutation test statistics 
            template<class D> void permutation_test(const Locus_data<D>&);

        private:
            // This function does the "permutation work"
            template<class D> void do_permutation( 
                    const Locus_data<D>&, 
                    const std::map<D, int>&);

            Test_pool<K, L> testPool_;
            T nCases_;
            std::vector<T> trait_;
            Matrix<T> caseFreqs_;//freqs for all permutations (one row per genotype)
            boost::ptr_vector<Perm_boost<T> > boosters_;
            std::vector<double> tMax_; //max test statistics 

            // For caching purpose
            std::vector<Con_tab<K, L> > tables_; 
            std::vector<int> index_[L+1];
            Bitset2 dummy_[L+1];
    };
    // ========================================================================
    // Dichotom implementations
    template<int K, int L, class T> inline 
        Dichotom<K, L, T>::Dichotom(
                const char* tests, 
                const Parameter& par,
                const std::vector<bool>& trait,
                const Permutation* pp
                ) : trait_(trait.begin(), trait.end()), testPool_(tests, par)
        {
            if (pp) {
                size_t n = par.nperm; //number of permutations
                nCases_ = std::accumulate(trait_.begin(), trait_.end(), 0); 

                tables_.resize(n);
                tMax_.resize(n);
                caseFreqs_.resize(L+1, n); //one extra row to account for missings 

                // Create and store permutations in matrix
                boost::shared_ptr<Perm_matrix<T> > pmat(
                        new Perm_matrix<T>(n, *pp, trait_, par.useBar));

                // Prepare accelerated permutation
                boosters_.reserve(L+1); 
                for (int i=0; i<L+1; i++) {
                    boosters_.push_back(new Perm_boost<T>(pmat));
                }
            }
            else {
                //no permutations
            }
        }
    template<int K, int L, class T> template<class D> inline void 
        Dichotom<K, L, T>::do_permutation(
                const Locus_data<D>& d, const std::map<D, int>& m)
        {
            std::pair<D, int> pairs_[L+1];
            copy(m.begin(), m.end(), pairs_);
            std::set<int> nonzero_freq; //permute only alleles with frequency > 0
            for (int i=0; i<L+1; i++) {
                D a = pairs_[i].first;       //allelic code
                int cnt = pairs_[i].second; //number of occurences of this code
                if (cnt > 0) {
                    nonzero_freq.insert(i);
                    index_[i] = index_code<D>(d.begin(), d.end(), a);
                    dummy_[i] = dummy_code<D>(d.begin(), d.end(), a);
                    pairs_[i] = boosters_[i].find_most_similar(dummy_[i], cnt);
                }
            }
            size_t max_i = distance(
                    pairs_, max_element(pairs_, pairs_ + (L+1),
                        comp_second<std::pair<D, int> >()));
            caseFreqs_[max_i] = nCases_; //init with marginal sum

            // Permute for each allelic code, but *not* the one with the maximal
            // frequency, which is derived most efficiently using the marginal sum
            nonzero_freq.erase(max_i);
            BOOST_FOREACH(int i, nonzero_freq) {
                boosters_[i].permute(
                        index_[i],          //index coded data
                        dummy_[i],          //dummy coded data
                        pairs_[i].first,    //allelic code
                        &caseFreqs_[i]);    //resulting case frequenices
                caseFreqs_[max_i] -= caseFreqs_[i]; //subtraction from marginal
            }
            // add data and result of skipped one to buffer of booster 
            boosters_[max_i].add_to_buffer(dummy_[max_i]);
            boosters_[max_i].add_to_buffer(caseFreqs_[max_i]);
        }

    template<int K, int L, class T> template<class D> inline void 
        Dichotom<K, L, T>::permutation_test(
                const Locus_data<D>& d) 
        {
            assert (trait_.size() == d.size());
            std::map<D, int> m = d.unique_with_counts();
            assert (m.size() == L+1);
            do_permutation(d, m);

            // Fill contingency tables 
            typename std::map<D, int>::const_iterator mi = m.begin(); 
            int j = 0; //row index of matrix with case frequency results
            int c = 0; //column index of contingency table
            for (mi = m.begin(); mi!=m.end(); mi++) {
                bool skipThis = mi->first == d.get_undef();
                if (!skipThis) {
                    int n = mi->second;      //n = #freqs of cases + controls
                    for (int i=0; i<tables_.size(); ++i) {  
                        tables_[i][0][c] = caseFreqs_[j][i];     //cases r[j]
                        tables_[i][1][c] = n - caseFreqs_[j][i]; //controls s[j]
                    }
                    c++; //always wanted to use that in C++, huh? ;-)
                }
                j++;
            }
            for_each_test_and_tab(tables_, testPool_, tMax_.begin());
            //tables_[0].print();//XXX
        }

    template<int K, int L, class T> template<class D> inline 
        boost::shared_ptr<Locus> Dichotom<K, L, T>::locus_test(
                const Locus_data<D>& d)
        {
            std::map<D, int> m = d.unique_with_counts();
            int nMiss = 0; //number of undefined/missing values
            if (d.hasMissings()) {
                nMiss = m[d.get_undef()];
                m.erase(d.get_undef()); 
            }

            // Create contingency table
            Con_tab<K, L>* tab;
            if (d.hasMissings()) { 
                // exclude missings from analysis
                std::vector<D> v(nMiss);
                remove_copy(d.begin(), d.end(), v, d.get_undef());
                tab = new Con_tab<K, L>(trait_.begin(), trait_.end(), 
                        v.begin(), v.end());
            }
            else {
                tab = new Con_tab<K, L>(trait_.begin(), trait_.end(), 
                        d.begin(), d.end());
            }
            std::vector<double> v(testPool_.size());
            for_each_test(*tab, testPool_.begin(), testPool_.end(), v.begin());
            delete tab;
            typename boost::shared_ptr<Locus> loc = d.get_locus();
            loc->add_test_stat<K, L>(v.begin(), v.end());
            return loc;
        }

} // namespace stat
} // namespace Permory

#endif // include guard

