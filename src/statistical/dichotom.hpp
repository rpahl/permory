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

#include "detail/config.hpp"
#include "contab.hpp"
#include "detail/functors.hpp" //pair_comp_2nd
#include "detail/parameter.hpp"
#include "gwas/locusdata.hpp"
#include "permutation/booster.hpp"  //Bitset_with_count,
#include "permutation/perm.hpp"
#include "statistical/testpool.hpp"

namespace Permory { namespace statistic {
    using namespace permutation;

    //
    // Analyze genotype data with binary/dichotomous trait
    //
    template<uint K, uint L, class T = unsigned short int> class Dichotom {
        public:

            // Iterator pass through
            typedef typename std::vector<double>::const_iterator const_iterator;
            const_iterator tmax_begin() const { return tMax_.begin(); }
            const_iterator tmax_end() const { return tMax_.end(); }

            Dichotom(
                    const detail::Parameter&, 
                    gwas::Gwas::const_inderator ind_begin,//individuals to create
                    gwas::Gwas::const_inderator ind_end,  // dichotomous trait from
                    const Permutation* pp=0);       // pre-stored permutations

            // Inspection
            size_t size() const { return testPool_.size(); }

            // Modification
            void renew_permutations(
                    const Permutation* pp,  //creates the permutations
                    size_t nperm,           //number of permutations
                    size_t tail_size);      //parameter of the permutation booster

            // Conversion
            // Compute test statistics for the data
            template<class D> std::vector<double> test(const gwas::Locus_data<D>&);
            // Compute permutation test statistics 
            template<class D> void permutation_test(const gwas::Locus_data<D>&);

        private:
            // This function does the "permutation work"
            template<class D> void do_permutation(const gwas::Locus_data<D>&);

            std::vector<T> prepare_trait(gwas::Gwas::const_inderator begin,
                    gwas::Gwas::const_inderator end);

            std::vector<T> trait_;
            Test_pool<Con_tab<K, L> > testPool_;
            T nCases_;
            detail::Matrix<T> caseFreqs_;   //freqs for all permutations 
                                            //(one row per genotype)
            boost::ptr_vector<Perm_boost<T> > boosters_;
            std::vector<double> tMax_; //max test statistics 

            // For caching purpose
            std::vector<Con_tab<K, L> > con_tabs_; 
            std::vector<int> index_[L+1];   //L+1 integer vectors
            Bitset_with_count dummy_[L+1];  //L+1 bitsets with cached bit counts
            bool useBitarithmetic_;
    };
    // ========================================================================
    // Dichotom implementations
    template<uint K, uint L, class T> inline 
        Dichotom<K, L, T>::Dichotom(
                const detail::Parameter& par,
                gwas::Gwas::const_inderator ind_begin,
                gwas::Gwas::const_inderator ind_end,
                const Permutation* pp) 
        : trait_(prepare_trait(ind_begin, ind_end)), testPool_(par),
        useBitarithmetic_(par.useBar) 
        {
            nCases_ = std::accumulate(trait_.begin(), trait_.end(), 0); 
            bool yesPermutation = (pp != 0);
            if (yesPermutation) {
                renew_permutations(pp, par.nperm_block, par.tail_size);
            }
        }

    template<uint K, uint L, class T> inline 
        void Dichotom<K, L, T>::renew_permutations(const Permutation* pp, size_t nperm, 
                size_t tail_size)
    {
        con_tabs_.resize(nperm);
        tMax_.clear();
        tMax_.resize(nperm);
        caseFreqs_.resize(L+1, nperm); //one extra row to account for missings 

        // Create and store permutations in matrix
        boost::shared_ptr<Perm_matrix<T> > pmat(
                new Perm_matrix<T>(nperm, *pp, trait_, useBitarithmetic_));

        // Prepare permutation booster
        boosters_.clear();
        boosters_.reserve(L+1); 
        for (uint i=0; i<L+1; i++) {
            boosters_.push_back(new Perm_boost<T>(pmat, tail_size));
        }
    }

    template<uint K, uint L, class T> template<class D> inline 
        std::vector<double> Dichotom<K, L, T>::test(const gwas::Locus_data<D>& data)
        {
            // Create contingency tab and analyze it with each test of the test pool
            Con_tab<K, L>* tab;
            if (data.hasMissings()) { //requires extra work to get rid of the missings
                std::vector<D> dd; dd.reserve(data.size());
                std::vector<T> tt; tt.reserve(data.size());

                typename gwas::Locus_data<D>::const_iterator it = data.begin(); 
                BOOST_FOREACH(T t, trait_) {
                    if (*it != data.get_undef()) { //only keep the valid values
                        dd.push_back(*it);
                        tt.push_back(t);
                    }
                }
                tab = make_Con_tab<T,D,K,L>(tt.begin(), tt.end(), dd.begin(), dd.end());
            }
            else {
                tab = make_Con_tab<T,D,K,L>(
                        trait_.begin(), trait_.end(), data.begin(), data.end());
            }
            std::vector<double> v(testPool_.size());
            for_each_test(*tab, testPool_.begin(), testPool_.end(), v.begin());
            delete tab;
            return v;
        }

    template<uint K, uint L, class T> template<class D> inline void 
        Dichotom<K, L, T>::permutation_test(const gwas::Locus_data<D>& data) 
        {
            assert (trait_.size() == data.size());
            bool yesPermutation = (not boosters_.empty());
            if (yesPermutation) {
                do_permutation(data);
            }

            // Fill contingency tables 
            uint j = 0; //row index of matrix with case frequency results
            uint c = 0; //column index of contingency table

            // the unique_iterator is defined in discretedata.hpp:
            // std::map<elem_type, count_type> unique_;//unique elements with counts
            typename gwas::Locus_data<D>::unique_iterator uniques = data.unique_begin();
            for (; uniques!=data.unique_end(); uniques++) {
                bool ok = not (uniques->first == data.get_undef());
                if (ok) {
                    uint n = uniques->second; //frequency of both (cases + controls)
                    for (uint t=0; t<con_tabs_.size(); ++t) {  
                        con_tabs_[t][0][c] = caseFreqs_[j][t];     //cases r[j]
                        con_tabs_[t][1][c] = n - caseFreqs_[j][t]; //controls s[j]
                    }
                    c++; 
                }
                j++;
            }
            // For each permutation i (i.e. for each obtained contingency 
            // table) compute the max over all test statistics, say max(i), and
            // then update tMax_[i] = max(tMax_[i], max(i))
            for_each_test_and_tab(con_tabs_, testPool_, tMax_.begin());
            //con_tabs_[0].print();//XXX
        }

    template<uint K, uint L, class T> template<class D> inline void 
        Dichotom<K, L, T>::do_permutation(const gwas::Locus_data<D>& data)
        {
            if (not (data.domain_cardinality() == L+1)) {
                throw std::runtime_error("Bad domain cardinality in permutation test.");
            }
            uint boost_index[L+1]; //see below
            size_t worst_idx = 0, maxcnt = 0;
            typename gwas::Locus_data<D>::unique_iterator it = data.unique_begin();
            for (uint i=0; i<L+1; i++) {
                size_t cnt = (size_t) it->second; //#occurences of the code
                D a = it->first; //allelic code
                index_[i].clear();
                index_[i].reserve(cnt);
                index_code<D>(index_[i], data.begin(), data.end(), a);
                dummy_[i] = dummy_code<D>(data.begin(), data.end(), a);

                // Determine index of most similar dummy code in the 
                // booster's buffer. If the smallest distance between this
                // dummy code and the most similar one is smaller than 'cnt', 
                // 'cnt' will contain this distance after the call.
                boost_index[i] = boosters_[i].find_most_similar(dummy_[i], cnt);

                // Keep track of the worst boostable element, which is the one 
                // that shows both highest occurences of the code and the 
                // highest distance to "neighboured" dummy codes
                if (cnt > maxcnt) { 
                    maxcnt = cnt;
                    worst_idx = i;
                }
                it++;
            }

            // Permute for each allelic code, except the one that can be least
            // optimized/boosted for permutation. For this, the frequency is 
            // simply derived via the marginal sum.
            caseFreqs_[worst_idx] = nCases_; //init with marginal sum
            for (uint i=0; i<L+1; i++) {
                if (i != worst_idx) {
                    boosters_[i].permute(
                            index_[i],          //index coded data
                            dummy_[i],          //dummy coded data
                            boost_index[i],     //index into booster's buffer
                            &caseFreqs_[i]);    //resulting case frequenices
                    caseFreqs_[worst_idx] -= caseFreqs_[i]; 
                }
            }
            // Since it was left out, the dummy code and the resulting case 
            // frequencies of the skipped code are added post hoc "by hand" 
            // to the buffer of the booster. 
            boosters_[worst_idx].add_to_buffer(dummy_[worst_idx]);
            boosters_[worst_idx].add_to_buffer(caseFreqs_[worst_idx]);
        }

    template<uint K, uint L, class T> std::vector<T>
        Dichotom<K, L, T>::prepare_trait(
                    gwas::Gwas::const_inderator ind_begin,
                    gwas::Gwas::const_inderator ind_end)
    {
        std::vector<T> result(ind_end - ind_begin);
        transform(ind_begin, ind_end, result.begin(),
                    std::mem_fun_ref(&Individual::isAffected));
        return result;
    }

} // namespace statistic
} // namespace Permory

#endif // include guard

