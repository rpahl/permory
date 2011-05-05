// Copyright (c) 2010 Roman Pahl
//               2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_teststat_hpp
#define permory_teststat_hpp

#include <vector>

#include "detail/config.hpp"
#include "detail/parameter.hpp"
#include "detail/pair.hpp"
#include "statistical/contab.hpp"
#include "gwas/locusdata.hpp"

namespace Permory { namespace statistic {

    using Permory::detail::Pair;
    using Permory::detail::make_pair;
    using Permory::gwas::Locus_data;

    template<class T> class Test_stat {
        public: 
            virtual ~Test_stat(){}
            double operator()(const T& tab) const {
                return do_operator(tab);
            }
        private:
            virtual double do_operator(const T&) const = 0;
    };

    //! @brief Standard trend test with pooled variance estimator
    class Trend : public Test_stat<Con_tab<2,3> > {
        private:
            double do_operator(const Con_tab<2,3>& ct) const; 
    };

    //! @brief Trend test extended for different weights and variance estimators
    class Trend_ext : public Test_stat<Con_tab<2,3> > {
        public:
            Trend_ext(const detail::Parameter& par) : ve(par.ve) { 
                w[0] = par.w[0]; 
                w[1] = par.w[1]; 
                w[2] = par.w[2]; 
            }
        private:
            double do_operator(const Con_tab<2,3>& ct) const;
            detail::Var_estimate ve;    //variance estimator
            double w[3];                //weights
    };

    class Chi_squ : public Test_stat<Con_tab<2,2> > {
        private:
            double do_operator(const Con_tab<2,2>& ct) const;
    };


    template<class T>
    class Trend_continuous : public Test_stat<std::vector<Pair<T> > > {
        public:
            // Ctors
            Trend_continuous(const std::vector<T>& trait)
                    : mu_y_(calculate_mu_y(trait)),
                      nomdenom_buf_(create_buffer(trait)),
                      sum_(accumulate(nomdenom_buf_.begin(),
                           nomdenom_buf_.end(), make_pair<T>(0,0)))
                { }

            // Inspectors
            T get_mu_y() const { return mu_y_; }
            const std::vector<Pair<T> >& get_buffer() const
                { return nomdenom_buf_; }
            T get_mu_j() const { return mu_j_; }
            T get_invariant() const { return invariant_; }
            T get_denom_invariant() const { return invariant_.second; }
            Pair<T> get_sum() const { return sum_; }

            // Modifiers
            template<class D> void update(const Locus_data<D>& data);

        protected:
            template<class D> Pair<T> calculate_invariant(
                    const Locus_data<D>& data) const;
            T calculate_mu_y(const std::vector<T>& trait) const;
            T calculate_mu_j(const Locus_data<uint>& data) const;
            template<class D> T calculate_mu_j(const Locus_data<D>& data) const;

        private:
            double do_operator(const std::vector<Pair<T> >& tab) const;

            std::vector<Pair<T> > create_buffer(const std::vector<T>& trait)
                    const;

            T mu_y_;            // mean of phenotypes:
                                //   1/N * \sum_{i=1}^{N} Y_i
            T mu_j_;            // mean of genotypes of marker j:
                                //   1/N * \sum_{i=1}^{M} X_ji
            std::vector<Pair<T> > nomdenom_buf_;
                                // nominator-denominator buffer:
                                //   first  = \sum_{i=1}^{N} (Y_i - \mu_y)
                                //   second = \sum_{i=1}^{N} (Y_i - \mu_y)^2
            const Pair<T> sum_; // Sum of all elements in nomdenom_buf_.
            Pair<T> invariant_; // invariant summand of nominator and
                                // denominator in marker j:
                                // nominator:
                                //   \mu_j * \sum_{i=1}^{N} (Y_i - \mu_y)
                                // denominator:
                                //   \mu_j^2 * \sum_{i=1}^{N} (Y_i - \mu_y)^2
    };

    // ========================================================================
    // Test_stat implementations
    inline double Trend::do_operator(const Con_tab<2,3>& ct) const
    {
        // Extract data from contingency table
        uint r[3], n[3]; //r[0] and n[0] are ignored since weight is 0 anyway
        for (uint i=1; i<3; ++i) {
            r[i] = ct(0, i);
            n[i] = ct.colsum(i);
        }
        double Swr = double(r[1] + 2*r[2]),    // weighted sum of the r[i]
               Swn = double(n[1] + 2*n[2]),    // weighted sums of the n[i]
               Swwn = double(n[1] + 4*n[2]);   
        double R = double(ct.rowsum(0)),
               N = R + double(ct.rowsum(1));

        // denominator, i.e. (pooled) variance estimator 
        double den = R*(N-R)*(N*Swwn - Swn*Swn);

        if (den > 0) {
            double nom = N*Swr - R*Swn;    //nominator, i.e. effect size parameter
            nom *= nom; 
            nom *= N;
            return nom/den; 
        }
        else {
            return 0;
        }
    }

    inline double Trend_ext::do_operator(const Con_tab<2,3>& ct) const
    { 
        using namespace detail;

        // Extract data from contingency table
        double r[3], s[3];
        for (uint i=0; i<3; ++i) {
            r[i] = double(ct(0, i));
            s[i] = double(ct(1, i));
        }
        double R = double(ct.rowsum(0)),
               S = double(ct.rowsum(1)),
               N = R + S;

        double Swr = 0, Swwr = 0, // weighted sums of the r[i]
               Sws = 0, Swws = 0, // weighted sums of the s[i]
               Swn = 0, Swwn = 0; // weighted sums of (r[i] + s[i]) 

        double d; //aux
        for(uint i=0; i<3; ++i) { 
            d = w[i]*r[i]; 
            Swr += d;
            Swwr += w[i]*d;
            d = w[i]*s[i]; 
            Sws += d;
            Swws += w[i]*d;
        }

        double nom = (S*Swr - R*Sws)/N; //nominator, i.e. effect size parameter
        nom *= nom; 

        double den; // denominator, i.e. variance estimator 
        switch (ve) //which kind of variance estimator to use?
        {
            case separately: // Cases and controls separately
                den = (S*S*(Swwr - Swr*Swr/R) + R*R*(Swws - Sws*Sws/S)) / (N*N);
                break;

            case controls: // Controls only 
                den = R*(Swws - Sws*Sws/S)/N;
                break;

            default: // Pooled variance
                Swn = Swr + Sws; Swwn = Swwr + Swws;
                den = R*S*(Swwn - Swn*Swn/N) / (N*N);
                break;
        }

        if (den > 0) 
            return nom/den; 
        else {
            return 0;
        }
    }

    inline double Chi_squ::do_operator(const Con_tab<2,2>& ct) const
    {
        // Extract data from contingency table
        double a = double(ct(0, 0)),
               b = double(ct(0, 1)),
               c = double(ct(1, 0)),
               d = double(ct(1, 1));
        double N = a + b + c + d;
        double den = (a+b)*(c+d)*(a+c)*(b+d);

        if (den > 0) {
            double nom = a*d - b*c;
            nom *= nom;
            nom *= N;
            return nom/den;
        }
        else {
            return 0;
        }
    }


    template<class T> template<class D> inline
    Pair<T> Trend_continuous<T>::calculate_invariant(
                const Locus_data<D>& data) const
    {
        if (data.hasMissings()) {
            Pair<T> result(0, 0);
            typename Locus_data<D>::const_iterator it = data.begin();
            BOOST_FOREACH(Pair<T> x, nomdenom_buf_) {
                if (not (*it++ == data.get_undef())) {
                    result += x;
                }
            }
            result.first *= mu_j_;
            result.second *= mu_j_ * mu_j_;
            return result;
        }
        else {
            return make_pair<T>(0, sum_.second * mu_j_ * mu_j_);
        }
    }

    template<class T> inline
    T Trend_continuous<T>::calculate_mu_y(const std::vector<T>& trait) const
    {
        return std::accumulate(trait.begin(), trait.end(), T(0)) / trait.size();
    }

    template<class T> inline std::vector<Pair<T> >
    Trend_continuous<T>::create_buffer(const std::vector<T>& trait) const
    {
        std::vector<Pair<T> > result(trait.size());
        for (size_t i = 0; i < trait.size(); ++i) {
            const T tmp = trait[i] - mu_y_;
            result[i] = make_pair(tmp, tmp*tmp);
        }
        return result;
    }

    template<class T> inline
    T Trend_continuous<T>::calculate_mu_j(const Locus_data<uint>& data) const
    {
        if (data.hasMissings()) {
            std::vector<uint> v;
            data.get_data_without_missings(&v);
            return std::accumulate(v.begin(), v.end(), 0.) / v.size();
        }
        else {
            return std::accumulate(data.begin(), data.end(), 0.) / data.size();
        }
    }
    template<class T> template<class D> inline
    T Trend_continuous<T>::calculate_mu_j(const Locus_data<D>& data) const
    {
        gwas::Locus_data<uint> *numeric = data.as_numeric();
        T temp = calculate_mu_j(*numeric);
        delete numeric;
        return temp;
    }

    template<class T> template<class D> inline
    void Trend_continuous<T>::update(const Locus_data<D>& data)
    {
        mu_j_ = calculate_mu_j(data);
        invariant_ = calculate_invariant(data);
    }

    template<class T> inline
    double Trend_continuous<T>::do_operator(const std::vector<Pair<T> >& tab) const
    {
        double nominator = -invariant_.first;
        double denominator = invariant_.second;
        // Skip index 0 as product in nominator and denominator will always
        // be 0.
        for (size_t i = 1; i < tab.size(); ++i) {
            detail::Pair<T> x = tab[i];
            nominator += i * x.first;
            denominator += (i * i - 2 * i * mu_j_) * x.second;
        }
        nominator *= nominator;
        return (nominator / denominator);
    }


} // namespace statistic
} // namespace Permory

#endif // include guard

