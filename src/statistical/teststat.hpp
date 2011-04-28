// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_teststat_hpp
#define permory_teststat_hpp

#include "detail/config.hpp"
#include "detail/parameter.hpp"
#include "statistical/contab.hpp"

namespace Permory { namespace statistic {

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

} // namespace statistic
} // namespace Permory

#endif // include guard

