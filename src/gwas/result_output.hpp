// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_write_result_hpp
#define permory_write_result_hpp

#include "gwas.hpp"
#include "io/file_out.hpp"
#include "io/output.hpp"
#include "statistical/pvalue.hpp"

namespace Permory { namespace gwas {

    //
    // Output results to console
    void result_to_console(detail::Parameter* par, io::Myout& myout, 
            const Gwas& study)
    {
        using namespace std;
        using namespace boost;
        using namespace io;
        size_t m = study.m();
        size_t nbad = 0;
        size_t nmin = 0;
        size_t nmax = 0;
        Gwas::const_iterator itLocus = study.begin(); 
        for (; itLocus!=study.end(); ++itLocus) {
            bool isPoly = itLocus->isPolymorph();
            double maf = itLocus->maf();
            nbad += not isPoly;
            nmin += maf < par->min_maf && isPoly;
            nmax += maf > par->max_maf && isPoly;
        }
        size_t sum = nbad + nmin + nmax;

        myout << verbose << stdpre << "Markers excluded from analyis:" << endl;
        myout << indent(4) << nbad << " non-polymorphic" << endl; 
        if (par->min_maf > 0) {
            myout << indent(4) << nmin << " with maf < " << par->min_maf << endl; 
        }
        if (par->max_maf < 0.5) {
            myout << indent(4) << nmax << " with maf > " << par->max_maf << endl; 
        }
        sum > m ? sum = m : sum;
        myout << indent(4) << sum << " in total" << endl;
        myout << normal << stdpre << m - sum << " (of " << m <<
            ") markers were used for analysis." << endl;
    }

    void result_to_file(
            detail::Parameter* par, 
            const Gwas& study, 
            const std::deque<size_t>& pval_cnts,
            const std::string& fn)
    {
        using namespace std;
        using namespace detail;
        using namespace io;
        using namespace statistic;
        size_t m = study.m();
        if (m < pval_cnts.size()) {
            throw std::length_error("More p-values than markers.");
        }
        std::deque<double> pp = single_step_pvalues(pval_cnts, par->nperm_total);

        // Handle existing file
        File_handle file(fn);
        if (bfs::exists(*file) && (par->interactive)) {
            char c;
            std::cout << "Results: overwrite `" << fn << "'? (y/n) ";
            std::cin >> c;
            if (c != 'y') {
                return;
            }
        }
        File_out out(fn);

        // Determine maximal widhts to be set for the different columns
        bool hasChr = false;
        size_t wrs = 0;     //max width for rsIDs
        size_t wgene = 0;   //max width of gene names
        Gwas::const_iterator itLoc;
        for (itLoc = study.begin(); itLoc != study.end(); itLoc++) {
            hasChr = hasChr || itLoc->chr() != Locus::none;
            wrs = std::max(wrs, itLoc->rs().size());
            wgene = std::max(wgene, itLoc->gene().size());
        }

        // Write header 
        std::string s = "Serial_no";
        size_t w = max(int(ceil(log10(m))), 9)+4;
        out << left << setw(w) << s;
        s = "rsID";
        if (wrs > 0) {
            out << left << setw(wrs+4) << s;
        }
        s = "chr";
        if (hasChr) {
            out << left << setw(7) << s;
        }
        s = "gene";
        if (wgene > 0) {
            out << left << setw(wgene+4) << s;
        }
        s = "maf";
        out << left << setw(9) << s;

        size_t ntest = par->tests.size();
        for (std::set<Test_type>::const_iterator it=par->tests.begin();
                it!=par->tests.end(); ++it)
        {
            s = test_type_to_string(*it);
            out << left << setw(13) << s;
        }
        if (ntest > 1) {
            out << left << setw(10) << "T_max";
        }
        out << left << setw(14) << "p_raw" << 
            left << setw(12) << "p_adjusted";
        if (par->pval_counts) {
            out << left << setw(10) << "p.counts";
        }
        out << endl;

        // Write results
        std::deque<double>::const_iterator itPadj = pp.begin();
        std::deque<size_t>::const_iterator itPcnt = pval_cnts.begin();
        for (itLoc = study.begin(); itLoc != study.begin()+pp.size(); itLoc++) {
            if (not itLoc->hasTeststat()) {
                continue;
            }
            size_t w = max(int(ceil(log10(m))), 9)+4;
            out << left << setw(w) << itLoc->id();
            if (wrs > 0) {
                out << left << setw(wrs+4) << itLoc->rs();
            }
            s = "chr";
            if (hasChr) {
                out << left << setw(7) << itLoc->chr();
            }
            s = "gene";
            if (wgene > 0) {
                out << left << setw(wgene+4) << itLoc->gene();
            }
            out << left << setw(9) << setprecision(3) << itLoc->maf();
            for (std::vector<double>::const_iterator it=itLoc->test_stats().begin();
                    it!=itLoc->test_stats().end(); ++it) {
                out << left << setw(13) << setprecision(5) << fixed << *it;
            }
            if (ntest > 1) {
                out << left << setw(10) << setprecision(5) << fixed << itLoc->tmax();
            }
            out << left << setw(14) << setprecision(4) << scientific << 
                1.0 - gsl_cdf_chisq_P(itLoc->tmax(), 1); 
            out << left << setw(12) << setprecision(ceil(log10(par->nperm_total))) <<
                fixed << *itPadj++;
            if (par->pval_counts) {
                out << left << setw(10) << *itPcnt++;
            }
            out << endl;
        }
    }

} // namespace gwas
} // namespace Permory

#endif // include guard


