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
        using namespace Permory::detail;
        size_t m = study.m();
        size_t nbad = 0;
        size_t nmin = 0;
        size_t nmax = 0;
        Gwas::const_iterator itLocus = study.begin(); 
        for (; itLocus!=study.end(); ++itLocus) {
            bool isPoly = itLocus->isPolymorph();
            double maf = itLocus->maf("pooled");
            nbad += not isPoly;
            nmin += maf < par->min_maf && isPoly;
            nmax += maf > par->max_maf && isPoly;
        }
        size_t sum = nbad + nmin + nmax;

        myout << verbose << stdpre << "Markers excluded from analysis:" << endl;
        myout << indent(4) << nbad << " non-polymorphic" << endl; 
        if (par->min_maf > 0) {
            myout << indent(4) << nmin << " with maf < " << par->min_maf << endl; 
        }
        if (par->max_maf < 0.5) {
            myout << indent(4) << nmax << " with maf > " << par->max_maf << endl; 
        }
        sum > m ? sum = m : sum;
        myout << indent(4) << sum << " in total" << endl;
        size_t num_analyzed_markers = m - sum;
        myout << normal << stdpre << num_analyzed_markers << " (of " << m <<
            ") markers were used for analysis." << endl;
        myout << normal << stdpre << "Effective number of independent tests (markers) = " <<
            study.meff(num_analyzed_markers) << endl;
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

        // Determine maximal widths to be set for the different columns
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
        size_t w = max(int(ceil(log10(m))), 9)+4;
        out << right << setw(w) << "Serial_no";
        if (wrs > 0) {
            out << right << setw(wrs+4) << "rsID";
        }
        if (hasChr) {
            out << right << setw(7) << "chr";
        }
        if (wgene > 0) {
            out << right << setw(wgene+4) << "gene";
        }
        out << right << setw(13) << "MAF_pooled";

        size_t ntest = par->tests.size();
        Enum_converter ec;
        for (std::set<Test_type>::const_iterator it=par->tests.begin();
                it!=par->tests.end(); ++it)
        {
            out << right << setw(13) << ec.key_to_string<Test_type>(*it);
        }
        if (ntest > 1) {
            out << right << setw(10) << "T_max";
        }
        out << right << setw(14) << "P_raw" << 
            right << setw(12) << "P_adjusted";
        if (par->pval_counts) {
            out << right << setw(10) << "P.counts";
        }
        out << endl;

        // Write results
        std::deque<double>::const_iterator itPadj = pp.begin();
        std::deque<size_t>::const_iterator itPcnt = pval_cnts.begin();
        for (itLoc = study.begin(); itLoc != study.begin()+pp.size(); itLoc++) {
            size_t w = max(int(ceil(log10(m))), 9)+4;
            out << right << setw(w) << itLoc->id();
            if (wrs > 0) {
                out << right << setw(wrs+4) << itLoc->rs();
            }
            if (hasChr) {
                out << right << setw(7) << itLoc->chr();
            }
            if (wgene > 0) {
                out << right << setw(wgene+4) << itLoc->gene();
            }
            out << right << setw(13) << setprecision(3) << itLoc->maf("pooled");


            if (itLoc->hasTeststat()) {
                for (std::vector<double>::const_iterator it=itLoc->test_stats().begin();
                        it!=itLoc->test_stats().end(); ++it) {
                    out << right << setw(13) << setprecision(5) << fixed << *it;
                }
                if (ntest > 1) {
                    out << right << setw(10) << setprecision(5) << fixed << itLoc->tmax();
                }
                out << right << setw(14) << setprecision(4) << scientific << 
                    1.0 - gsl_cdf_chisq_P(itLoc->tmax(), 1); 
                out << right << setw(12) << setprecision(ceil(log10(par->nperm_total))) <<
                    fixed << *itPadj;
                if (par->pval_counts) {
                    out << right << setw(10) << *itPcnt;
                }
            }
            else {
                for (size_t i=0; i<ntest; ++i) {
                    out << right << setw(13) << "NA";
                }
                if (ntest > 1) {
                    out << right << "NA";
                }
                out << right << setw(14) << "NA";
                out << right << setw(12) << "NA";
                if (par->pval_counts) {
                    out << right << setw(10) << "NA";
                }
            }
            out << endl;
            itPadj++;
            itPcnt++;
        }
    }

} // namespace gwas
} // namespace Permory

#endif // include guard


