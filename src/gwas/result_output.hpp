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
        for (itLocus; itLocus!=study.end(); ++itLocus) {
            bool isPoly = itLocus->isPolymorph();
            double maf = itLocus->maf();
            nbad += not isPoly;
            nmin += maf < par->min_maf && isPoly;
            nmax += maf > par->max_maf && isPoly;
        }
        size_t sum = nbad + nmin + nmax;

        myout << verbose << stdpre << "Loci filtered out:" << endl;
        myout << indent(4) << "non-polymorphic: " << nbad << endl; 
        if (par->min_maf > 0) {
            myout << indent(4) << "maf < " << par->min_maf << ": " << nmin << endl; 
        }
        if (par->max_maf < 1) {
            myout << indent(4) << "maf > " << par->max_maf << ": " << nmax << endl; 
        }
        myout << indent(4) << "Total: " << sum << endl;
        myout << normal << stdpre << m - sum << " of " << m <<
            " loci were used for analysis." << endl;
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
            std::cerr << "Logfile: overwrite `" << fn << "'? (y/n) ";
            std::cin >> c;
            if (c != 'y') {
                return;
            }
        }
        File_out out(fn);

        Gwas::const_iterator itLoc;
        bool hasChr = false;
        size_t wrs = 0,
               wgene = 0;
        for (itLoc = study.begin(); itLoc != study.end(); itLoc++) {
            hasChr = hasChr || itLoc->chr() != Locus::none;
            wrs = std::max(wrs, itLoc->rs().size());
            wgene = std::max(wgene, itLoc->gene().size());
        }

        // Header 
        std::string s = "Serial no.";
        out << left << setw(12) << s;
        s = "rsID";
        if (wrs > 0) {
            out << left << setw(wrs+3) << s;
        }
        s = "chr";
        if (hasChr) {
            out << left << setw(7) << s;
        }
        s = "gene";
        if (wgene > 0) {
            out << left << setw(wgene+3) << s;
        }
        size_t ntest = par->tests.size();
        for (std::set<Test_type>::const_iterator it=par->tests.begin();
                it!=par->tests.end(); ++it)
        {
            s = test_type_to_string(*it);
            out << left << setw(16) << s;
        }
        if (ntest > 1) {
            out << left << setw(12) << "T_max";
        }
        out << left << setw(15) << "p.unadjusted" << 
            left << setw(12) << "p.adjusted";
        if (par->pval_counts) {
            out << left << setw(10) << "p.counts";
        }
        out << endl;

        //std::deque<double>::const_iterator itP = pvals_unadj.begin();
        std::deque<double>::const_iterator itPadj = pp.begin();
        std::deque<size_t>::const_iterator itPcnt = pval_cnts.begin();
        for (itLoc = study.begin(); itLoc != study.begin()+pp.size(); itLoc++) {
            out << left << setw(12) << itLoc->id();
            if (wrs > 0) {
                out << left << setw(wrs+3) << itLoc->rs();
            }
            s = "chr";
            if (hasChr) {
                out << left << setw(7) << itLoc->rs();
            }
            s = "gene";
            if (wgene > 0) {
                out << left << setw(wgene+3) << itLoc->gene();
            }
            for (std::vector<double>::const_iterator it=itLoc->test_stats().begin();
                    it!=itLoc->test_stats().end(); ++it) {
                out << left << setw(16) << setprecision(8) << fixed << *it;
            }
            if (ntest > 1) {
                out << left << setw(12) << setprecision(6) << fixed << itLoc->tmax();
            }
            out << left << setw(15) << setprecision(6) << scientific << 0.0;//*itP++;TODO
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


