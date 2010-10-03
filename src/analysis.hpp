// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_analysis_hpp
#define permory_analysis_hpp

#include <set>
#include <string>
#include <vector>

#include <boost/progress.hpp>

#include "detail/config.hpp"
#include "detail/parameter.hpp"
#include "gwas.hpp"
#include "locusdata.hpp"
#include "io/output.hpp"
#include "permutation/perm.hpp"
#include "read_phenotype_data.hpp"
#include "read_locus_data.hpp"
#include "statistical/dichotom.hpp"
#include "statistical/pvalue.hpp"
#include "write_result.hpp"

namespace Permory 
{
    std::vector<Individual> create_study_sample(
            detail::Parameter* par, io::Myout& myout)
    {
        using namespace std;
        using namespace boost;
        using namespace detail;
        using namespace io;
        vector<Individual> v;
        string fn = par->fn_trait;
        if (not fn.empty()) {
            myout << normal << stdpre << "Read trait from file..." << endl;
            Datafile_format dff = detect_phenotype_data_format(fn);
            par->phenotype_data_format = dff;
            myout << verbose << indent(4) << "`" << par->fn_trait << "' -> assumed " <<
                "file format: " << datafile_format_to_string(dff) << endl;
            read_individuals(*par, fn, &v);
        }
        else {
            // No file so create individuals by specified numbers
            size_t n = par->ncase + par->ncontrol;
            if (n == 0) {
                throw runtime_error("Sample size must not be 0");
            }
            myout << normal << stdpre << "Create trait via given numbers..." << endl;
            v = case_control_sample(n, double(par->ncase)/double(n)); //individual.hpp
        }
        return v;
    }

    void scan_loci(Gwas* the_study, detail::Parameter* par, io::Myout& myout)
    {
        using namespace std;
        using namespace boost;
        using namespace io;
        myout << normal << stdpre << "Scanning marker data..." << endl;
        BOOST_FOREACH(string fn, par->fn_marker_data) {
            Datafile_format dff = detect_marker_data_format(fn, par->undef_allele_code);
            myout << verbose << indent(4) << "`" << fn << "' -> assumed file format: " << 
                detail::datafile_format_to_string(dff);
            if (dff == unknown) {
                par->fn_bad_files.insert(fn);
                myout << " - will be ignored." << endl;
                continue;
            }
            myout << endl;
            read_loci(dff, fn, the_study->pointer_to_loci());
        }
        if (not par->fn_bad_files.empty()) {
            BOOST_FOREACH(std::string s, par->fn_bad_files) {
                par->fn_marker_data.erase(s);
            }
            cerr << warnpre << "Warning: " << par->fn_bad_files.size() << " data file(s) " <<
                "will NOT be processed due to unrecognized data format." << endl;
        }
    }

    void present_results(detail::Parameter* par, io::Myout& myout,
            Gwas& study, std::deque<double> tmax)
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

        deque<double> t_orig(study.m());
        transform(study.begin(), study.end(), t_orig.begin(), 
                mem_fun_ref(&Locus::tmax));
        sort(t_orig.begin(), t_orig.end(), greater<double>());

        detail::print_seq(t_orig.begin(), t_orig.begin()+50);
        //cout << endl;
        cout << endl;
        //detail::print_seq(tmax.begin(), tmax.begin()+50);
        //deque<size_t> counts = statistic::single_step_counts(t_orig, tmax);
        //detail::print_seq(counts.begin(), counts.end());
        cout << endl;
        deque<double> pvals = statistic::single_step_pvalues(t_orig, tmax);
        size_t offset = std::min(pvals.size(), size_t(50));
        detail::print_seq(pvals.begin(), pvals.begin()+offset);
        cout << endl;
    }

    //
    //  Compute test statistics and perform permutation test
    //
    template<uint K, uint L> void analyze_dichotom(
            detail::Parameter* par, io::Myout& myout, 
            const std::vector<bool>& trait, Gwas* study,
            std::set<char> the_domain=std::set<char>())
    {
        using namespace std;
        using namespace boost;
        using namespace io;

        // Prepare progress bar
        size_t m = study->m();
        size_t perm_todo = par->nperm_total;        //remaining permutations
        scoped_ptr<progress_display> pprogress;
        bool show_progress = not par->quiet;
        if (show_progress) {
            double d = ceil(double(perm_todo)/double(par->nperm_block));
            pprogress.reset(new progress_display(m*size_t(d), std::cout, "","",""));
        }

        deque<double> tmax;                         //permutation Tmax statistics
        permutation::Permutation pp(par->seed);     //permutation factory
        Gwas::iterator itLocus = study->begin();    //points to current locus
        statistic::Dichotom<K,L> stat(*par, trait); //handles the statistic stuff
        bool isFirstRound = true;
        while (perm_todo > 0) {
            size_t nperm = par->nperm_block;
            if (nperm > perm_todo) {
                nperm = perm_todo;
            }
            stat.renew_permutations(&pp, nperm, par->tail_size); //fresh random numbers
            itLocus = study->begin(); 

            BOOST_FOREACH(string fn, par->fn_marker_data) { 
                Datafile_format dff = detect_marker_data_format(fn, par->undef_allele_code);
                Locus_data_reader<char> loc_reader(fn, par->undef_allele_code);

                while (loc_reader.hasData()) {
                    Locus_data<char> data(loc_reader.get_next(), par->undef_allele_code);
                    if (data.size() == 0) {
                        continue;
                    }

                    if (par->gen_type == genotype) {
                        // Handle haplotype-based data formats 
                        if (dff == presto || dff == plink_tped) { 
                            data = data.merge_alleles_to_genotypes(2);
                        }
                        //ensure that all possible genotypes appear in the domain
                        data.add_to_domain(the_domain); 
                    }

                    if (data.size() != study->sample_size()) {
                        // Provide some more information in case of this error
                        // as it may be hard to spot in large data sets
                        std::string s = "Sample size (";
                        s.append(boost::lexical_cast<string>(study->sample_size()));
                        s.append(") does not match with length of marker data (");
                        s.append(boost::lexical_cast<string>(data.size()));
                        s.append(").\n");
                        throw std::runtime_error(s);
                    }

                    if (isFirstRound) { 
                        // The non-permutation stuff needs only to be done once 
                        bool isPoly = data.isPolymorph();
                        itLocus->set_polymorph(isPoly); 
                        double maf = data.maf(genotype);
                        itLocus->set_maf(maf); 
                        bool mafOk = maf > par->min_maf && maf < par->max_maf;
                        if (isPoly && mafOk) {
                            itLocus->add_test_stats(stat.test(data)); 
                        }
                    }
                    if (itLocus->hasTeststat() && perm_todo > 0) {
                        stat.permutation_test(data);   //permute and compute
                    }
                    if (show_progress) {
                        ++(*pprogress);
                    }
                    itLocus++;
                }
            }
            perm_todo -= nperm;
            copy(stat.tmax_begin(), stat.tmax_end(), back_inserter(tmax));
            isFirstRound = false;
        }
        myout << endl;
        present_results(par, myout, *study, tmax);
    }

    void gwas_analysis(detail::Parameter* par, io::Myout& myout)
    {
        using namespace std;
        using namespace io;

        // Phenotype data (either read in or create)
        Gwas study(create_study_sample(par, myout)); 
        if (not par->fn_trait.empty()) {
            myout << normal << indent(4) << "Found ";
        }
        else {
            myout << normal << indent(4) << "Created ";
        }
        myout << study.ncase() << " cases and " << study.ncontrol() << 
            " controls." << endl;

        // Create *dichotomous* trait/phenotype data
        vector<bool> trait(study.sample_size());
        transform(study.ind_begin(), study.ind_end(), trait.begin(),
                mem_fun_ref(&Individual::isAffected));

        // Loci information read in
        scan_loci(&study, par, myout);
        size_t m = study.m();
        if (study.m() == 0) {
            throw runtime_error("No marker loci found.");
        }
        myout << normal << stdpre << "Found " << m << " loci." << endl << endl; 

        myout << normal << stdpre << "Starting analysis..." << endl;
        std::set<char> data_domain;
        switch (par->gen_type) {
            case genotype:  //2x3 contingency table analysis
                data_domain.insert('0');
                data_domain.insert('1');
                data_domain.insert('2');
                analyze_dichotom<2,3>(par, myout, trait, &study, data_domain);
                break;
            case haplotype: //2x2 contingency table analysis
                analyze_dichotom<2,2>(par, myout, trait, &study);
                break;
        }
    }

} // namespace Permory

#endif // include guard


