// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_gwas_hpp
#define permory_gwas_hpp

#include <set>
#include <string>
#include <vector>

#include <boost/progress.hpp>

#include "detail/config.hpp"
#include "detail/parameter.hpp"
#include "individual.hpp"
#include "locus.hpp"
#include "locusdata.hpp"
#include "io/output.hpp"
#include "permutation/perm.hpp"
#include "read_phenotype_data.hpp"
#include "read_locus_data.hpp"
#include "statistical/dichotom.hpp"
#include "statistical/pvalue.hpp"

namespace Permory 
{
    //using namespace ::boost::multi_index;
    //namespace bmi = boost::multi_index;

    // Genome-wide association study
    class Gwas {
        public:
            //typedef typename boost::ptr_deque<data_type>::iterator iterator;
            typedef std::deque<Locus>::iterator iterator;
            typedef std::deque<Locus>::const_iterator const_iterator;
            typedef std::vector<Individual>::const_iterator const_inderator;

            // Iterator pass through
            iterator begin() { return loci_.begin(); }
            iterator end() { return loci_.end(); }
            const_iterator begin() const { return loci_.begin(); }
            const_iterator end() const { return loci_.end(); }
            const_inderator ind_begin() const { return ind_.begin(); }
            const_inderator ind_end() const { return ind_.end(); }

            // Ctor
            Gwas(const std::vector<Individual>& ind);   
            Gwas(const_inderator start, const_inderator end);   

            // Inspection
            size_t m() const { return loci_.size(); }
            size_t ncase() const;
            size_t ncontrol() const { return (sample_size() - ncase()); }
            size_t sample_size() const { return ind_.size(); }
            bool has_unique_loci() const;

            // Modification
            std::deque<Locus>* pointer_to_loci() { return &loci_; }
            void add_loci(const_iterator start, const_iterator end);

        private:
            std::vector<Individual> ind_;   //recruited individuals
            std::deque<Locus> loci_;
    };

    // Gwas implementation
    // ========================================================================
    inline Gwas::Gwas(
            const std::vector<Individual>& ind)
       : ind_(ind)
    {
        if (ind_.size() < 2) {
            throw std::invalid_argument("Gwas must have at least 2 individuals.");
        }
    }

    inline Gwas::Gwas(
            const_inderator start, const_inderator end)
        : ind_(start, end)
    {
        if (ind_.size() < 2) {
            throw std::invalid_argument("Gwas must have at least 2 individuals.");
        }
    }

    inline size_t Gwas::ncase() const
    {
        return std::count_if(ind_.begin(), ind_.end(), 
                    std::mem_fun_ref(&Individual::isAffected));
    }

    inline bool Gwas::has_unique_loci() const
    {
        std::deque<Locus> loc(loci_);
        std::sort(loc.begin(), loc.end());
        return (std::adjacent_find(loc.begin(), loc.end()) == loc.end());
    }

    inline void Gwas::add_loci(
            const_iterator start, const_iterator end) 
    {
        std::copy(start, end, std::back_inserter(loci_));
    }


    std::vector<Individual> create_study_sample(
            detail::Parameter* par, io::Out_log& myout)
    {
        using namespace std;
        using namespace boost;
        using namespace detail;
        using namespace io;
        vector<Individual> v;
        string fn = par->fn_trait;
        if (not fn.empty()) {
            myout << normal << "...Read trait from file `" << par->fn_trait << "'" << myendl;

            // Read from file
            par->phenotype_data_format = detect_phenotype_data_format(fn);
            myout << verbose << "...Detected file format of trait file: " << 
                detail::datafile_format_to_string(par->phenotype_data_format) << myendl;
            read_individuals(*par, fn, &v);
        }
        else {
            // No file so create individuals by specified numbers
            size_t n = par->ncase + par->ncontrol;
            if (n == 0) {
                throw runtime_error("Sample size must not be 0");
            }
            myout << normal << "...Creating trait by number..." << myendl;
            v = case_control_sample(n, double(par->ncase)/double(n)); //individual.hpp
        }
        return v;
    }

    template<uint K, uint L> void analyze_dichotom_haplotype(
            detail::Parameter* par, io::Out_log& myout)
    {
        //TODO
    }

    template<uint K, uint L> void analyze_dichotom_genotype(
            detail::Parameter* par, io::Out_log& myout)
    {
        using namespace std;
        using namespace boost;
        using namespace io;
        std::set<uint> genotype_domain;
        for (uint i=0; i<L; i++) {
            genotype_domain.insert(i);
        }

        //
        // Phenotype data (either read in or create)
        //
        vector<Individual> v = create_study_sample(par, myout);
        Gwas study(v); 
        if (not par->fn_trait.empty()) {
            myout << normal << "...Found ";
        }
        else {
            myout << normal << "...Created ";
        }
        myout << study.ncase() << " cases and " << study.ncontrol() << 
            " controls." << myendl;

        //
        // Loci information read in
        //
        myout << normal << "...Scanning marker data - please be patient..." << myendl;
        BOOST_FOREACH(string fn, par->fn_marker_data) {
            Datafile_format dff = detect_marker_data_format(fn, par->undef_allele_code);
            myout << verbose << "\t" << fn << " - file format: " << 
                detail::datafile_format_to_string(dff);
            if (dff == unknown) {
                par->fn_bad_files.insert(fn);
                myout << " - will be ignored." << myendl;
                continue;
            }
            myout << myendl;
            read_loci(dff, fn, study.pointer_to_loci());
        }
        if (not par->fn_bad_files.empty()) {
            BOOST_FOREACH(std::string s, par->fn_bad_files) {
                par->fn_marker_data.erase(s);
            }
            cerr << "! Warning: " << par->fn_bad_files.size() << " data file(s) " <<
                "will NOT be processed due to unrecognized data format." << endl;
        }
        size_t m = study.m();
        if (m == 0) {
            throw runtime_error("No marker loci found.");
        }
        myout << normal << "...Found " << m << " loci." << myendl << myendl; 

        // Create dichotomous trait/phenotype data
        vector<bool> trait(study.sample_size());
        transform(study.ind_begin(), study.ind_end(), trait.begin(),
                mem_fun_ref(&Individual::isAffected));

        //
        //  Compute test statistics of original data and permutation test
        //
        myout << normal << "..." << "Starting analysis ..." << myendl;
        size_t m_maf = 0;   //number of markers passing maf criterion
        size_t m_valid = 0; //number of markers passing every filter
        deque<double> tmax; //permutation Tmax statistics
        statistic::Dichotom<K, L> stat(*par, trait); //handles the statistic stuff
        Gwas::iterator itLocus = study.begin(); //points to current locus

        size_t perm_todo = par->nperm_total;
        permutation::Permutation pp(par->seed);
        double d = ceil(double(perm_todo)/double(par->nperm_block));
        scoped_ptr<progress_display> pprogress;
        bool show_progress = not par->quiet;
        if (show_progress) {
            pprogress.reset(new progress_display(m*size_t(d), std::cout, "","",""));
        }

        bool isFirstRound = true;
        while (perm_todo > 0) {
            size_t nperm = par->nperm_block;
            if (nperm > perm_todo) {
                nperm = perm_todo;
            }
            stat.renew_permutations(&pp, nperm, par->tail_size); //fresh random numbers
            itLocus = study.begin(); 

            BOOST_FOREACH(string fn, par->fn_marker_data) { 
                Datafile_format dff = detect_marker_data_format(fn, par->undef_allele_code);
                Locus_data_reader<char> loc_reader(fn, par->undef_allele_code);

                while (loc_reader.hasData()) {
                    Locus_data<char> data(loc_reader.get_next(), par->undef_allele_code);
                    scoped_ptr<Locus_data<uint> > data_ptr;

                    // Re-format data depending on data file format
                    if (dff == permory_data || dff == slide) {   
                        data_ptr.reset(data.as_numeric());
                    }
                    else if (dff == presto || dff == plink_tped) { 
                        data_ptr.reset(data.merge_alleles_to_genotypes(2));
                    }
                    //ensure that all possible genotypes appear in the domain
                    data_ptr->add_to_domain(genotype_domain); 

                    // Stuff of the orginal data needs only be processed once
                    if (isFirstRound) { 
                        bool isPoly = data_ptr->isPolymorph();
                        itLocus->set_polymorph(isPoly); 
                        bool mafOk = (data_ptr->maf(par->gen_type) > par->min_maf);
                        m_maf += mafOk;
                        if (isPoly && mafOk) {
                            m_valid++;
                            itLocus->add_test_stats(stat.test(*data_ptr)); 
                        }
                    }
                    if (itLocus->hasTeststat() && perm_todo > 0) {
                        stat.permutation_test(*data_ptr);   //permute and compute
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
        myout << myendl;
        size_t m_poly = count_if(study.begin(), study.end(), 
                mem_fun_ref(&Locus::isPolymorph));
        myout << normal << "..." << m - m_poly << " loci were non-polymorphic." << myendl; 
        myout << "..." << m_poly - m_maf << " polymorphic loci had minor " <<
            "allele frequency lower than " << par->min_maf << myendl;
        myout << normal << "..." << m_valid << " of " << m << " loci remained " <<
            "for analysis." << myendl << myendl;

        deque<double> t_orig(study.m());
        transform(study.begin(), study.end(), t_orig.begin(), 
                mem_fun_ref(&Locus::tmax));
        sort(t_orig.begin(), t_orig.end(), greater<double>());

        //detail::print_seq(t_orig.begin(), t_orig.end());
        //cout << endl;
        //detail::print_seq(tmax.begin(), tmax.end());
        cout << endl;
        //sort(study.begin(), study.end(), mem_fun_ref(&Locus::tmax));
        //detail::print_seq(tmax.begin(), tmax.end());
        //deque<size_t> counts = statistic::single_step_counts(t_orig, tmax);
        //detail::print_seq(counts.begin(), counts.end());
        cout << endl;
        deque<double> pvals = statistic::single_step_pvalues(t_orig, tmax);
        size_t offset = std::min(pvals.size(), size_t(50));
        detail::print_seq(pvals.begin(), pvals.begin()+offset);
        cout << endl;
    }

    void gwas_analysis(detail::Parameter* par, io::Out_log& myout)
    {
        switch (par->gen_type) {
            case genotype:
                analyze_dichotom_genotype<2,3>(par, myout);
                break;
            case haplotype:
                analyze_dichotom_haplotype<2,2>(par, myout);
                break;
        }
    }

} // namespace Permory

#endif // include guard

/*
   bmi::multi_index_container<
   Locus,
   indexed_by<
   ordered_unique<identity<Locus> >,
   ordered_non_unique<
   BOOST_MULTI_INDEX_CONST_MEM_FUN(Locus, double, tmax)
   >
   > 
   > loci_;
   */


