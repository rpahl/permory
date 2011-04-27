// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_analysis_hpp
#define permory_analysis_hpp

#include <set>
#include <string>
#include <vector>
#include <iomanip>

#include <boost/progress.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "detail/config.hpp"
#include "detail/parameter.hpp"
#include "detail/exception.hpp"
#include "gwas.hpp"
#include "locusdata.hpp"
#include "locus_filter.hpp"
#include "io/output.hpp"
#include "permutation/perm.hpp"
#include "read_phenotype_data.hpp"
#include "read_locus_data.hpp"
#include "result_output.hpp"
#include "statistical/dichotom.hpp"
#include "statistical/quantitative.hpp"
#include "statistical/pvalue.hpp"

namespace Permory { namespace gwas {

    class Analyzer {
        public:
            // Ctor
            Analyzer(
                    detail::Parameter* par, io::Myout& out,
                    Gwas* study, std::set<char> the_domain=std::set<char>())
                : par_(par), out_(out), study_(study), domain_(the_domain)
                { init_filters(); }

            // Dtor
            virtual ~Analyzer() { }

            // Modifiers
            template<class S, class T> void analyze_dichotom();

        protected:
            virtual void output_results(std::deque<double>& tmax);
            virtual void init_filters(); // Define locus filter (e.g. maf filter)

            detail::Parameter* par_;
            io::Myout& out_;
            Gwas* study_;
            std::set<char> domain_;

            boost::ptr_vector<Locus_filter> locus_filters_;
    };

    // Factories
    // ========================================================================
    class Abstract_analyzer_factory {
        public:
            // Ctor
            // Dtor
            virtual ~Abstract_analyzer_factory() { }

            boost::shared_ptr<Analyzer> operator()(
                    detail::Parameter* par, io::Myout& out, Gwas* study,
                    std::set<char> the_domain=std::set<char>())
            {
                return get_analyzer(par, out, study, the_domain);
            }

        private:
            virtual boost::shared_ptr<Analyzer> get_analyzer(
                    detail::Parameter* par, io::Myout& out, Gwas* study,
                    std::set<char> the_domain=std::set<char>()) = 0;
    };
    class Default_analyzer_factory : public Abstract_analyzer_factory {
        public:
            // Ctor
            // Dtor
            virtual ~Default_analyzer_factory() { }

        private:
            boost::shared_ptr<Analyzer> get_analyzer(
                    detail::Parameter* par, io::Myout& out, Gwas* study,
                    std::set<char> the_domain=std::set<char>())
            {
                boost::shared_ptr<Analyzer> ptr(new Analyzer(par, out, study, the_domain));
                return ptr;
            }
    };

    // Analyzer implementation
    // ========================================================================

    //
    //  Compute test statistics and perform permutation test
    template<class S, class T> void Analyzer::analyze_dichotom()
    {
        using namespace std;
        using namespace boost;
        using namespace io;

        // Prepare progress bar
        size_t m = study_->m();
        size_t perm_todo = par_->nperm_total;        //remaining permutations
        scoped_ptr<progress_display> pprogress;
        bool show_progress = not par_->quiet;
        if (show_progress) {
            double d = ceil(double(perm_todo)/double(par_->nperm_block));
            pprogress.reset(new progress_display(m*size_t(d), std::cout, "","",""));
        }

        vector<Individual> trait(study_->ind_begin(), study_->ind_end());
        // If allelic analyis, double trait status for each of the two alleles
        if (par_->marker_type == allelic) {
            vector<Individual> tmp;
            tmp.reserve(trait.size()*2);
            BOOST_FOREACH(Individual b, trait) {
                tmp.push_back(b);
                tmp.push_back(b);
            }
            trait.swap(tmp);
        }

        deque<double> tmax;     //holds the maximum test statistic per permutation
        S stat(*par_, trait.begin(), trait.end());   //computes all statistic stuff
        permutation::Permutation pp(par_->seed);     //does the shuffling/permutation
        Gwas::iterator itLocus = study_->begin();    //points to current locus

        bool isFirstRound = true;
        while (perm_todo > 0) {
            // By default analysis is done in blocks of 10000 permutations. In
            // each block each marker is analyzed one by one
            size_t nperm = par_->nperm_block;
            if (nperm > perm_todo) {
                nperm = perm_todo;
            }
            stat.renew_permutations(&pp, nperm, par_->tail_size); //fresh random numbers
            itLocus = study_->begin(); 

            BOOST_FOREACH(string fn, par_->fn_marker_data) { 
                Locus_data_reader<char> loc_reader(fn, par_->undef_allele_code);

                while (loc_reader.hasData()) {
                    std::vector<char> v;
                    loc_reader.get_next(v);
                    Locus_data<char> data(v, par_->undef_allele_code);
                    if (data.size() == 0) {
                        continue;
                    }

                    if (par_->marker_type == genotype) {
                        if (data.size() == 2*trait.size()) {
                            data = data.condense_alleles_to_genotypes(2);
                        }
                        // Ensure that all possible genotypes appear in the domain
                        data.add_to_domain(domain_);
                    }

                    if (data.size() != trait.size()) {
                        // Provide some more information in case of this error
                        // as it may be hard to spot in large data sets
                        throw Data_length_mismatch_error(
                                itLocus->id(), trait.size(), data.size());
                    }

                    // The non-permutation stuff needs only to be done once 
                    if (isFirstRound) { 
                        itLocus->set_polymorph(data.isPolymorph()); 
                        itLocus->set_maf(data.maf(par_->marker_type)); 

                        bool hasPassed = true;
                        ptr_vector<Locus_filter>::iterator itFilter = 
                            locus_filters_.begin();
                        for (; itFilter != locus_filters_.end(); ++itFilter) {
                            if (not (*itFilter)(*itLocus)) {
                                hasPassed = false;
                                break;
                            }
                        }

                        if (hasPassed) {
                            itLocus->add_test_stats(stat.test(data)); 
                        }
                    }
                    if (itLocus->hasTeststat()) {
                        stat.permutation_test(data);   
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

        output_results(tmax);
    }

    void Analyzer::output_results(std::deque<double>& tmax)
    {
        using namespace std;
        using namespace io;
        using namespace statistic;

        #define TIME(X, Y) t.restart(); Y; out_ << all << stdpre << X << t.elapsed() << " s" << endl;
        boost::timer t;

        out_ << endl;
        result_to_console(par_, out_, *study_);

        // Get tmax of original data and use permutation tmax to derive p-values
        out_ << normal << stdpre << "Writing results." << endl;
        deque<double> t_orig(study_->m());

        TIME("Runtime transform all: ",
            transform(study_->begin(), study_->end(),
                    t_orig.begin(), mem_fun_ref(&Locus::tmax)));

        TIME("Runtime single_step_counts all: ",
            deque<size_t> counts = single_step_counts(t_orig, &tmax));
        std::string fn = par_->out_prefix;
        fn.append(".all");
        TIME("Runtime result_to_file all: ",
            result_to_file(par_, *study_, counts, fn));

        // The same but this time just for the top p-values
        TIME("Runtime sort top: ",
            sort(study_->begin(), study_->end(), Locus_tmax_greater()));
        t_orig.resize(par_->ntop);
        TIME("Runtime transform top: ",
            transform(study_->begin(), study_->begin()+par_->ntop,
                    t_orig.begin(), mem_fun_ref(&Locus::tmax)));
        TIME("Runtime single_step_counts top: ",
            counts = single_step_counts(t_orig, &tmax));
        fn = par_->out_prefix;
        fn.append(".top");
        TIME("Runtime result_to_file top: ",
            result_to_file(par_, *study_, counts, fn));
    }

    void Analyzer::init_filters()
    {
        locus_filters_.push_back(new Maf_filter(par_->min_maf, par_->max_maf));
        locus_filters_.push_back(new Polymorph_filter());
    }

    // None-Class Functions
    // ========================================================================
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
            myout << normal << stdpre << "Reading trait from file..." << endl;
            Datafile_format dff = detect_phenotype_data_format(fn);
            par->phenotype_data_format = dff;
            myout << verbose << indent(4) << "`" << par->fn_trait << "' -> assuming " <<
                "file format " << datafile_format_to_string(dff) << endl;
            read_individuals(*par, fn, &v);
        }
        else {
            if (not (par->val_type == Record::dichotomous)) {
                throw invalid_argument(
                        "Need trait file for non dichotomous phenotypes.");
            }
            // No file so create individuals by specified numbers
            size_t n = par->ncase + par->ncontrol;
            if (n == 0) {
                throw runtime_error("Sample size must not be 0");
            }
            myout << normal << stdpre << "Creating trait via given numbers..." << endl;
            v = case_control_sample(n, double(par->ncase)/double(n)); //individual.hpp
        }
        return v;
    }

    //
    // Scan all marker data files and store loci
    void scan_loci(Gwas* the_study, detail::Parameter* par, io::Myout& myout)
    {
        using namespace std;
        using namespace boost;
        using namespace io;
        myout << normal << stdpre << "Scanning marker data..." << endl;
        BOOST_FOREACH(string fn, par->fn_marker_data) {
            Datafile_format dff = detect_marker_data_format(fn, par->undef_allele_code);
            myout << verbose << indent(4) << "`" << fn << "' -> assuming file format " << 
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

    void gwas_analysis(detail::Parameter* par, io::Myout& myout,
                        Abstract_analyzer_factory& factory)
    {
        using namespace std;
        using namespace io;

        // Phenotype data (either read in or create)
        Gwas study(create_study_sample(par, myout)); 
        if (par->val_type == Record::dichotomous) {
            if (not par->fn_trait.empty()) {
                myout << normal << stdpre << "Found ";
            }
            else {
                myout << normal << stdpre << "Created ";
            }
            myout << study.ncase() << " cases and " << study.ncontrol() <<
                " controls." << endl;
        }

        // Loci information read in
        scan_loci(&study, par, myout);
        size_t m = study.m();
        if (study.m() == 0) {
            throw runtime_error("No marker available.");
        }
        myout << normal << stdpre << m << " markers found." << endl << endl; 

        myout << normal << stdpre << "Starting analysis..." << endl;
        std::set<char> data_domain;
        switch (par->marker_type) {
            case genotype:  //2x3 contingency table analysis
            {
                data_domain.insert('0');
                data_domain.insert('1');
                data_domain.insert('2');
                boost::shared_ptr<Analyzer> analyzer = factory(par, myout, &study, data_domain);
                if (par->val_type == Record::continuous) {
                    analyzer->analyze_dichotom<statistic::Quantitative<3>, double>();
                }
                else {
                    analyzer->analyze_dichotom<statistic::Dichotom<2,3>, bool>();
                }
                break;
            }
            case allelic: //2x2 contingency table analysis
                boost::shared_ptr<Analyzer> analyzer = factory(par, myout, &study);
                analyzer->analyze_dichotom<statistic::Dichotom<2,2>, bool>();
                break;
        }
    }

} // namespace gwas
} // namespace Permory

#endif // include guard


