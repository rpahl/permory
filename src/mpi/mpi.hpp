// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_mpi_hpp
#define permory_mpi_hpp

#include <set>
#include <deque>
#include <vector>

//#include <boost/progress.hpp>
//#include <boost/ptr_container/ptr_vector.hpp>
//#include <boost/shared_ptr.hpp>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/progress.hpp>   //timer

#include <boost/mpi/operations.hpp>

#include "detail/config.hpp"
#include "detail/parameter.hpp"
#include "detail/functors.hpp"
#include "gwas/gwas.hpp"
#include "gwas/analysis.hpp"
#include "io/output.hpp"

// Make deque_concat() commulative to gain more speed.
namespace boost { namespace mpi {

  template<>
  struct is_commutative<Permory::detail::deque_concat<double>, std::deque<double> >
    : mpl::true_ { };
} } // end namespace boost::mpi

namespace Permory { namespace gwas {
    namespace mpi = boost::mpi;

    class Mpi_analyzer : public Analyzer {
        public:
            // Ctor
            Mpi_analyzer(
                    mpi::environment& env, mpi::communicator& world,
                    detail::Parameter* par, io::Myout& out,
                    const std::vector<bool>& trait, Gwas* study,
                    std::set<char> the_domain=std::set<char>())
                : Analyzer(par, out, trait, study, the_domain),
                  env_(env), world_(world)
            {
                // memorize old nperm_total and set per process nperm_total
                orig_nperm_total_ = par_->nperm_total;
                par_->nperm_total = nperm_per_process();

                // different seed for all processes
                par_->seed += world_.rank();

                if (world_.rank() > 0) {
                    par_->quiet = true;
                    par_->verbose = false;
                }
            }

            // Dtor
            virtual ~Mpi_analyzer() { }

            // Modifiers
            template<uint K, uint L> void analyze_dichotom();

        protected:
            size_t nperm_per_process() const;

            virtual void output_results(std::deque<double> tmax);

        private:
            mpi::environment& env_;
            mpi::communicator& world_;

            size_t orig_nperm_total_;
    };

    class Mpi_analyzer_factory : public Abstract_analyzer_factory {
        public:
            // Ctor
            Mpi_analyzer_factory(int& argc, char **& argv)
                : Abstract_analyzer_factory(argc, argv), env_(argc, argv)
                { }

            // Dtor
            virtual ~Mpi_analyzer_factory() { }

        private:
            boost::shared_ptr<Analyzer> get_analyzer(
                    detail::Parameter* par, io::Myout& out,
                    const std::vector<bool>& trait, Gwas* study,
                    std::set<char> the_domain=std::set<char>())
            {
                boost::shared_ptr<Analyzer> ptr(new Mpi_analyzer(env_, world_, par, out, trait, study, the_domain));
                return ptr;
            }

            mpi::environment env_;
            mpi::communicator world_;
    };

    // Analyzer implementation
    // ========================================================================

    size_t Mpi_analyzer::nperm_per_process() const
    {
        size_t per_process = par_->nperm_total / world_.size();
        return world_.rank() == 0
                     ? par_->nperm_total - size_t(per_process * (world_.size() - 1))
                     : per_process;
    }

    //
    //  Compute test statistics and perform permutation test
    void Mpi_analyzer::output_results(std::deque<double> tmax)
    {
        using namespace std;
        using namespace boost::mpi;
        using namespace detail;

        if (world_.rank() == 0) {
            boost::timer t;
            deque<double> tmax_result;
            reduce(world_, tmax, tmax_result, deque_concat<double>(), 0);
            out_ << io::stdpre << "Runtime reduce: " << t.elapsed() << " s" << endl;
            t.restart();
            par_->nperm_total = orig_nperm_total_;  // reset nperm_total for correct output calculations
            Analyzer::output_results(tmax_result);
            out_ << io::stdpre << "Runtime output: " << t.elapsed() << " s" << endl;
        }
        else {
            reduce(world_, tmax, deque_concat<double>(), 0);
        }
    }

} // namespace gwas
} // namespace Permory

#endif // include guard


