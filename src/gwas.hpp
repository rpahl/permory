// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_individual_hpp
#define permory_individual_hpp

#include <string>
#include <vector>

#include "config.hpp"
#include "individual.hpp"
#include "locus.hpp"

namespace Permory 
{
    class Study {
        public:
            // Iterator pass through
            typedef std::set<Individual>::iterator iterator;
            typedef std::set<Individual>::const_iterator const_iterator;
            iterator begin() { return sample_.begin(); }
            iterator end() { return sample_.end(); }
            const_iterator begin() const { return sample_.begin(); }
            const_iterator end() const { return sample_.end(); }

            // Ctor
            Study() {}
            Study(const std::set<Individual>& sample) 
                : sample_(sample)
            {}

            // Inspection
            const std::set<Individual>& sample() const { return sample_; }

            // Modification
            void recruit(const Individual& i) { sample_.insert(i); }
            void exclude(const Individual& i) { sample_.erase(i); }

        private:
            std::set<Individual> sample_;
    };

    // Genome-wide association study
    class Gwas {
        public:
            // Ctor

            // Inspection
            size_t sample_size() const { return sample_.sample().size(); }
            size_t m() const { return loci_.size(); }

        private:
            Study sample_;
            std::vector<Locus> loci_;
    };
} // namespace Permory

#endif // include guard


