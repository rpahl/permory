// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_individual_hpp
#define permory_individual_hpp

#include <string>
#include <vector>

#include "config.hpp"
#include "locus.hpp"

namespace Permory 
{
    // forward declarations
    class Measurement;

    class Individual { 
        public:
            enum Sex {nosex=0, male, female};

            // Ctor
            Individual(size_t id, 
                    std::string name="", 
                    Sex sex=nosex) 
                : id_(id), name_(name), sex_(sex) 
            { }

            // Inspection
            size_t id() const { return id_; }
            std::string name() const { return name_; }
            Sex sex() const { return sex_; }
            const std::vector<Measurement>& measurements() const { return v_; }
            bool operator<(const Individual& i) const { return id_ < i.id(); }
            bool operator>(const Individual& i) const { return id_ > i.id(); }
            bool operator==(const Individual& i) const { return id_ == i.id(); }

            // Modification
            void add_measurement(const Measurement& m) { v_.push_back(m); }

        private:
            size_t id_;  //individual ID
            std::string name_;    
            Sex sex_;
            std::vector<Measurement> v_;
    };


    // Data measured for an individual at a particular experiment/study
    struct Measurement {
        Measurement(
                Individual* i, 
                bool aff=false, 
                double phen=0.0)
            : pi(i), affected(aff), phenotype(phen)
        {}

        bool affected;
        double phenotype;
        Individual* pi;
    };

} // namespace Permory

#endif // include guard


