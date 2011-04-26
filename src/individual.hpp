// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_individual_hpp
#define permory_individual_hpp

#include <string>
#include <vector>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include "detail/config.hpp"
#include "gwas/locus.hpp"

namespace Permory 
{
    // Data record for an individual at a particular experiment/study
    struct Record {
            enum Value_type {undefined=0, dichotomous, continuous};

            Record(double d=0.0, Value_type type = continuous)
                : val(d), theType(type)
            {}
            bool operator<(const Record& r) const { return val < r.val; }
            bool operator==(const Record& r) const { return val == r.val; }
            bool as_dichotom() const { return (val != 0); }
            bool is_defined() const { return theType != undefined; }

            double val;
            Value_type theType;

        private:
            // serialization stuff
            friend class boost::serialization::access;
            template<class Archive>
            void serialize(Archive & ar, const unsigned int version)
            {
                ar & val;
                ar & theType;
            }
    };

    //
    // An individual as part of a clinical trial or association study, etc...
    //
    class Individual {
        public:
            enum Sex {nosex=0, male, female};

            typedef std::vector<Record>::iterator iterator;
            typedef std::vector<Record>::const_iterator const_iterator;
            iterator begin() { return r_.begin(); }
            iterator end() { return r_.end(); }
            const_iterator begin() const { return r_.begin(); }
            const_iterator end() const { return r_.end(); }

            // Ctors
            Individual(size_t id, std::string name="", Sex sex=nosex, double c=0) 
                : id_(id), name_(name), sex_(sex), cost_(c) 
            {}

            // Inspection
            bool operator<(const Individual& i) const { return id_ < i.id(); }
            bool operator>(const Individual& i) const { return id_ > i.id(); }
            bool operator==(const Individual& i) const { return id_ == i.id(); }
            bool hasRecord() const { return (!r_.empty()); }
            bool isAffected() const { return hasRecord() && begin()->as_dichotom(); }
            size_t id() const { return id_; }
            size_t nrecord() const { return r_.size(); }
            double costs() const { return cost_; }
            std::string name() const { return name_; }
            Sex sex() const { return sex_; }

            // Modification
            void set_id(size_t id) { id_ = id; }
            void add_measurement(Record& r) { r_.push_back(r); }
            void set_cost(double c) { cost_ = c; }

        private:
            size_t id_;         //individual ID
            std::string name_;  //the individual's name
            Sex sex_;
            std::vector<Record> r_;
            double cost_; //for example, cost of recruitment or analysis

            // serialization stuff
            friend class boost::serialization::access;

            Individual() {}

            template<class Archive>
            void serialize(Archive & ar, const unsigned int version)
            {
                ar & id_;
                ar & name_;
                ar & sex_;
                ar & r_;
                ar & cost_;
            }
    };

    //
    // Sample of individuals with dichotomous phenotypes consisting of 
    // (1-r) controls and r cases  
    std::vector<Individual> case_control_sample(
            size_t n, double r=0.5, 
            double c=0) //costs
    {
        std::vector<Individual> v;
        if (r <= 0 || r >= 1)
            throw std::invalid_argument("The ratio must be in (0,1).");
        v.reserve(n);

        Record rec0(0.0, Record::dichotomous);
        Record rec1(1.0, Record::dichotomous);

        // Fill first (1-r) with controls and rest with cases
        for (size_t i=0; i<n; i++) {
            Individual ind(i, "", Individual::nosex, c);
            if (double(i)/double(n) < r) {
                ind.add_measurement(rec0);
            }
            else {
                ind.add_measurement(rec1);
            }
            v.push_back(ind);
        }
        return v;
    }

    //
    // Sample of individuals with continuous phenotypes
    std::vector<Individual> contiunous_sample(
            const std::vector<double>& phenotypes, double c=0)
    {
        std::vector<Individual> v;
        v.reserve(phenotypes.size());
        for (size_t i=0; i<phenotypes.size(); i++) {
            Individual ind(i, "", Individual::nosex, c);
            Record rec(phenotypes[i], Record::continuous);
            ind.add_measurement(rec);
            v.push_back(ind);
        }
        return v;
    }

} // namespace Permory

// Boost Serialization API Version Information
// ========================================================================
BOOST_CLASS_VERSION(Permory::Record, 1)
BOOST_CLASS_VERSION(Permory::Individual, 1)

#endif // include guard


