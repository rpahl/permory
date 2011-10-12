// Copyright (c) 2011 Volker Steiß
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#define PERMORY_TEST gwas_test

#include <fstream>

#include "detail/parameter.hpp"
#include "gwas/read_phenotype_data.hpp"
#include "test.hpp"

using namespace std;
using namespace boost;
using namespace unit_test;
using namespace Permory;
using namespace Permory::io;
using namespace Permory::detail;
using namespace Permory::gwas;

void check_individuals(const vector<Individual>& individuals) {
    BOOST_CHECK_EQUAL(size_t(5), individuals.size());

    vector<Individual>::const_iterator it = individuals.begin();
    BOOST_CHECK_EQUAL(0, it++->begin()->val);
    BOOST_CHECK_EQUAL(-1, it++->begin()->val);
    BOOST_CHECK_EQUAL(-0.42, it++->begin()->val);
    BOOST_CHECK_EQUAL(3.721, it++->begin()->val);
    BOOST_CHECK_EQUAL(423, it++->begin()->val);
}

void read_individuals_from_tfam_test()
{
    const string filename = "test/data/pheno_quant_test.tfam";
    Parameter par;
    par.phenotype_domain = Record::continuous;
    par.phenotype_data_format = plink_tfam;
    vector<Individual> individuals;

    read_individuals_from_tfam(par, filename, &individuals);

    check_individuals(individuals);
}

void read_individuals_test() {
    const string filename = "test/data/pheno_quant_test.bgl";
    Parameter par;
    par.phenotype_domain = Record::continuous;
    par.phenotype_data_format = presto;
    vector<Individual> individuals;

    read_individuals(par, filename, &individuals);

    check_individuals(individuals);
}

void determine_phenotype_domain_test() {
    { // test compact format with dichotomous phenotypes
        const string filename = "test/data/pheno_binary_test.compact";
        Parameter par;
        par.phenotype_data_format = compact;
        Record::Value_type expected = Record::dichotomous;

        BOOST_CHECK_EQUAL(determine_phenotype_domain(par, filename), expected);
    }

    { // test plink_tfam format with quantitative phenotypes
        const string filename = "test/data/pheno_quant_test.tfam";
        Parameter par;
        par.phenotype_data_format = plink_tfam;
        Record::Value_type expected = Record::continuous;

        BOOST_CHECK_EQUAL(determine_phenotype_domain(par, filename), expected);
    }
    { // test plink_tfam format with dichotomous phenotypes
        const string filename = "test/data/pheno_binary_test.tfam";
        Parameter par;
        par.phenotype_data_format = plink_tfam;
        Record::Value_type expected = Record::dichotomous;

        BOOST_CHECK_EQUAL(determine_phenotype_domain(par, filename), expected);

        const string filename2 = "test/data/pheno_binary_test2.tfam";
        BOOST_CHECK_EQUAL(determine_phenotype_domain(par, filename2), expected);
    }
    { // test plink_tfam format with undefined phenotypes
        const string filename = "test/data/pheno_binary_undef_test.tfam";
        Parameter par;
        par.phenotype_data_format = plink_tfam;
        par.undef_phenotype_code = "?";

        BOOST_CHECK_THROW(determine_phenotype_domain(par, filename),
                          std::domain_error);
    }

    { // test presto format with quantitative phenotypes
        const string filename = "test/data/pheno_quant_test.bgl";
        Parameter par;
        par.phenotype_data_format = presto;
        Record::Value_type expected = Record::continuous;

        BOOST_CHECK_EQUAL(determine_phenotype_domain(par, filename), expected);
    }
    { // test presto format with dichotomous phenotypes
        const string filename = "test/data/pheno_binary_test.bgl";
        Parameter par;
        par.phenotype_data_format = presto;
        Record::Value_type expected = Record::dichotomous;

        BOOST_CHECK_EQUAL(determine_phenotype_domain(par, filename), expected);
    }
}

test_suite* init_unit_test_suite( int argc, char* argv[] )
{
    test_suite *test = BOOST_TEST_SUITE("Functions and classes from src/gwas");

    test->add(BOOST_TEST_CASE(&read_individuals_from_tfam_test));
    test->add(BOOST_TEST_CASE(&read_individuals_test));
    test->add(BOOST_TEST_CASE(&determine_phenotype_domain_test));

    return test;
}

