// Copyright (c) 2010-2011 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_enums_hpp
#define permory_detail_enums_hpp

#include <boost/bimap.hpp>

namespace Permory { namespace detail {

    // Supported data formats
    enum datafile_format{
        unknown=0, compact, slide, presto, plink_tfam, plink_tped };

    enum Marker_type {allelic, genotype};

    // Supported statistical tests
    enum Test_type { trend, trend_extended, chisq };

    // Variance estimation
    enum Var_estimate {pooled=0, controls=1, separately=2};

    // Verbosity level
    enum Verbosity {all=0, verbose=1, normal=2, muted=3};

    //
    // Provides conversion of enumeration types to string and vice versa
    //
    class Enum_converter {
        public:
            typedef boost::bimap<datafile_format, std::string> dff_bimap;
            typedef boost::bimap<Marker_type, std::string> marker_bimap;
            typedef boost::bimap<Test_type, std::string> test_bimap;
            typedef boost::bimap<Var_estimate, std::string> ve_bimap;
            typedef boost::bimap<Verbosity, std::string> verb_bimap;

            // Ctor
            Enum_converter();

            template<class T> std::string key_to_string(T);
            template<class T> T string_to_key(const std::string&);

        private:
            // Member variables
            dff_bimap file_formats_;
            marker_bimap markers_;
            test_bimap tests_;
            ve_bimap var_estimates_;
            verb_bimap verbs_;

            // Member functions
            template<class T> 
                std::string key_to_string(const boost::bimap<T, std::string>&, T);
            template<class T> T 
                string_to_key(const boost::bimap<T, std::string>&, const std::string&);
            template<class enum_t> 
                std::string type2str();
    };

    // Enum_converter implementation
    // ========================================================================
    Enum_converter::Enum_converter()
    {
        file_formats_.insert(dff_bimap::value_type(compact, "compact (*.comp)"));
        file_formats_.insert(dff_bimap::value_type(plink_tfam, "PLINK (*.tfam)"));
        file_formats_.insert(dff_bimap::value_type(plink_tped, "PLINK (*.tped)"));
        file_formats_.insert(dff_bimap::value_type(presto, "PRESTO (*.bgl)"));
        file_formats_.insert(dff_bimap::value_type(slide, "SLIDE (*.slide)"));
        file_formats_.insert(dff_bimap::value_type(unknown, "unknown"));

        markers_.insert(marker_bimap::value_type(allelic, "allelic"));
        markers_.insert(marker_bimap::value_type(genotype, "genotype"));

        tests_.insert(test_bimap::value_type(trend, "T_trend"));
        tests_.insert(test_bimap::value_type(trend_extended, "T_trend_ext"));
        tests_.insert(test_bimap::value_type(chisq, "T_chi^2"));

        var_estimates_.insert(ve_bimap::value_type(pooled, "pooled"));
        var_estimates_.insert(ve_bimap::value_type(controls, "only controls"));
        var_estimates_.insert(ve_bimap::value_type(separately, "separately cases and controls"));

        verbs_.insert(verb_bimap::value_type(all, "all"));
        verbs_.insert(verb_bimap::value_type(verbose, "verbose"));
        verbs_.insert(verb_bimap::value_type(normal, "normal"));
        verbs_.insert(verb_bimap::value_type(muted, "muted"));
    }

    // 
    // Conversion implementations
    template<class T> std::string Enum_converter::key_to_string(
            const boost::bimap<T, std::string>& map, T key)
    {
        typedef boost::bimap<T, std::string> bimap_t;
        typename bimap_t::left_const_iterator it = map.left.find(key);
        if (it != map.left.end()) {
            return it->second;
        }
        else {
            throw std::invalid_argument("Unknown type.");
        }
    }

    template<class T> T Enum_converter::string_to_key( 
            const boost::bimap<T, std::string>& map, const std::string& key)
    {
        typedef boost::bimap<T, std::string> bimap_t;
        typename bimap_t::right_const_iterator it = map.right.find(key);
        if (it != map.right.end()) {
            return it->second;
        }
        else {
            std::string s(key);
            s.append(": Unknown ");
            s.append(this->type2str<T>());
            throw std::invalid_argument(s);
        }
    }


    // 
    // "Key-to-string" specializations
    template<> std::string Enum_converter::key_to_string<datafile_format>(
            datafile_format key) {
        return key_to_string(file_formats_, key); 
    }
    template<> std::string Enum_converter::key_to_string<Marker_type>(
            Marker_type key) {
        return key_to_string(markers_, key); 
    }
    template<> std::string Enum_converter::key_to_string<Test_type>(
            Test_type key) {
        return key_to_string(tests_, key); 
    }
    template<> std::string Enum_converter::key_to_string<Var_estimate>(
            Var_estimate key) {
        return key_to_string(var_estimates_, key); 
    }
    template<> std::string Enum_converter::key_to_string<Verbosity>(
            Verbosity key) {
        return key_to_string(verbs_, key); 
    }

    // 
    // "String-to-key" specializations
    template<> datafile_format Enum_converter::string_to_key<datafile_format>(const std::string& s) {
        return string_to_key(file_formats_, s);
    }
    template<> Marker_type Enum_converter::string_to_key<Marker_type>(const std::string& s) {
        return string_to_key(markers_, s);
    }
    template<> Test_type Enum_converter::string_to_key<Test_type>(const std::string& s) {
        return string_to_key(tests_, s);
    }
    template<> Var_estimate Enum_converter::string_to_key<Var_estimate>(const std::string& s) {
        return string_to_key(var_estimates_, s);
    }
    template<> Verbosity Enum_converter::string_to_key<Verbosity>(const std::string& s) {
        return string_to_key(verbs_, s);
    }


    //
    // "Type-to-string" specializations
    template<>
        std::string Enum_converter::type2str<datafile_format>() { return "data file format"; }
    template<>
        std::string Enum_converter::type2str<Marker_type>() { return "marker type"; }
    template<>
        std::string Enum_converter::type2str<Test_type>() { return "statistical test"; }
    template<>
        std::string Enum_converter::type2str<Var_estimate>() { return "variance estimation"; }
    template<>
        std::string Enum_converter::type2str<Verbosity>() { return "verbosity level"; }
} //namespace detail
} //namespace Permory

#endif

