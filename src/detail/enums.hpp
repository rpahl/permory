// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_detail_enums_hpp
#define permory_detail_enums_hpp

namespace Permory { namespace detail {

    // Supported data formats
    enum Datafile_format{
        unknown=0, 
        compact, //permory_meta, 
        slide, 
        presto, 
        plink_tfam, plink_tped
    };
    const std::string datafile_format_to_string(const Datafile_format df) 
    {
        switch(df) {
            case compact: 
                return "compact (*.comp)"; 
                break;
                /*
            case permory_meta: 
                return "PERMORY (*.met)"; 
                break;
                */
            case slide: 
                return "SLIDE (*.slide)"; 
                break;
            case presto: 
                return "PRESTO (*.bgl)"; 
                break;
            case plink_tfam: 
                return "PLINK (*.tfam)"; 
                break;
            case plink_tped: 
                return "PLINK (*.tped)"; 
                break;
            case unknown: 
                return "unknown"; 
                break;
            default:
                return "undefined";
        }
    }

    enum Marker_type {allelic, genotype};
    const std::string marker_type_to_string(const Marker_type mt) 
    {
        switch(mt) {
            case allelic: 
                return "allelic"; 
                break;
            case genotype: 
                return "genotype"; 
                break;
            default:
                return "undefined";
        }
    }

    // Supported statistical tests
    enum Test_type { trend, trend_extended, chisq };
    const std::string test_type_to_string(const Test_type tt) 
    {
        switch(tt) {
            case trend: 
                return "T_trend"; 
                break;
            case trend_extended: 
                return "T_trend_ext"; 
                break;
            case chisq: 
                return "T_chi^2"; 
                break;
            default:
                return "undefined";
        }
    }

    // Variance estimation
    enum Var_estimate {pooled=0, controls=1, separately=2};
    const std::string var_estimate_to_string(const Var_estimate ve) 
    {
        switch(ve) {
            case pooled: 
                return "pooled"; 
                break;
            case controls: 
                return "only controls"; 
                break;
            case separately: 
                return "separately cases and controls"; 
                break;
            default:
                return "undefined";
        }
    }

    // Verbosity level
    enum Verbosity {all=0, verbose=1, normal=2, muted=3};
    const std::string verbosity_to_string(const Verbosity verb) 
    {
        switch(verb) {
            case all: 
                return "all"; 
                break;
            case verbose: 
                return "verbose"; 
                break;
            case normal: 
                return "normal"; 
                break;
            case muted: 
                return "muted"; 
                break;
            default:
                return "undefined";
        }
    }


} //namespace detail
} //namespace Permory

#endif

