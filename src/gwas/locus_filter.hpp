// Copyright (c) 2010 Roman Pahl
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

#ifndef permory_locus_filter_hpp
#define permory_locus_filter_hpp

#include "detail/config.hpp"
#include "locus.hpp"

namespace Permory { namespace gwas {

    struct Locus_filter : boost::noncopyable {
        public:
            bool operator()(const Locus& loc) {
                return do_operator(loc);
            }
        private:
            virtual bool do_operator(const Locus& loc) = 0;
    };

    struct Maf_filter : public Locus_filter {
        public:
            Maf_filter(const std::string& s, double min=0.0, double max=0.5)
                : type_(s), min_maf(min), max_maf(max)
            {}
        private:
            bool do_operator(const Locus& loc) {
                double maf = loc.maf(type_);
                return (maf > min_maf && maf < max_maf);
            }
            const std::string type_;
            double min_maf;
            double max_maf;
    };

    struct Polymorph_filter : public Locus_filter {
        private:
            bool do_operator(const Locus& loc) {
                return loc.isPolymorph();
            }
    };

} // namespace gwas
} // namespace Permory

#endif // include guard
