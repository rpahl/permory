

/**
 * @author Roman Pahl
 */

#ifndef permory_locus_hpp
#define permory_locus_hpp

#include <string>
#include <vector>

#include "config.hpp"

namespace Permory 
{
    class Individual; // forward declaration

    class Locus {
        public: 
            enum Chr {
                none=0, chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, 
                chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, 
                chr17, chr18, chr19, chr20, chr21, chr22, X, Y};

            // Members
            size_t id;     //unique id
            string name;     //rs-id or other Locus identifier
            Chr chr;            //chr1-chr22, X, Y or na if undefined
            size_t bp;          //base pair position in bp units
            double cm;          // cM map position
            std::vector<string> allele;
            std::vector<Permory::Individual>* individuals; 

            // Ctor
            Locus(
                    size_t id,
                    string name="",
                    Chr chr=none
                 ) 
                : id(id), name(name) 
            {
                //this->individuals = 0;
            }

            // Inspector 
            bool operator<(size_t i) const { return id < i; }
            bool operator<(double d) const { return tsMax_ < d; }
            bool operator<(const Locus& x) const { 
                return (chr < x.chr || (chr == x.chr && bp < x.bp) );
            }
            bool operator==(const Locus& x) const { return ( id == x.id ); }

            // Modifier
            template<int K, int L> void add_test_stat(
                    std::vector<double>::const_iterator start,
                    std::vector<double>::const_iterator end);
        private:
            std::vector<double> ts_;    //test statistics
            double tsMax_;              //max(ts_)
    };

    // ========================================================================
    // Locus implementations
    template<int K, int L> inline void Locus::add_test_stat(
            std::vector<double>::const_iterator start,
            std::vector<double>::const_iterator end)
    {
        copy (start, end, back_inserter(ts_));
        tsMax_ = std::max(tsMax_, *std::max_element(ts_.begin(), ts_.end()));
    }


} // namespace Permory

#endif // include guard
