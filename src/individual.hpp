

/**
 * @author Roman Pahl
 */

#ifndef permory_individual_hpp
#define permory_individual_hpp

#include <string>
#include <vector>

#include "config.hpp"

namespace Permory 
{
    class Locus; // forward declaration

    // Information referring to the individual that is not a measurement 
    class Individual { 
        public:
            enum Sex {nosex=0, male, female};

            // Members
            size_t id;     //individual ID
            string name;
            Sex sex;       
            string famid;     //family ID:
            string patid;     //paternal ID
            string matid;     //maternal ID
            std::vector<Locus>* loci;

            // Ctor
            Individual(size_t id, Sex sex=nosex) : id(id), sex(sex) 
            {
                this->loci = 0;
            }

            // Operator for use with containers
            bool operator<(const Individual& i) const { return id < i.id; }
    };

    // Data measured for an individual at a particular experiment/study
    struct Measurement {
        bool affected;
        double phenotype;

        Individual* pi;
    };

    template<class C> struct Data_per_individual {
        std::vector<C> dat;
        Individual* pi;
    };
} // namespace Permory

#endif // include guard


