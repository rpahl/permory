PERMORY version 1.1 -- Nov 7, 2010
------------------------------------
    PERMORY is a free, open-source software, designed to perform efficient 
    permutation tests for large-scale genetic data sets (e.g. genome-wide
    association studies (GWAS)).


License notes
-------------
    PERMORY Copyright (c) 2010-2011 Roman Pahl, IMBE Marburg
                               2011 Volker Steiß
    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)


Package content
---------------
    - data (folder)     # example data of different formats
        tiny.bgl.gz     # marker data (Beagle format as used by PRESTO)
        tiny.tfam       # trait data (transposed fileset PLINK format)
        tiny.tped       # marker data (transposed fileset PLINK format)
        tinyA.slide     # allelic marker data (SLIDE format)
        tinyG.slide     # genotype marker data (SLIDE format)
    - example (folder)  # example scripts
        mpi_marc.sh     # example MPI submission script
    compactor           # data compactor executable
    LICENSE_1_0.txt     # license information
    permory             # permory executable
    permory.cfg         # example configuration file 
    QUICKSTART.txt      # general usage and application examples
    QUICKSTART_MPI.txt  # parallel computing using message passing interface
    README.txt          # this file


Getting started
---------------
    For a quick start and to see PERMORY in action, try the application examples
    in the QUICKSTART.txt that came with the package.


Note
----
The X and Y chromosome are not handled by PERMORY and must be removed from the
data set before starting the analysis.


---
Roman Pahl
http://www.permory.org
