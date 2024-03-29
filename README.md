# Limited interspecific gene flow in the evolutionary history of the icefish _Chionodraco_ spp.

This repository contains scripts and accessory files used for the paper "Limited interspecific gene flow in the evolutionary history of the icefish _Chionodraco_ spp."

- [command_lines_bash.sh](command_lines_bash.sh) contains the all command lines to be used in a Linux terminal.
- [command_lines_R.R](command_lines_R.R) contains the all command lines to be used in R.
- [popmap.txt](popmap.txt) is the popmap used in Stacks. This file contains all the samples included in the RAD library. Some individuals were excluded from the final analysis because they are technical replicates or because the have too many mmissing data. These samples are preceded by a hashtag (Stacks does not read lines starting with a hashtag in the popmap file).

- [TreeMix](TreeMix) contains the vcf file (dataset2) and accessory files for the analyses at species and population level with TreeMix, and the modified plotting functions.
- [Dsuite](Dsuite) contains the accessory files for the analyses at species and population level with Dsuite. The vcf file (dataset1) is too large to upload here, but can be provided on request.
- [SNAPP](SNAPP)  contains the vcf file (dataset3) and the accessory files to run SNAPP.
- [fastsimcoal](fastsimcoal) contains the vcf file (dataset4) and the files specifying the models that were simulated (no-flow: no gene flow allowed; continous: continous gene flow; on-off: gene flow allowed only during interglacial times).


Raw sequence reads analysed in this study are available under the NCBI BioProject Accession number: PRJNA1041981 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1041981/).
