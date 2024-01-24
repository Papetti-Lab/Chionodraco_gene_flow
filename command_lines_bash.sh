### ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ###
### This file contains the command lines used in a linux terminal to run the analysed presented in the paper "Limited interspecific gene flow in the evolutionary history of the icefish genus Chionodraco spp." ###
### ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ###


	### Quality control with FastQC ###

IN=/path/to/raw/reads/folder
OUT=/path/to/cleaned/reads/folder
# Define an array with fastq file names
FILE=($IN/*.fastq.gz)
echo ${FILE[@]} # check the file names
# Run FastQC
for file in ${FILE[@]}
do
    fastqc -t 2 ${file} -o $OUT
done





	### Demultiplex samples, discard low quality reads and trim them with Stacks ##

# Define input and output folders
IN=/path/to/raw/reads/folder
OUT=/path/to/output/the/processed/files/folder
# Command line for the first library
process_radtags -1 $IN/2791_S31_L002_R1_001.fastq.gz -2 $IN/2791_S31_L002_R2_001.fastq.gz -o $OUT -b $IN/barcodes_library_01.txt -e sbfI -r -c -q -t 143
# Command line for the second library
process_radtags -1 $IN/sample_S0_L006_R1_001.fastq.gz -2 $IN/sample_S0_L006_R2_001.fastq.gz -o $OUT -b $IN/barcodes_library_02.txt -e sbfI -r -c -q -t 143





	### Remove PCR clones with Stacks ###

# Define input and output folders
IN=/path/to/cleaned/reads/generated/by/process_radtags/folder
OUT=/path/to/output/the/processed/files/folder

# Define an array with first input file names in a set of paired-end sequences
PAIR_1=($IN/*.1.fq.gz)
echo ${PAIR_1[@]} # check the file names
# Define an array with second input file names in a set of paired-end sequences
PAIR_2=($IN/*.2.fq.gz)
echo ${PAIR_2[@]} # check the file names

# Clone filter
for i in "${!PAIR_1[@]}" ; do
	clone_filter -P -1 ${PAIR_1[i]} -2 ${PAIR_2[i]} -i gzfastq -o $OUT > clone-filter-log.txt
 done





	### Align to reference genome with BWA ###

GENOME=/path/to/your/genome/file
INDEX=/path/to/your/index/bwa_index # the path to the folder where you want the index to be placed followed by the name of the index that will be produced ("bwa_index" in this case)
# Generate the index
bwa index $GENOME -p $INDEX > bwa_index.log

# Define path to input and output folders
IN=/path/to/reads/to/be/aligned/folder
OUT=/path/to/output/the/processed/files/folder

# Move into the input folder
cd $IN

# Define an array with all file names without extension
NAMES+=($(echo "$(ls *.fq.gz)" | cut -f 1 -d '.' | sort -u))
echo ${NAMES[@]} # check the file names

# Align
for i in "${!NAMES[@]}" ; do
	bwa mem -t 8 -M $INDEX $IN/${NAMES[i]}.1.1.fq.gz $IN/${NAMES[i]}.2.2.fq.gz | samtools view -b -h | samtools sort --threads 8 > $OUT/${NAMES[i]}.bam
done





	### Calculate alignment statistics ###
IN=/path/to/bam/files/folder
OUT=/path/to/output/the/processed/files/folder
cd $IN
# Define an array with all file names without extension
NAMES+=($(echo "$(ls *.bam)" | cut -f 1 -d '.'))
echo ${NAMES[@]} # check the file names
# Run samtools flagstat
for i in "${!NAMES[@]}" ; do
    samtools flagstat $IN/${NAMES[i]}.bam -O tsv --threads 8 >> $OUT/${NAMES[i]}_flagstat.txt
done





	### SNP calling with Stacks ###

# gstacks
IN=/path/to/bam/files/folder
OUT=/path/to/output/the/processed/files/folder
POPMAP=/path/to/popmap/file
gstacks -I $IN -M $POPMAP -O $OUT -t 8

# populations (note that the arguments -e and --merge-sites are no lognger available at the moment this document was written - October 2023, Stacks version 2.65)
IN=/path/to/gstacks/folder
OUT=/path/to/output/the/processed/files/folder
POPMAP=/path/to/popmap/file
populations -P $IN -M $POPMAP -O $OUT -p 4 -r 0.75 --max-obs-het 0.75 -e sbfI --merge-sites --ordered-export --vcf  -t 8

# populations was first run with all samples and technical replicates.
# Then, consistency of SNP calling between replicates was checked with the script SNPs_error.R, available at https://github.com/AliciaMstt/RAD-error-rates/blob/master/SNPs_error.R.
# One sample per replicate pair and individuals with more than 40% of missing data were removed and the function populations was run again. 





	### Filtering ###

cd /path/to/folder/where/you/want/the/output

# Command line to produce dataset 1 used for fineRADstructure and Dsuite
INPUT_0=/path/to/populations.snps.vcf	# the VCF file produced by stacks
vcftools --vcf $INPUT_0 --remove-indels --minGQ 30 --min-meanDP 10 --max-meanDP 40 --minDP 10 --maxDP 40 --mac 2 --max-missing 0.75 --recode --stdout --out dataset1 > dataset1.vcf

# Command line to produce dataset 2 used to compute Fst and PCoA, and run Admixture and TreeMix
# Convert from vcf to plink
./plink --vcf dataset1.vcf --const-fid --allow-extra-chr --make-bed --out dataset1
#Prune for LD
./plink --bfile dataset1 --const-fid --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out dataset1-pruned
# Remove from dataset1 linked SNPs identified by plink
vcftools --vcf dataset1.vcf --remove-indels --exclude dataset1-pruned.prune.out --recode --stdout --out dataset1-pruned > dataset1-pruned.vcf
# Remove any SNP that is in Hardy-Weinberg disequilibrium in at least one population
./filter_hwe_by_pop.pl -v dataset1-pruned.vcf -p popmap.txt -h 0.01 -c 0 -o dataset2





	### fineRADstructure ###

IN=/path/to/dataset1/dataset1.vcf # filtered vcf file
OUT=/path/to/output/the/processed/files/folder
POPMAP=/path/to/popmap/file
# Convert VCF file to radpainter format
populations -V $IN -O $OUT -M $POPMAP --radpainter
# Retain no more than ten SNPs per RAD locus
cd $OUT
python Stacks2fineRAD.py -i dataset1-haplotypes.tsv -n 10

# Calculate the co-ancestry matrix
FILE=dataset1 # name of the input file, not the input file itself
RADpainter paint $FILE.radpainter
# Assign individuals to populations
finestructure -x 100000 -y 100000 -z 1000 ${FILE}_chunks.out ${FILE}_chunks.mcmc.xml
# Tree building
finestructure -m T -x 100000 ${FILE}_chunks.out ${FILE}_chunks.mcmc.xml ${FILE}_chunks.mcmcTree.xml





		### Dsuite ###
cd /path/to/folder/where/you/want/the/output		
FILE=/path/to/dataset1.vcf # filtered vcf file
TREE=/path/to/tree/file # tree in Newick format
SETS=/path/to/sets/file # species/population map indicating the outgroup
# Run Dsuite
Dsuite Dtrios -t $TREE $FILE $SETS





	### Admixture ###

# Convert from vcf to plink
cd /path/to/folder/where/you/want/the/output	
IN=/path/to/dataset2/file
./plink --vcf $IN --const-fid --allow-extra-chr --make-bed --out dataset2
# Prepare the input for Admixture
FILE=dataset2 # name of the file, not the file itself
awk '{$1=0;print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim

# Run Admixture
for K in {2..11} # adjust the number of interation of the loop according to the number of clusters you want to test
do
	admixture32 --cv $FILE.bed $K | tee log${K}.out
done

# Visualize cross-validation (CV) errors
grep -h CV log*.out
# Save cross-validation errors in a convenient format to plot them in R
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' | sort -k 1 -n > $FILE-cv_errors.txt





	### TreeMix ###

# Conver VCF file to TreeMix format
IN=/path/to/dataset2/file
OUT=/path/to/output/the/processed/files/folder
POPMAP=/path/to/popmap/file
populations -V $IN -O $OUT -M $POPMAP --treemix
cd $OUT
echo "$(tail -n +2 dataset2.p.treemix)" > dataset2
bgzip dataset2
# Run TreeMix
FILE=dataset2 # name of the input file, not the input file itself
for i in {0..5} # adjust the number of interation of the loop according to the number of migration edges you want to test
do
    treemix -i $FILE.gz -m $i -o $FILE.$i -root C.wilsoni -bootstrap -k 1000 > treemix_${i}_log &
done





	### RagTag ###

# RagTag of Aethotaxis genome with D. mawsoni as reference
REF=/path/to/reference/genome/file
QUERY=/path/to/query/genome/file
UNIMAP=/path/to/unimap/unimap # aligner chosen for this analysis
OUT=/path/to/output/the/processed/files/folder

conda activate myenv # activate the python environment where ragtag is installed
ragtag.py scaffold $REF $QUERY -C -u -t 8 --aligner $UNIMAP -o $OUT



	### Run Dsuite for each chromosome separately ###

# Use the command lines described in the previous sections to align fastq file to the new genome generated with RagTag and perform SNP calling

# Filter vcf file
cd /path/to/folder/where/you/want/the/output
INPUT_0=/path/to/populations.snps.vcf	# the VCF file produced by stacks
vcftools --vcf $INPUT_0_ragtag --remove-indels --minGQ 30 --min-meanDP 10 --max-meanDP 40 --minDP 10 --maxDP 40 --mac 2 --max-missing 0.75 --recode --stdout --out dataset1-ragtag > dataset1-ragtag.vcf

# Run Dsuite 
FILE=/path/to/dataset1-ragtag.vcf # vcf file
TREE=/path/to/tree/file # tree in Newick format
SETS=/path/to/sets/file # species map indicating the outgroup

for i in {0,1,3,4,5,6,8,10,11,13,21,22,23,24,25,26,28,30,31,32,33,34,35,38,41} # these numbers are derived from how the chromosomes were numbered in the reference genome used in ragtag
do
  mkdir chr$i
  vcftools --vcf $FILE --chr HiC_scaffold_${i}_RagTag --recode --stdout --out rag$i > rag$i.vcf
  mv rag$i.vcf chr$i"/"
  mv rag$i.log chr$i"/"
  cp SETS.txt chr$i"/"
  cp tree-species.nwk chr$i"/"
  cd chr$i
  Dsuite Dtrios -t $TREE rag$i.vcf $SETS
  cd ../
done





	### SNAPP ###

# Filter the vcf file produced by Stacks
cd /path/to/folder/where/you/want/the/output
INPUT_0=/path/to/populations.snps.vcf	# the VCF file produced by stacks
INDV=/path/to/text/file/containg ten random individuals per species
vcftools --vcf $INPUT_0 --remove-indels --minGQ 40 --min-meanDP 10 --max-meanDP 40 --minDP 10 --maxDP 40 --max-missing 0.15 --thin 10000 --keep $INDV --recode --stdout --out snapp0 > snapp0.vcf
vcftools --vcf snapp0.vcf --keep $INDV --recode --stdout > snapp.vcf

# Generate the SNAPP input file from the filtered vcf file
INPUT_3=/path/to/snapp.vcf # filtered vcf file 
STARTING_TREE=/path/to/starting.tre # file with a starting tree, that is starting.tre
SPECIES=/path/to/species.txt # file with a table of samples and corresponding species IDs
CONSTRAINT=/path/to/constraint.txt # file specifying the age constraint 

ruby snapp_prep.rb -v $INPUT_3 -t $SPECIES -c $CONSTRAINT -x dataset_3.xml -o dataset_3 -s $STARTING_TREE

# Once the xml file is produced, you have to decide whether you want to use SNAPP with a molecular clock model, to obtain a time-calibrated tree but imposing a constant effective population size (Ne) along the branches, or allowing for different Ne per branch (combining the two is possible but it is even more computational demanding).
# If you want to allow different population sizes for the different branches of the tree, you’ll have to open the XML file and uncomment the GammaMover and RateMixer operators.
# To uncomment some text, delete the symbols <!-- and --> surrounding the text. Example: <!-- commented text -->
# These are the two lines to comment/uncomment depndig on the type on analysis you want to perform

#	<!--
#	<operator id="GammaMover" spec="snap.operators.GammaMover" coalescenceRate="@coalescenceRate" scale="0.5" weight="8.0"/>
#   <operator id="RateMixer" spec="snap.operators.RateMixer" coalescenceRate="@coalescenceRate" scaleFactors="0.25" tree="@tree" weight="1.0"/>
#   -->


# If you want to perform both kind of analysis, you need two separate xml file. Duplicate the xml file, comment/uncomment the GammaMover and RateMixer parts in the two copies and give two different names to the two copies.
# Then, you need two change two fields in the Logger part of the xml file in order to match the names of the output files with the name of the input files.
# Change the file name in the following two lines in the "Loggers" section of the xml file.

# < logger fileName="filename.log" logEvery="250">
# < logger fileName="filename.trees" logEvery="250" mode="tree">


# Run SNAPP
INPUT_SNAPP=/path/to/dataset_3.xml # the xml file generated by the ruby script 
beast -threads 20 -overwrite $INPUT_SNAPP

# If you are interested in population size, note that the output file does not contain the estimated value of population size but Theta estimates.
# To convert Theta estimates into Ne, use the equation Ne = Theta ÷ (4 × r ÷ ng), where r is the substitution-rate estimate (= the rate of the strict clock estimated by SNAPP) and ng is the number of generations per time unit
# More details are given in Stange et al. 2018, https://doi.org/10.1093/sysbio/syy006.





	### fastsimcoal ###
# Filter the vcf file produced by Stacks
cd /path/to/folder/where/you/want/the/output
INPUT_0=/path/to/populations.snps.vcf	# the VCF file produced by stacks
vcftools --vcf $INPUT_0 --remove-indels --minGQ 40 --min-meanDP 10 --max-meanDP 40 --minDP 10 --maxDP 40 --max-missing 0.85 --thin 1000 --recode --stdout --out fastsimcoal > fastsimcoal.vcf
# Convert VCF to SFS
# Open the file VCFtoSFS_downsample.sh with a text editor and modify the settings to define the VCF file to convert, the tag for the resulting files and the minimum sample size per population
# The script VCFtoSFS_downsample.sh is availbale at http://cmpg.unibe.ch/software/fastsimcoal27/additionalScripts.html, in the folder buildSFSFromVCF.zip
VCFtoSFS_downsample.sh

# Run fastsimcoal
# Gene flow estimation for model 2  and model 3 (100 replicates)
PREFIX=tag used to identify the different models to be tested

for i in {1..100}
 do
   mkdir run$i
   cp ${PREFIX}.tpl run$i"/"
   cp ${PREFIX}.est run$i"/"
   cp${PREFIX}_jointMAFpop1_0.obs run$i"/"
   cp ${PREFIX}_jointMAFpop2_0.obs run$i"/"
   cp ${PREFIX}_jointMAFpop2_1.obs run$i"/"
   cd run$i
   fsc26 -t ${PREFIX}.tpl -e ${PREFIX}.est -m -0 -C 10 -M -n 100000 -L 50 -c 10 -B 10 -q
   cd ..
done

# Extract the files of the best run among the 100 replicates and copies them into a new folder (called bestrun) with the script available at https://raw.githubusercontent.com/speciationgenomics/scripts/master/fsc-selectbestrun.sh
fsc-selectbestrun.sh

# Simulation of the three scenarios (100 replicates per scenario)
PREFIX=name used to identify the input for the different models to be simulated

mkdir par_runs
cd par_runs

for i in {1..100}
 do
   mkdir run$i
   cp ../${PREFIX}_maxL.par run$i
   cp ../${PREFIX}_jointMAFpop1_0.obs run$i/${PREFIX}_maxL_jointMAFpop1_0.obs
   cp ../${PREFIX}_jointMAFpop2_0.obs run$i/${PREFIX}_maxL_jointMAFpop2_0.obs
   cp ../${PREFIX}_jointMAFpop2_1.obs run$i/${PREFIX}_maxL_jointMAFpop2_1.obs
   cd run$i
 fsc26 -i ${PREFIX}_maxL.par -n 1000000 -m -0 -c 10 -B 10 -q
# Fastsimcoal will generate a new folder called ${model} and write files in there
   cd ..
 
# collect the lhood values
sed -n '2,3p' run$i/${PREFIX}_maxL/${PREFIX}_maxL.lhoods  >> ${PREFIX}.lhoods

done

cat ${PREFIX}.lhoods | awk '{print NR,$1}' | sort -k 2 >> sorted.txt


# Make a boxplot for the likelihoods present in the file sorted.txt for all the different models
# Plot observed and expected site frequency spectra in R


