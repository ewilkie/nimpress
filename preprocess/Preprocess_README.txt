###################################
## Running Nimpress_preprocess.R ##
###################################

Located in the same dir are the following dependencies:
- Nimpress_preprocess_functions.R script containing dependent functions
- Suppl data is located and downloaded to the Suppl folder 
- Output dir Nimpress_preprocess_Output 

###########
## Input ##
###########

Minimum input is a file to be processed in the following format containing header and columns:

rsID,Risk_allele,Freq,<OR or BETA>,P-value

<rsID> = the dbSNP Reference SNP cluster ID
<Risk_allele> = the allele which has been determined to be associated with the trait/ disease investigated
<Freq> = popultation frequency of the risk allele
<OR> or <BETA> = Odds ratio or beta respectively, numbers will be treated according to header, all value in the header need to be consistent with the type specified (can't mix OR and BETA)

############
## Output ##
############

Finalised data to be used as input to Nimpress is located in the output folder (Nimpress_preprocess_Output) with suffix NIMPRESS_input.txt

Intermediate file that can be used for debugging has suffix Intermediate_results.csv. This file contains the following information:
- dbSNP location, REF and ALT alleles. 
- the GRCh37 allele based on the above location from which strand is inferred
- FLAG.RM indicated whether that rsID has been excluded from the output
- FLAG.AMBIGUOUS indicates whether the REF and ALT alleles are complementary
- BED.COVERAGE whether the rsID is located in the blacklisted region
- INPUT.RISK.TYPE whether risk allele is the REF or the ALT allele
- FLAG.RISK.FLIPPED whether the risk allele is complementary to the REF or ALT and therefore strand flipping has occured
- FLAG.LDPROXY whether a substituted SNV will be used in the output
- LDPROXY.rsID	rsID of substituted SNV
- LDPROXY.CHR	chromosome of substituted SNV
- LDPROXY.START	location of substituted SNV
- LDPROXY.REF.ALLELE	of substituted SNV
- LDPROXY.ALT.ALLELE	of substituted SNV
- INPUT.ALLELE.FREQ	input frequency
- INPUT.BETA input beta (converted from OR if that is supplied instead)

###################
## Functionality ##
###################

-Only works on genome build version GRCh37

Notes on Linkage disequalibrium proxy substitution:
SNPs can only be substituted if they are in the 1000 Genome reference panel and have GRCh37 coordinates. Therefore if GRCh38 is selected and LDproxy parameters specified, the LDproxy functionality will be turned off

monoallelic in the population selected will be removed

dbSNP version 151 (Database of Single Nucleotide Polymorphisms [DBSNP], 2007) is used to match query RS numbers with the genomic coordinates (GRCh37) of the SNPs of interest

LDproxy_pop: a 1000 Genomes Project population: 
Select one from: https://www.internationalgenome.org/faq/which-populations-are-part-your-study

The underlying preprocessing script has the following functionality

- convert OR to Beta via log transformation
- query dbSNP to extract reference allele and genomic location. rsIDs which don't represent SNVs will be treated as unusable
- if remove "black-listed genomic regions" flag is set or bed file is provided, but if LDpoxy parameters are not set rsIDs that fall in those regions will be removed. If LDproxy parameters are set, rsIDS without coverage will be substituted with alternative rsIDs that have an LD value > 0.9

- check for strand flipping and define correct alleles if strand flipping has occured

## NCBI information on strand orientation
## ncbi.nlm.nih.gov/books/NBK44476/#Reports.how_do_i_determine_orientation_o


###########
## Usage ##
###########

Minimal command to run preprocessing: Rscript Nimpress_preprocess.R --file=Example/Example_File_to_process.csv 

Optional parameters:
--offset:   DEFAULT: 0.0

If LDproxy is desired both of the following parameters are required:
--LDproxy_pop: a 1000 Genomes Project population more details below: DEFAULT: GBR
--LDproxy_token <generate token via https://ldlink.nci.nih.gov/?tab=apiaccess>

If removal of rsID in difficult to sequence genomic regions is desired or substitueted with rsIDs in LD if above options are set use:
--blacklisted_bed: bed file containng regions. Will be downloaded and places in /Suppl folder if Nimpress_preprocess_setup.sh is run

The --file input needs to be a .csv file, where each line represents a file to be preprocessed in the following format:

<GWAS summary statistic file and path in csv format>,<description>,<citation>

The first two parameters are required, description and citiation can be left blank however are encouraged to be used for ease of recordkeeping. The name of the input file will be used as the basis for the name of the output file and therefore should be used to distinguish different subtypes if desired. 

#########################
## Built and tested on ##
#########################

R version 4.0.0 (2020-04-24)

Libraries:
pacman_0.5.1 
docopt_0.7.1
data.table_1.12.8
GenomicRanges_1.40.0
rentrez_1.2.2
BSgenome.Hsapiens.UCSC.hg19_1.4.3

###############
## Citations ##
###############

Default black-listed genomic regions:
Krusche, P., Trigg, L., Boutros, P.C. et al. 
Best practices for benchmarking germline small-variant calls in human genomes. 
Nat Biotechnol 37, 555–560 (2019). https://doi.org/10.1038/s41587-019-0054-x

LDlinkR:
Myers Timothy A., Chanock Stephen J., Machiela Mitchell J.
LDlinkR: An R Package for Rapidly Calculating Linkage Disequilibrium Statistics in Diverse Populations  
Frontiers in Genetics 11, 157 (2020). https://doi.org/10.3389/fgene.2020.00157   

