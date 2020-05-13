###################################
## Nimpress Preprocessing README ##
###################################

Built and tested on R version 4.0.0 (2020-04-24)

###########
## Setup ##
###########

Run script: Nimpress_preprocess_setup.sh in nimpress directory to install R packages and download blacklisted bed files

###################################
## Running Nimpress_preprocess.R ##
###################################

To generate the required input format for NIMPRESS a preprocessing script can be used. This script, Nimpress_preprocess.R, requires the following minimum input parameters:

--file <input file in format described below>
Either --GRCh37 or --GRCh38

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

The GWAS summary statistic file requires the following header and columns:

rsID,Risk_allele,Freq,<OR or BETA>,P-value

<rsID> = the dbSNP Reference SNP cluster ID
<Risk_allele> = the allele which has been determined to be associated with the trait/ disease investigated
<Freq> = popultation frequency of the risk allele
<OR> or <BETA> = Odds ratio or beta respectively, numbers will be treated according to header, all value in the header need to be consistent with the type specified (can't mix OR and BETA)
<P-value> = significance value obtained from the GWAS study for that rsID

The underlying preprocessing script has the following functionality

- convert OR to Beta via log transformation
- query dbSNP to extract reference allele and genomic location. rsIDs which don't represent SNVs will be treated as unusable. If LDproxy is enabled, alternative rsIDs will be identified
- if remove "black-listed genomic regions" flag is set or bed file is provided, if LDpoxy parameters are not set rsIDs that fall in those regions will be removed or if LDproxy parameters are set substituted with alternative rsIDs that have an LD value > 0.9 
- check for strand flipping and define correct alleles if strand flipping has occured

Citatione for default black-listed genomic regions:
Krusche, P., Trigg, L., Boutros, P.C. et al. 
Best practices for benchmarking germline small-variant calls in human genomes. 
Nat Biotechnol 37, 555–560 (2019). https://doi.org/10.1038/s41587-019-0054-x

LDproxy_pop: a 1000 Genomes Project population: 
Select one from: https://www.internationalgenome.org/faq/which-populations-are-part-your-study


 
