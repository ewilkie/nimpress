##########################
## Preprocessing README ##
##########################

To generate the required input format for NIMPRESS a preprocessing script can be used. This script, Nimpress_preprocess.R, requires the following minimum input parameters:

--file <input file in format described below>
Either --GRCh37 or --GRCh38

Optional parameters:
--offset:   DEFAULT: 0.0

If LDproxy is desired both of the following parameters are required:
--LDproxy_pop: background population DEFAULT: GBR>
--LDproxy_token <generate token via https://ldlink.nci.nih.gov/?tab=apiaccess>

If removal of rsID in difficult to sequence genomic regions is desired use either of the following:
-r: removal of SNPS in GIAB v2.0 stratification BED files for relevant genome version
--blacklisted_bed: private region in bed file 

The files_to_process file needs to be a .csv file, where each line represents a file to be preprocessed in the following format:

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
- query dbSNP to extract reference allele and genomic location
- if remove "black-listed genomic regions" flag is set or bed file is provided, rsIDs that fall in those regions will be removed if LDpoxy parameters are not set or substituted with alternative rsIDs that have an LD value > 0.9 if LDproxy parameters are set
- check for strand flipping and define correct alleles if strand flipping has occured

Citatione for default black-listed genomic regions:
Krusche, P., Trigg, L., Boutros, P.C. et al. 
Best practices for benchmarking germline small-variant calls in human genomes. 
Nat Biotechnol 37, 555–560 (2019). https://doi.org/10.1038/s41587-019-0054-x

 
