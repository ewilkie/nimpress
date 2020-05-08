##########################
## Preprocessing README ##
##########################

To generate the required input format for NIMPRESS a preprocessing script can be used. This script, Nimpress_preprocess.R, requires the following input parameters:

--files_to_process <input file in format described below>
--genome_version <GRCh37 or GRCh38>

Optional parameters:
--offset <DEFAULT 0.0>
--LDproxy_pop <background population DEFAULT: GBR>
--LDproxy_token <generate token via https://ldlink.nci.nih.gov/?tab=apiaccess>
--blacklisted_bed <bed file containing blacklisted regions>

The files_to_process file needs to be a .csv file, where each line represents a file to be preprocessed in the following format:

<GWAS summary statistic file and path in csv format>,<description>,<citation>

The first two parameters are required, while description and citiation can be left blank however are encouraged to be used for ease of recordkeeping. The name of the input file will be used as the basis for the name of the output file and therefore should be used to distinguish different subtypes if desired. 

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
- if a "black-listed genomic regions" bed file is provided, rsIDs that fall in those regions will be substituted with alternative rsIDs that have an LD value > 0.9
- check for strand flipping and define correct alleles if flipping has occured

 
