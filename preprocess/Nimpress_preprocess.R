###############################
## TO DO - might be optional ##
###############################

## with and without out master file, fo running single and multiple files
## does Freq need to be < 1?? 
## verbose turn off and on
## some statisitcs

##########################################
## GWAS summary data curations pipeline ##
##########################################

## turn warnings off
oldw <- getOption("warn")
options(warn = -1)

## load command option library
pacman::p_load(docopt)

#############
## docopts ##
#############

## example command:  Rscript Nimpress_preprocess.R --file /Users/ewilkie/Documents/Polygenic/nimpress/preprocess/Example/Example_File_to_process.txt --GRCh37 --blacklisted_bed /Users/ewilkie/Documents/CCI_general_data_files/GRCh37_alldifficultregions.bed --LDproxy_pop GRB --LDproxy_token cbe1b45bc8be


'NIMPRESS preprocess
Usage:
  Nimpress_preprocess.R --file=<file_to_process> (--GRCh37 | --GRCh38) [(--LDproxy_pop=<BG population> --LDproxy_token=<token>) --blacklisted_bed=<bed> --outpath=<outpath> --offset=<offset>]  
  Nimpress_preprocess.R --file=<file_to_process> --GRCh37 --blacklisted_bed=<bed>
  Nimpress_preprocess.R --file=<file_to_process> --GRCh37 --blacklisted_bed=<bed> --LDproxy_pop=<BG population> --LDproxy_token=<token>
  Nimpress_preprocess.R --file=<file_to_process> --GRCh37 --LDproxy_pop=<BG population> --LDproxy_token=<token>
  Nimpress_preprocess.R (-h | --help)
  Nimpress_preprocess.R --version
Arguments:
    --file=<file_to_process>     See Example folder for completed template
                          Text file containing risk loci file in template format:
                          <GWAS summary statistic file and path in csv format>,<description>,<citation>
                          GWAS summary statistic file contents:
                          rsID,Risk_allele,Freq,<OR or BETA>,P-value
                          <rsID> = the dbSNP Reference SNP cluster ID
                          <Risk_allele> = the allele associated with the trait or disease investigated
                          <Freq> = popultation frequency of the risk allele
                          <OR> or <BETA> = Odds ratio or beta respectively, all value need to be of same type
                          <P-value> = significance value obtained from the GWAS study for that rsID

    --GRCh37              use genome version GRCh37
    --GRCh38              use genome version GRCh38
     
Options:
  -h --help                         Show this screen.
  --version                         Show version.
  --outpath=<outpath>               Path to output location [DEFAULT: ./Nimpress_preprocess_Output]
  --offset=<offset>                 Offset for NIMPRESS [DEFAULT: 0.0]
  --blacklisted_bed=<bed>           Blacklisted bed region
  --LDproxy_pop=<BG population>     Background populations for LDproxy
  --LDproxy_token=<token>           Generate token via https://ldlink.nci.nih.gov/?tab=apiaccess

    
' -> doc

arguments <- docopt(doc, version = 'NIMPRESS Preprocess for R version 4.0.0 (2020-04-24)\n')
print(arguments)
stop("just checking arguments")

## for testing only 
setwd("/Users/ewilkie/Documents/Polygenic/nimpress/preprocess/")
arguments <- list()
arguments$GRCh37 = TRUE
arguments$blacklisted_bed = "/Users/ewilkie/Documents/CCI_general_data_files/GRCh37_alldifficultregions.bed"
arguments$file = "./Example/Example_File_to_process.csv"
arguments$LDproxy_pop="GRB"
arguments$LDproxy_pop="cbe1b45bc8be"

##########################################
## testing that still needs to be done: ##
##########################################


## GRCh38 - but need different file for that

###################
## Initial setup ##
###################

###############
## libraries ##
###############

message("[1/..] Loading libraries...")
pacman::p_load(data.table,GenomicRanges,rentrez)
## rentrez,bedr,LDlinkR,stringr,GenomicRanges

###################
## Genomic files ##
###################

message("[2/..] Setting up genomic environment...")

## get the correct assembly file
if(arguments$GRCh37 == TRUE){
  gv <- "GRCh37"
  assembly_file <- "./Suppl/GCF_000001405.13_GRCh37_assembly_report.txt"
}else if (arguments$GRCh38 == TRUE){
  gv <- "GRCh38"
  assembly_file <- "./Suppl/GCF_000001405.26_GRCh38_assembly_report.txt"
}

ass <- read.table(assembly_file, header=F, sep="\t",stringsAsFactors = F)
ass_sub <- ass[,c(3,7)]
assembly <- ass_sub[grep("NC_", ass_sub[,2]),]
colnames(assembly) <- c("CHR", "NC_CHR")
#print(assembly)

## setup blacklisted genome bed file

if(!is.null(arguments$blacklisted_bed)){
  bedfile = arguments$blacklisted_bed
  ovlp <- fread(bedfile, header = FALSE, stringsAsFactors = FALSE)
  colnames(ovlp) <- c("chr", "start", "end")
  gr <- makeGRangesFromDataFrame(ovlp, keep.extra.columns = TRUE)
}else{
  bedfile = NULL
}


#################
## Input files ##
#################

message("[3/..] Reading master file...")

master_file <- read.table(arguments$file, sep=",", header=T)
master_file.list <- split(master_file, seq(nrow(master_file)))

## load processing functions files
source("./Nimpress_preprocess_functions.R")

message("[4/..] Starting file processing...")

## set up loop when single run is finished
#for(f in 1:length(master_file.list)){
#}
f <- 2

message(paste("[", 3+f ,"/..] Processing file:", master_file.list[[f]]$GWAS_summary_statistic_file_and_path, sep="" ))
input <- master_file.list[[f]]$GWAS_summary_statistic_file_and_path

## check and format input file
gf_ok <- check_gwas_file(input)
#print(gf_ok)

## extract unique rsID
rsID_ind <- grep("rsID", colnames(gf_ok))
rsID <- gf_ok[,rsID_ind]
rsIDu <- as.vector(unique(rsID))
  
##########################################
## get rdID genomic location from dbSNP ##
##########################################

## list containing "rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE" - for Risk allele check
## this loop takes about 28.85696 / 13 = 2.219766 secs per SNPS
ela <- Sys.time()
rsID_loc <- list()

for (rsid in 1: length(rsIDu)){
  message(paste("Getting info for : ", rsIDu[rsid],sep=""))
  res <- getrsID_info(rsIDu[rsid])
  rsID_loc[[rsid]] <- res
}

ela <- Sys.time() - ela
print(ela)

rsID_loc_df <- do.call(rbind, rsID_loc)
print(rsID_loc_df)

##################
## Get coverage ##
##################

## tghis obviously doesn't work with deletions
get_cov <- function(urer,s){
  snp <- GRanges(seqnames=as.numeric(urer[s,"CHR"]), ranges=IRanges(start=as.numeric(urer[s,"START"])-1, end=as.numeric(urer[s,"START"]), starts.in.df.are.0based))
  hits <- findOverlaps(gr,snp)
  if(length(hits@from) > 0){
    bedcov = FALSE
  }else{
    bedcov = TRUE
  }
  return(bedcov)
}


ela <- Sys.time()


rsID_loc_df <- do.call(rbind, rsID_loc)
colnames(rsID_loc_df) <- c("CHR","START","rsID")
urer <- rsID_loc_df



if(bedfile == "NULL"){
  bedcov <- FALSE
  urercov <- cbind(urer, cov)
}else{
  urer <- rsID_loc_df
  cov <- vector()
  for (s in 1:nrow(urer)){
    nc <- get_cov(urer,s)
    cov <- c(cov, nc)
    
  }
  urercov <- cbind(urer, cov)
}

ela <- Sys.time() - ela
print(ela)

########!!!!!!!! swapped true with false since cov file type has changed
## this is for getLDproxy function
original_SNP <- urercov[,"rsID"]






#############
## LDproxy ##
#############

## check if LDproxy is on

## check if population is one of the allowed - is population needed??
## check if proxy is correct



## if off and bedfiles is provided - remove those that fall inside the Bedfile region
## if on, check for for overlap 




