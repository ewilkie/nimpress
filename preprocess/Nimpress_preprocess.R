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
  Nimpress_preprocess.R --file=<file_to_process> (--GRCh37 | --GRCh38) [(--LDproxy_pop=<BG population> --LDproxy_token=<token>) (--blacklisted_regions_file=<bed> | --remove_blacklisted_regions) --outpath=<outpath> --offset=<offset>]  
  Nimpress_preprocess.R --file=<file_to_process> --GRCh37 --remove_blacklisted_regions
  Nimpress_preprocess.R --file=<file_to_process> --GRCh37 --remove_blacklisted_regions --LDproxy_pop=<BG population> --LDproxy_token=<token>
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
  --LDproxy_pop=<BG population>     Background populations for LDproxy
  --LDproxy_token=<token>           Generate token via https://ldlink.nci.nih.gov/?tab=apiaccess
  --remove_blacklisted_regions      Difficult to sequence regions from GIAB for rsID removal/substitution
  --blacklisted_regions_file=<bed>   Provide different bedfile for rsID removal/substitution 

    
' -> doc

arguments <- docopt(doc, version = 'NIMPRESS Preprocess for R version 4.0.0 (2020-04-24)\n')
print(arguments)
stop("just checking arguments")

## for testing only 
setwd("/Users/ewilkie/Documents/Work/CCI/Polygenic/nimpress/preprocess")
arguments <- list()
arguments$GRCh37 = TRUE
arguments$GRCh38 = FALSE
#arguments$remove_blacklisted_regions = TRUE
arguments$blacklisted_regions_file = "/Users/ewilkie/Documents/Work/CCI/CCI_general_data_files/GRCh37_alldifficultregions.tier3.sorted.merged.sorted.bed"
arguments$remove_blacklisted_regions = FALSE
arguments$file = "Example/Example_File_to_process.csv"
arguments$LDproxy_pop ="GBR"
arguments$LDproxy_token ="cbe1b45bc8be"

##########################################
## testing that still needs to be done: ##
##########################################

## preprocess setup file needs to be completed
## implement error catching for blacklist file download

## double check that bedcove TRUE/FALSE does the right thing since I changed the coding of this variable 

## set up look to run master file list or composite of file??

## GRCh38 - but need different file for that - don't need to these ldproxy with this since it doesn't require 

## LDproxy

## is there a better way to do the token check than to run a test query and capture the error as string? 
## is LDproxy batch faster? - well be a mess to implement to check for which ones are kept and coverage
## to do about errors from LDproxy: 
## error: rs965506592 is not in 1000G reference panel.,
## error: rs334 is monoallelic in the GBR population.,

## from readME: "rsIDs which don't represent SNVs will be treated as unusable" -> has this been implemented? 

## might need to rearrange order of operation so that it errors are caught at the start and not after massive proprocessing 
## such as LDproxy flag = On but Bedfile=NULL

###################
###################
## Initial setup ##
###################
###################

## load processing functions files
source("Nimpress_preprocess_functions.R")

###############
## libraries ##
###############

## LDlinkR gets loaded later if used since it takes a while to install
message("[1/..] Loading libraries...")
pacman::p_load(data.table,GenomicRanges,rentrez,BSgenome.Hsapiens.UCSC.hg19)
## rentrez,bedr,LDlinkR,stringr,GenomicRanges

###################
## Genomic files ##
###################

message("[2/..] Setting up genomic environment...")

## get the correct assembly file
if(arguments$GRCh37 == TRUE){
  gv <- "GRCh37"
  assembly_file <- "Suppl/GCF_000001405.13_GRCh37_assembly_report.txt"
}else if (arguments$GRCh38 == TRUE){
  gv <- "GRCh38"
  assembly_file <- "Suppl/GCF_000001405.26_GRCh38_assembly_report.txt"
}

ass <- read.table(assembly_file, header=F, sep="\t",stringsAsFactors = F)
ass_sub <- ass[,c(3,7)]
assembly <- ass_sub[grep("NC_", ass_sub[,2]),]
colnames(assembly) <- c("CHR", "NC_CHR")
#print(assembly)

#######################################
## setup blacklisted genome bed file ##
#######################################

if(arguments$remove_blacklisted_regions == TRUE ){
  ## add error catching for blacklist file download - this is in preprocessing script
  #bedfile <- path to dowloaded bed
  #bedfile_to_Granges()
}else if(!is.null(arguments$blacklisted_regions_file)){
  ## hardcode file in for now
  bedfile <- arguments$blacklisted_regions_file
  gr <- bedfile_to_Granges(bedfile)
}else{
  bedfile = NULL
}

#############
## LDproxy ##
#############

## check if LDproxy is on
LDproxy_flag = "OFF"
if(!is.null(arguments$LDproxy_pop) & !is.null(arguments$LDproxy_token) & arguments$GRCh38 == TRUE){
  message("Warning: LDproxy does not support genome version GRCh38 and therefore LDproxy will be disabled")
} else if(!is.null(arguments$LDproxy_pop) & !is.null(arguments$LDproxy_token) && arguments$GRCh37 == TRUE){
  message("[3/..] Testing LDproxy parameters... ")
  pacman::p_load(LDlinkR)
  ## to check LDproxy pop input - get all available pops from package
  pop <- list_pop()
  bgpc <- pop$pop_code
  ## check valid background population
  if(arguments$LDproxy_pop %!in% bgpc){
    stop(paste(arguments$LDproxy_pop, " is not a valid background population. Select one from: ", paste(bgpc, collapse=",") , sep="" ))
  }
  ## test query
  testproxy <- LDproxy("rs456", arguments$LDproxy_pop, "r2", token = arguments$LDproxy_token)
  if(testproxy[1,1] != "  error: Invalid or expired API token. Please register using the API Access tab: https://ldlink.nci.nih.gov/?tab=apiaccess,"){
    message("LDproxy parameters ok")
    LDproxy_flag = "ON"
  }
}


#################
## Input files ##
#################

message("[4/..] Reading master file...")

master_file <- read.table(arguments$file, sep=",", header=T)
master_file.list <- split(master_file, seq(nrow(master_file)))


#####################
#####################
## File processing ##
#####################
#####################

message("[5/..] Starting file processing...")

## set up loop when single run is finished
#for(f in 1:length(master_file.list)){
#}
f <- 1

message(paste("[", 5+f ,"/..] Processing file: ", master_file.list[[f]]$GWAS_summary_statistic_file_and_path, sep="" ))
input <- master_file.list[[f]]$GWAS_summary_statistic_file_and_path

## check and format input file
gf_ok <- check_gwas_file(input)
#print(gf_ok)


## get indext of rsID column
rsID_ind <- grep("rsID", colnames(gf_ok))
## extract rsID and make unique
rsID <- gf_ok[,rsID_ind]
rsIDu <- as.vector(unique(rsID))
  
##########################################
## get rdID genomic location from dbSNP ##
##########################################

## list containing "rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE" - for Risk allele check
## this loop takes about 28.85696 / 13 = 2.219766 secs per SNPS
#ela <- Sys.time()

rsID_loc <- list()
for (rsid in 1: length(rsIDu)){
  message(paste("Getting info for : ", rsIDu[rsid],sep=""))
  res <- getrsID_info(rsIDu[rsid])
  rsID_loc[[rsid]] <- res
}

#ela <- Sys.time() - ela
#print(ela)

rsID_loc_df <- do.call(rbind, rsID_loc)
#print(rsID_loc_df)

##################################
## Check REF allele with genome ##
##################################

## create range of snps
snp_gr <- GRanges(seqnames=paste("chr", as.numeric(rsID_loc_df$CHR), sep=""), ranges=IRanges(start=as.numeric(rsID_loc_df$START)), starts.in.df.are.0based=F)

## results returned are in forward strand
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, snp_gr)
seqs_df <- as.data.frame(seqs)
colnames(seqs_df) <- "SEQ"

SNP_info <- cbind(rsID_loc_df, seqs_df)

## check for strand and ambigiouty between REF and ALT
strand <- list()
ambig <- list()

stp <- strsplit(SNP_info$ALT.ALLELE, "|")
for(rs in 1:nrow(SNP_info)){
  ## check whether REF is forwards or reverse strand (same as SEQ)
  if(SNP_info[rs,"REF.ALLELE"] == SNP_info[rs,"SEQ"]){
    strand[[rs]] <- "+"
  } else if(SNP_info[rs,"REF.ALLELE"] == complement(SNP_info[rs,"SEQ"])){
    strand[[rs]] <- "-"
  } else{
    strand[[rs]] <- NA
  }
  
  if(length(stp[[rs]]) == 1 && complement(SNP_info[rs,"REF.ALLELE"]) == stp[[rs]]){
    ambig[[rs]] <- "Y"
  }else if (length(stp[[rs]]) > 1 && complement(SNP_info[rs,"REF.ALLELE"]) %in% stp[[rs]]){
    ambig[[rs]] <- "D"
  }else{
    ambig[[rs]] <- "N"
  }
}

SNP_is <- cbind(SNP_info, strand=unlist(strand), ambiguous= unlist(ambig))

##############################
## Check for Ambiguous SNPS ##
##############################

## ambiguity is defined as when there are multiple plausible scenarios.
## for example:
#rs334  11  5248232          T      A|C|G   T      +           A
# is A the ALT or a complement of the REF? 
# rs4948492  10 63719739          C        G|T   C      +           C
# is C the ref or the complement of ALT

# rs13034020   2 61043834          A          G   A      +           T
# in this scenario the Risk T is clearly the complement of REF

## Check whether Risk ALlele is REF or ALT or neither (wrong)/ Flipped 

rsm <- merge(SNP_is, gf_ok[,1:2], by="rsID")
stp2 <- strsplit(rsm$ALT.ALLELE, "|")

risk_type <- list()
flipped <- list()

for(rll in 1:nrow(rsm)){
  ## if REF and ALT are ambiguous, then the risk_type and flipped are NA 
  if(rsm[rll, "ambiguous"] == "Y"){
    risk_type[[rll]] = NA
    flipped[[rll]] = NA
  }
  ## if ambigious flag is D or N - first check for ambigious then for check type and flipped
  else{
    if(rsm[rll,"Risk_allele"] == complement(rsm[rll,"REF.ALLELE"]) & rsm[rll,"Risk_allele"] %in% stp2[[rll]]){
      risk_type[[rll]] = NA
      flipped[[rll]] = NA
      rsm[rll,"ambiguous"] = "Y"
    }
    # Risk == REF and  == complemnt(ALT)
    # same as complment(risk) == complment(REF) & ALT
    else if(complement(rsm[rll,"Risk_allele"]) == complement(rsm[rll,"REF.ALLELE"]) & complement(rsm[rll,"Risk_allele"]) %in% stp2[[rll]]){
      risk_type[[rll]] = NA
      flipped[[rll]] = NA
      rsm[rll,"ambiguous"] = "Y"
    }
    else if(rsm[rll,"Risk_allele"] == rsm[rll,"REF.ALLELE"]){
      risk_type[[rll]] = "REF"
      flipped[[rll]] = "N"
    }else if (rsm[rll,"Risk_allele"] %in% stp2[[rll]]){
      risk_type[[rll]] = "ALT"
      flipped[[rll]] = "N"
      ## complement means flipped
    }else if(complement(rsm[rll,"Risk_allele"]) == rsm[rll,"REF.ALLELE"]){
      risk_type[[rll]] = "REF"
      flipped[[rll]] = "Y"
    }else if (complement(rsm[rll,"Risk_allele"]) %in% stp2[[rll]]){
      risk_type[[rll]] = "ALT"
      flipped[[rll]] = "Y"
    }
  }
}

SNP_fin <- cbind(rsm, risk_type=unlist(risk_type), flipped=unlist(flipped))
#SNP_fin[order(SNP_fin$ambiguous),]


##################
## Get coverage ##
##################

## can speed this up by extracting unique and the merging? Might not make much difference 
##1.623069 / 12 = 0.13 secs per SNP

#ela <- Sys.time()
if(bedfile == "NULL"){
  bedcov <- FALSE
  urercov <- cbind(SNP_fin, bedcov)
}else{
  bedcov <- vector()
  for (s in 1:nrow(SNP_fin)){
    nc <- get_cov(SNP_fin[s,])
    bedcov <- c(bedcov, nc)
  }
  urercov <- cbind(SNP_fin, bedcov)
}

#ela <- Sys.time() - ela
#print(ela)


#urercov[order(urercov$ambiguous,urercov$bedcov),]
########!!!!!!!! since cov file type has changed FALSE now means that no removal - check this correct 
## updated input file contains some FALSE and SOME TRUE 
#urercov[, "bedcov"] <- TRUE

#############
## LDproxy ##
#############

## what to do about: error: rs965506592 is not in 1000G reference panel.,

## ldproxy takes 8.791638 mins for 10 SNPS, so a bit under a minute for 1


## if no bedfile all bedcov will equal FALSE, if no fall in region no sub needed
if(length(unique(urercov$bedcov)) == 1 & unique(urercov$bedcov) == FALSE){
  ## write output
  
  #LDproxy_out <- rsID CHR.x  START.x REF.ALLELE ALT.ALLELE bedcov RSID_Proxy CHR.y  START.y  REF  ALT
  
## if LDproxy not used and bedfile provided, remove those without coverage. 
}else if (LDproxy_flag == "OFF" & bedfile != "NULL"){
  message("Remove blacklisted regions enabled, but LDproxy disabled...")
  message("rsID without coverage will be removed")
  message("If however you would like those rsIDs without coverage to be substituted, add --remove_blacklisted_regions in commandline parameters")
  #-------> LOG FILE those that are removed should be put into "log file"
  
  #LDproxy_out <- rsID CHR.x  START.x REF.ALLELE ALT.ALLELE bedcov RSID_Proxy CHR.y  START.y  REF  ALT
  
## if ldproxy is on but bedile not provided, don't need to do resub -- however LDproxy only needs to be done for thise falling in the bed region. So with bedfile both sub and resub are needed without bedfile don't know which to sub or whether its necesery so put an error flag to say, that LDproxy requires bedfile 
}else if (LDproxy_flag == "ON" & bedfile == "NULL"){
  message("LDproxy required bedfile for sub")
  #LDproxy_out <- rsID CHR.x  START.x REF.ALLELE ALT.ALLELE bedcov RSID_Proxy CHR.y  START.y  REF  ALT
    
## do ldproxy and exlude those that fall in the bed regions  - do sub and resub
}else if (LDproxy_flag == "ON" & bedfile != "NULL"){
  
  ela <- Sys.time()
  
  ## only run LDproxy on unique results for those that don't have coverage or are ambigious
  run_ind <- which(urercov$bedcov == TRUE & urercov$ambiguous != "Y")
  ldproxy_input <- unique(urercov[run_ind,"rsID"])
  ## to keep track of which RSids are already in DF, without expanding loop 
  SNP_kept <- unique(urercov[-run_ind,"rsID"])
  
  ldproxy_ls <- list()
  for(s in 1:length(ldproxy_input)){
    ## need to implement collapse of duplicate alt alleles - 
    ldproxy_res <- getLDproxy(ldproxy_input[s], arguments$LDproxy_pop, arguments$LDproxy_token, SNP_kept)
    
    ## if results is NA add it anyways for later removal 
    if(is.na(ldproxy_res$RSID_Proxy)[1]){
      ldproxy_res_keep <- ldproxy_res
    }else{
      ## check which not in dataset 
      wni <- which(ldproxy_res$RSID_Proxy %!in% SNP_kept)
      if(length(wni) == 0){
        ldproxy_res_keep <- LDproxy_NA_res(ldproxy_input[s])
      }else{
        ## keep one result only that is not already in dataset
        ldproxy_res_keep <- ldproxy_res[wni[1],]
        ## this is so that no duplicates appear in the data
        SNP_kept <- c(SNP_kept,ldproxy_res_keep$RSID_Proxy)
      }
    }
    
    ## combine ldproxy_res with orignal input - automatically only uses bedcov == TRUE data
    mres <- merge(urercov, ldproxy_res_keep, by.x="rsID", by.y="RSID_input")
    ldproxy_ls[[s]] <- mres
  }
  
  #ela <- Sys.time() - ela
  
  LDproxy_df <- do.call(rbind, ldproxy_ls)
  LDproxy_in <- LDproxy_df[!is.na(LDproxy_df$RSID_Proxy),]
  
  ## format results that don't have bedcov
  nocov_padd <- cbind(urercov[urercov$bedcov == FALSE,],NA,NA,NA,NA,NA)
  colnames(nocov_padd) <- c("rsID", "CHR.x", "START.x", "REF.ALLELE.x", "ALT.ALLELE.x", "SEQ", "strand", "ambiguous", "Risk_allele", "risk_type", "flipped", "bedcov", "RSID_Proxy", "CHR.y", "START.y", "REF.ALLELE.y", "ALT.ALLELE.y")
  
  ## final output
  LDproxy_out <- rbind(nocov_padd, LDproxy_in)
}

## to get colnames 
#pz <- paste0('"', paste(colnames(LDproxy_in), collapse='", "'), '"')
#print(pz, quote=F)


############################
## DEFINE CORRECT ALLELES ##
############################


