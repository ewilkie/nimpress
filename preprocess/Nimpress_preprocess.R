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
## is LDproxy batch faster? 
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
pacman::p_load(data.table,GenomicRanges,rentrez)
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
f <- 2

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

##################
## Get coverage ##
##################

## can speed this up by extracting unique and the merging? Might not make much difference 
##1.623069 / 12 = 0.13 secs per SNP

#ela <- Sys.time()
if(bedfile == "NULL"){
  bedcov <- FALSE
  urercov <- cbind(rsID_loc_df, bedcov)
}else{
  bedcov <- vector()
  for (s in 1:nrow(rsID_loc_df)){
    nc <- get_cov(rsID_loc_df[s,])
    bedcov <- c(bedcov, nc)
  }
  urercov <- cbind(rsID_loc_df, bedcov)
}

#ela <- Sys.time() - ela
#print(ela)

########!!!!!!!! since cov file type has changed FALSE now means that no removal - check this correct 
## updated input file contains some FALSE and SOME TRUE 
#urercov[, "bedcov"] <- TRUE

#############
## LDproxy ##
#############

## what to do about: error: rs965506592 is not in 1000G reference panel.,

## ldproxy takes 8.791638 mins for 10 SNPS, so a bit under a minute for 1

snp <- "rs965506592"
pop <- arguments$LDproxy_pop
token <- arguments$LDproxy_token
n <- 1
snp <- as.vector(urercov[n,"rsID"])

snp <- "rs3731217"


## if no bedfile all bedcov will equal FALSE, if no fall in region no sub needed
if(length(unique(urercov$bedcov)) == 1 & unique(urercov$bedcov) == FALSE){
  ## write output
  
## if LDproxy not used and bedfile provided, remove those without coverage. 
}else if (LDproxy_flag == "OFF" & bedfile != "NULL"){
  message("Remove blacklisted regions enabled, but LDproxy disabled...")
  message("rsID without coverage will be removed")
  message("If however you would like those rsIDs without coverage to be substituted, add --remove_blacklisted_regions in commandline parameters")
  #-------> LOG FILE those that are removed should be put into "log file"
  ## need to know output format for downstream
  
## if ldproxy is on but bedile not provided, don't need to do resub -- however LDproxy only needs to be done for thise falling in the bed region. So with bedfile both sub and resub are needed without bedfile don't know which to sub or whether its necesery so put an error flag to say, that LDproxy requires bedfile 
}else if (LDproxy_flag == "ON" & bedfile == "NULL"){
  message("LDproxy required bedfile for sub")

    
## do ldproxy and exlude those that fall in the bed regions  - do sub and resub
}else if (LDproxy_flag == "ON" & bedfile != "NULL"){
  
  ## need to expand getLDproxy and ldpoxy_res_keep to take into account bedfile coverage - resub
  ela <- Sys.time()
  
  ## only run LDproxy on unique results. 
  ldproxy_input <- unique(urercov$rsID)
  ## to keep track of which RSids are already in DF, without expanding loop 
  SNP_kept <- ldproxy_input
  ldproxy_ls <- list()
  for(s in 1:length(ldproxy_input)){
    ## last argument is to retain this snps that are already in the list
    ldpoxy_res <- getLDproxy(ldproxy_input[s], arguments$LDproxy_pop, arguments$LDproxy_token, SNP_kept)
    ## if result is not NA and account for results length of 1 and more
    if(!is.na(ldpoxy_res$RSID_Proxy[1])){
      ## keep one result only that is not already in dataset
      ldpoxy_res_keep <- ldpoxy_res[which(ldpoxy_res$RSID_Proxy %!in% ldproxy_input)[1],]
      ## this is so that no duplicates appear in the data
      SNP_kept <- c(SNP_kept,ldpoxy_res_keep$RSID_Proxy)
    }else{
      ldpoxy_res_keep <- ldpoxy_res
    }
    ##### combine ldproxy_res with orignal input 
    mres <- merge(urercov, ldpoxy_res_keep, by.x="rsID", by.y="RSID_input")
    ## remove duplicate alt alleles
    ldproxy_ls[[s]] <- mres
  } 
  ela <- Sys.time() - ela

  do.call(rbind, ldproxy_ls)
  
}
  
########################## FROM OLD CODE --- DELETE once used ##############################


## from old function. Undecided what to do with it   
for(i in 1:length(coord_check)){
  if(arguments$bed == "NULL"){
    getALT <- getrsID_info(my_proxies_keep3[i,1])
    if(!is.na(getALT)){
      ALTcol <- paste(sort(as.vector(getALT$ALT.ALLELE)), collapse=",")
      Alleles <- paste("(", unique(as.vector(getALT$REF.ALLELE)) , "/", ALTcol, ")", sep="")
      new_rd <- c(as.matrix(my_proxies_keep3[i,1:2]), Alleles)
      rsid_keep <- c(rsid_keep,as.vector(my_proxies_keep3[i,"RS_Number"]))
      break
    }else{
      next
    }
  }else{
    coord <- coord_check[i]
    snp_chr <- sub("chr","" , sub('\\:.*', '', coord))
    start <- as.numeric(sub("\\-.*", "", sub('.*\\:', '', coord)))
    snp <- GRanges(seqnames=snp_chr, ranges=IRanges(start=start, end=start))
    hits <- findOverlaps(gr,snp)
    if(length(hits@from) > 0){
      if(as.vector(my_proxies_keep3[i,"RS_Number"]) %!in% rsid_keep){
        ## need to get all alternative allelse fot this new rsID 
        getALT <- getrsID_info(as.vector(my_proxies_keep3[i,1]))
        if(!is.na(getALT)){
          ALTcol <- paste(sort(as.vector(getALT$all.alleles$ALT.ALLELE)), collapse=",")
          Alleles <- paste("(", unique(as.vector(getALT$all.alleles$REF.ALLELE)) , "/", ALTcol, ")", sep="")
          new_rd <- c(as.matrix(my_proxies_keep3[i,1:2]), Alleles)
          rsid_keep <- c(rsid_keep,as.vector(my_proxies_keep3[i,"RS_Number"]))
          break
        }else{
          next
        }
      }
    }else{
      new_rd <- NA
    } 
  }
}




## just for testting
urercov[1,"bedcov"] <- TRUE

rsid_keep <- vector()
LDres <- list()
for(n in 1:nrow(urercov)){
  if(urercov[n,"bedcov"] == "TRUE"){
    LD <- getLDproxy(as.vector(urercov[n,"rsID"]),arguments$LDproxy_pop, arguments$LDproxy_token)
    if (length(LD) == 3){
      coord <- strsplit(as.vector(LD[2]), ":")[[1]]
      j <- as.vector(LD[3])
      all <- regmatches(j, gregexpr("(?<=\\().*?(?=\\))", j, perl=T))[[1]]
      ref <- strsplit(all, "/")[[1]][1]
      alt <- strsplit(all, "/")[[1]][2]
      res <- c(as.vector(LD[1]), gsub("chr", "", coord[1]), coord[2], ref, alt)
      LDres[[n]] <- res
      rsid_keep <- c(rsid_keep, as.vector(LD[1]))
    }else if(is.na(LD)) {
      LDres[[n]] <- NA
    }
  }else{
    LDres[[n]] <- NA
  }
}


LDres_mat <- do.call(rbind, LDres)
if(ncol(LDres_mat) == 5){
  colnames(LDres_mat) <- c("ALT.rsID", "ALT.chr", "ALT.start","REF.new", "ALT.new")
}else{
  ## when there are no LDproxy results, either all cov == TRUE or when cov == FALSE, but no proxy exists 
  ## this is so code below doesn't break
  LDres_mat <- cbind(LDres_mat, LDres_mat,LDres_mat,LDres_mat,LDres_mat)
  colnames(LDres_mat) <- c("ALT.rsID", "ALT.chr", "ALT.start","REF.new", "ALT.new")
}

## combine original rsID and alternative
comb <- as.data.frame(cbind(urercov,LDres_mat))









