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
arguments$LDproxy_pop ="GBR"
arguments$LDproxy_token ="cbe1b45bc8be"

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
pacman::p_load(data.table,GenomicRanges,rentrez, )
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
  gr <- makeGRangesFromDataFrame(ovlp, keep.extra.columns = TRUE,starts.in.df.are.0based=TRUE)
}else{
  bedfile = NULL
}


## check LDproxy input 
bgpc <- c("CHB","JPT","CHS","CDX","KHV","CEU","TSI","FIN","GBR","IBS","YRI","LWK","GWD","MSL","ESN","ASW","ACB","MXL","PUR","CLM","PEL","GIH","PJL","BEB","STU","ITU")
## check if LDproxy is on
LDproxy_flag = "OFF"
if(!is.null(arguments$LDproxy_pop) && !is.null(arguments$LDproxy_token)){
  message("[3/..] Testing LDproxy parameters... ")
  pacman::p_load(LDlinkR)
  ## check valid background population
  if(arguments$LDproxy_pop %!in% bgpc){
    stop(paste(arguments$LDproxy_pop, " is not a valid background population. Select one from: https://www.internationalgenome.org/faq/which-populations-are-part-your-study"))
  }
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

## load processing functions files
source("./Nimpress_preprocess_functions.R")

message("[5/..] Starting file processing...")

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

## can speed this up by extracting unique and the merging? Might not make much difference 
##1.623069 / 12 = 0.13 secs per SNP

ela <- Sys.time()
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

ela <- Sys.time() - ela
print(ela)

########!!!!!!!! since cov file type has changed FALSE now means that no removal

#############
## LDproxy ##
#############

## if no bed file and no LDSUP skip this step
## if ldproxy off and bedfiles is provided - remove those that fall inside the Bedfile region

## if no bedfile all bedcov will equal FALSE, if no fall in region no sub needed
if(length(unique(urercov$bedcov)) == 1 & unique(urercov$bedcov) == FALSE){
  ## write output
## if LDproxy not used and bedfile provided, remove those without coverage. 
}else if (LDproxy_flag == "OFF" & bedfile != "NULL"){
  #-------> LOG FILE those that are removed should be put into "log file"
## if ldproxy is on  
}else if (LDproxy_flag == "ON"){
  ## do LDproxy & Sub 
}




###################
## Dlinkpipeline ##
###################

## coords are for GRCh37

## certain rsIDs can be perfectly linked. 
## in this case the LDproxy rsID turned out to be the same
##https://ldlink.nci.nih.gov/?var1=rs4948492&var2=rs4245597&pop=GBR&tab=ldpair
## could potentially also arise due to different issues, but the solution is to keep track of the new rsIDs and if that is already obtained, find another one, if no exist, drop the original rsID
## how to do this?

getLDproxy <- function(snp){
  ## for breaking
  new_rd <- TRUE
  ## run Query
  my_proxies <- LDproxy(snp, pop = "GBR", r2d = "r2", token = "cbe1b45bc8be", file = FALSE)
  ## extract only those with R2 >= 0.9
  my_proxies_keep <- my_proxies[which(my_proxies$R2 >= 0.9),c(1,2,3)]
  ## remove those without rs number
  my_proxies_keep2 <- my_proxies_keep[grep("rs", my_proxies_keep$RS_Number),]
  ## remove those that are already in the dataset
  my_proxies_keep3 <- my_proxies_keep2[-which(my_proxies_keep2$RS_Number %in% original_SNP),]
  ## if no proxies
  if(nrow(my_proxies_keep3) == 0){
    new_rd <- NA
  }else{
    ## split coords
    coord <- sub(".*:(\\d+).*", "\\1", my_proxies_keep3$Coord)
    coord_check <- paste(my_proxies_keep3$Coord , "-" , coord, sep="")
    for(i in 1:length(coord_check)){
      if(arguments$bed == "NULL"){
        getALT <- getrsID_info(my_proxies_keep3[i,1])
        if(!is.na(getALT)){
          ALTcol <- paste(sort(as.vector(getALT$all.alleles$ALT.ALLELE)), collapse=",")
          Alleles <- paste("(", unique(as.vector(getALT$all.alleles$REF.ALLELE)) , "/", ALTcol, ")", sep="")
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
  }
  return(new_rd)
}  


rsid_keep <- vector()
LDres <- list()
for(n in 1:nrow(urercov)){
  if(urercov[n,4] == "FALSE"){
    LD <- getLDproxy(as.vector(urercov[n,3]))
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


