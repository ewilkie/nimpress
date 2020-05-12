######################################
## Nimpress preprocessing functions ##
######################################

## Not in function
'%!in%' <- function(x,y)!('%in%'(x,y));

###########################
## processing input file ##
###########################

check_gwas_file <- function(input){
  
  gwas_file <- read.table(input, sep=",", header=T, colClasses=c("character","character","numeric", "numeric", "numeric"))
  
  ## remove blank rows
  blank <- which(gwas_file[,1] == "")
  na <- which(is.na(gwas_file[,1]))
  rm_blank <- c(blank, na)
  if(length(rm_blank) > 0){
    gwas_file <- gwas_file[-rm_blank,]
  }
  
  ## remove blank columns
  if(ncol(gwas_file) > 5){
    gwas_file <- gwas_file[,-(6:ncol(gwas_file))]
  }
  
  ## remove all leading and trailing blank spaces
  gwas_file <- as.data.frame(apply(gwas_file,2,function(x)gsub('\\s+', '',x)))
  
  ## check if all rsIDs are in propper format
  sw <- startsWith(gwas_file$rsID, "rs")
  sww <- which(sw !=TRUE)
  if(length(sww) == 0){
    message("All rsID ok")
  }else{
    stop(paste("Line ", sww, " does not contain a valid rsID", sep=""))
  }
  
  ## check if valid nucleotides
  nva <- which(gwas_file$Risk_allele %!in% c("A", "T", "G", "C"))
  if(length(nva) == 0){
    message("All Risk_allele ok")
  }else{
    stop(paste("Line ", nva, " does not contain a valid risk allele", sep=""))
  }
  
  ## Format Effect size
  or.ind <- grep("OR", colnames(gwas_file))
  beta.ind <- grep("Beta", colnames(gwas_file))
  
  if(length(or.ind) == 1){
    OR <- as.vector(gwas_file[,or.ind])
    gwas_file[,or.ind] <- log(as.numeric(sub('\\(.*', '', OR)))
    colnames(gwas_file)[or.ind] <- "Effect.size"
    message("Coverting OR to BETA")
  }else if(length(beta.ind) == 1){
    colnames(gwas_file)[beta.ind] <- "Effect.size"
  }
  
  return(gwas_file)
}  

##########################################
## get rdID genomic location from dbSNP ##
##########################################

## extracts information from rsID dbSNP lookup: CHR, start, REF.Allele, ALT.Allele
format_dbSNP <- function(x){
  ## extract assembly, genome position and variant details 
  g1 <- gsub("^[A-Za-z]*=", "", x)
  ## remove everything after "|"
  g2 <- strsplit(gsub("\\|.*$", "", g1), ",")
  ## split chr from location
  g3 <- unlist(strsplit(g2[[1]], ":"))
  ## format into df
  odd <- seq_along(g3) %% 2 == 1
  NC_CHR <- g3[odd]
  nt <- gsub("^g\\.", "", g3[!odd])
  START <- gsub("[A-Z]>[A-Z]", "", nt)
  nt2 <- unlist(lapply(strsplit(nt,"[0-9]"),function(x) x[length(x)]))
  nt3 <- do.call(rbind, strsplit(nt2, ">"))
  colnames(nt3) <- c("REF", "ALT")
  df <- data.frame(NC_CHR,START, nt3)
  df2 <- unique(merge(assembly,df,by="NC_CHR"))
  df3 <- df2[,-1]
  return(df3)
}

## this function lookes up information in the dbSNP database for each rsID
## same results as a list it handle if multiple ALT alleles exist

getrsID_info <- function(rsid_input){
  snp_term <- paste(rsid_input, "[RS]", sep="")
  r_search <- entrez_search(db="snp", term=snp_term)
  ## if no result returned
  if(length(r_search$id) == 0){
    final_snp <- cbind(rsid_input, NA,NA,NA,NA)
    colnames(final_snp) <- c("rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE")
  }else{
    ## info associated with snp_id
    multi_summs <- entrez_summary(db="snp", id=r_search$id)
    uid <- unique(extract_from_esummary(multi_summs, c("snp_id","snp_class")))
    all_recs <- entrez_fetch(db="snp", id=uid, rettype="xml")
    tax_list <- XML::xmlToList(all_recs)
    
    if(tax_list$DocumentSummary$SNP_CLASS == "snv"){
      dbSNP_res <- format_dbSNP(tax_list$DocumentSummary$DOCSUM)
      ## multiple Alt Alleles are on different lines 
      final_snp <- cbind(rsid_input, dbSNP_res)
      colnames(final_snp) <- c("rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE")
    ## if rsID doesn't represent SNV
    }else{
      final_snp <- cbind(rsid_input, NA,NA,NA,NA)
      colnames(final_snp) <- c("rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE")
    }
  }
  return(final_snp)
}



###################################################
## Get coverage with blacklisted bad file region ##
###################################################

get_cov <- function(snp_info){
  snp <- GRanges(seqnames=as.numeric(snp_info$CHR), ranges=IRanges(start=as.numeric(snp_info$START), end=as.numeric(snp_info$START)+1),starts.in.df.are.0based=TRUE)
  hits <- findOverlaps(gr,snp)
  if(length(hits@from) > 0){
    bedcov = TRUE
  }else{
    bedcov = FALSE
  }
  return(bedcov)
}
