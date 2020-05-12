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

## this function lookes up information in the dbSNP database for each rsID
## extracts information for GRCh37: CHR, start, REF.Allele, ALT.Allele
## same results as a list it handle if multiple ALT alleles exist


getrsID_info <- function(rsid_input){
  
  snp_term <- paste(rsid_input, "[RS]", sep="")
  r_search <- entrez_search(db="snp", term=snp_term)
  
  
  if(length(r_search$id) == 0){
    final_snp <- cbind(rsid_input, NA,NA,NA,NA)
    colnames(final_snp) <- c("rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE")
  }else{
    multi_summs <- entrez_summary(db="snp", id=r_search$id)
    
    ## get the unique SNP_ID
    uid <- unique(extract_from_esummary(multi_summs, c("snp_id")))
    all_recs <- entrez_fetch(db="snp", id=uid, rettype="xml")
    tax_list <- XML::xmlToList(all_recs)
    
    ## extract assembly, genome position and variant details 
    g1 <- gsub("^[A-Za-z]*=", "", tax_list$DocumentSummary$DOCSUM)
    g2 <- strsplit(g1,"\\|")
    g3 <- as.data.frame(strsplit(g2[[1]][1],",")[[1]])
    colnames(g3) <- "ID"
    
    g4 <- apply(g3, 2, function(x) strsplit(x,":g\\."))
    g5 <- as.data.frame(do.call(rbind, g4$ID),stringsAsFactors=FALSE)
    colnames(g5) <- c("assembly", "START")
    
    g6 <- strsplit(g5$START,"[0-9]")
    g7 <- unlist(lapply(g6,function(x) x[length(x)]))
    g8 <- do.call(rbind, strsplit(g7, ">"))
    g9 <- cbind(g5,g8)
    rmg <- gsub("[^0-9]", "", g9$START)
    g9$START <- rmg
    
    ## can't remember why this is in here
    #if(!is.null(g9$g8)){
    #  g9 <- g9[-which(g9$g8 == "del"),]
    #}
    
    ## check if there are no entires for that rsID
    g10 <- g9
    colnames(g10)[3:4] <- c("REF_Allele","ALT_Allele")
    
    ## subset to only NC
    g11 <- g10[grep("^NC",  g10[,1]),]
    
    ## extract CHR
    g12 <- gsub("NC_0+","",g11[,1])
    CHR <- gsub("\\.[0-9]*", "", g12)
    
    g13 <- data.frame(g11$assembly, CHR, g11$START, g11$REF_Allele, g11$ALT_Allele)
    
    ## Get right Assembly from SNP
    inter <- intersect(g13[,1], assembly[,2])
    g14 <- unique(g13[grep(inter,g13[,1]),])
    g15 <- g14[,-1]
    
    ## multiple Alt Alleles are on different lines 
    final_snp <- cbind(rsid_input, g15)
    colnames(final_snp) <- c("rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE")
  }
  return(final_snp)
}
