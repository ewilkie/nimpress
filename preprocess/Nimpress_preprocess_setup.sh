##############################
## Nimpress prepocess setup ##
##############################

path="$(pwd)/Suppl"
##get blacklisted regions - place in Suppl folder
wget --directory-prefix=$path ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v2.0/GRCh37/union/GRCh37_alldifficultregions.bed.gz
unzip "${path}/GRCh37_alldifficultregions.bed.gz"







