library(tidyverse)
library(readr)

# read metadata.csv with only "Sample_ID" and "*_sr" features
metadata <- read_csv("./metadata.csv") %>% 
  select(Sample_ID,(matches("_sr$")))

# read antibiotic resistance associated unitigs
unitigs_azm <- read_table2("azm_sr_gwas_filtered_unitigs.Rtab")
unitigs_cfx <- read_table2("cfx_sr_gwas_filtered_unitigs.Rtab")
unitigs_cip <- read_table2("cip_sr_gwas_filtered_unitigs.Rtab")

#Function to transpose unitigs table, join metadata, then create .csv file
unitigs <- function(unitig){
  
  # deparse + substitute combo to extract object name
  unitig_name <- deparse(substitute(unitig))
  
  '%>%' <- tidyr::'%>%'
  
  unitig %>%
    # transpose
    tidyr::pivot_longer(-pattern_id, "Sample_ID", "value") %>% 
    tidyr::pivot_wider(Sample_ID, pattern_id) %>%
    # join relevant antibiotic metadata
    dplyr::inner_join(dplyr::select(metadata, "Sample_ID", dplyr::matches(stringr::str_extract(unitig_name, "[[:alpha:]]+$"))), by = "Sample_ID") %>% 
    dplyr::select('Sample_ID', dplyr::last_col(), dplyr::everything()) %>% 
    # rename unitigs to format "geneX" where X is a positive integer
    data.table::setnames(., old = names(.[3:length(.)]), new = paste0('gene', seq_along(.[3:length(.)])), skip_absent = T) %>% 
    tidyr::drop_na() %>% 
    # create .csv file
    readr::write_csv(path=paste0('./', stringr::str_extract(unitig_name,"[[:alpha:]]+[_][[:alpha:]]+"), '.csv'))
}

# Manual Transpose, join, and create .csv files of azm unitigs
# unitigs_azm %>%
#   pivot_longer(-pattern_id, "Sample_ID", "value") %>%
#   pivot_wider(Sample_ID, pattern_id) %>%
#   inner_join(select(metadata, "Sample_ID", matches('^azm')), by = "Sample_ID") %>%
#   select('Sample_ID', last_col(), everything()) %>%
#   data.table::setnames(., old = names(.[2:length(.)]), new = paste0('gene', seq_along(.[2:length(.)])), skip_absent = T) %>% 
#   drop_na() %>%
#   write_csv(path="./unitigs_azm.csv")

# Use 'unitigs' function to create .csv files
unitigs(unitigs_azm)
unitigs(unitigs_cfx)
unitigs(unitigs_cip)
