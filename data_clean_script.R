library(tidyverse)
library(readr)
# read metadata.csv with only "Sample_ID" and "*_sr" features
metadata <- read_csv("./metadata.csv") %>% 
  select(Sample_ID,(matches("_sr$")))

# read antibiotic resistance associated unitigs
unitigs_azm <- read_table2("azm_sr_gwas_filtered_unitigs.Rtab")
unitigs_cfx <- read_table2("cfx_sr_gwas_filtered_unitigs.Rtab")
unitigs_cip <- read_table2("cip_sr_gwas_filtered_unitigs.Rtab")

detach("package:readr", unload = T)

#Function to transpose and create .csv files of unitigs
unitig_trans <- function(unitig){
  # requires loads packages necessary for function
  lapply(c("magrittr", "tidyr", "readr", "stringr"), require, character.only = TRUE)
  
  # deparse + substitute combo to extract object name
  unitig_name <- deparse(substitute(unitig))
  
  unitig %>% 
    tidyr::pivot_longer(-pattern_id, "Sample_ID", "value") %>% 
    tidyr::pivot_wider(Sample_ID, pattern_id) %>%
    readr::write_csv(path=paste0('./', stringr::str_extract(unitig_name,"[[:alpha:]]+[_][[:alpha:]]+"), '.csv'))
}

# # Transpose and create .csv file of "azm_sr" unitigs
# unitigs_azm %>%
#   pivot_longer(-pattern_id, "Sample_ID", "value") %>%
#   pivot_wider(Sample_ID, pattern_id) %>%
# write_csv(path="./unitigs_azm.csv")