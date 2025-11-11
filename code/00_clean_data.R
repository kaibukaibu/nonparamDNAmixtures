# Code for cleaning 1p and mixture data for paper "Non-parametric simulation of DNA mixture profiles from one-person trace profiles"
# We shall turn the data into a format needed for DNAmixtures and remove replicate alleles (for some profiles, there are some alleles
# present which have 2 or more peak height data points, we only choose the highest peak).
# Author: Kai Budrikas (kaib@itu.dk)





# -----------------------------------------------------------------------------------------------

# Load packages and set data locations

library(tidyverse)
library(readxl)


setwd("C:/Users/kaib/OneDrive - ITU/Documents/GitHub/repDNAmixtures/") #Project location
input_provedit <- "data/data_provedit/PROVEDIt_1-5-Person CSVs UnFiltered_3500_GF29cycles/" # Provedit mixtures folder
output_folder <- "data/data_provedit_cleaned/" # Folder for saving the results, must end with "/"
source("code/helping_functions.R") # Some functions for transforming the profile data from long to wide format and vice versa




# -----------------------------------------------------------------------------------------------




########## GENOTYPES ################################################

genotypes00 <- read_xlsx(paste0(input_provedit, "PROVEDIt_RD14-0003 GF Known Genotypes.xlsx")) %>% 
  select(-"Research ID") %>% 
  rename(ID = "Sample ID") %>% 
  mutate(ID = paste0("c", ID)) %>% 
  pivot_longer(cols=2:25, names_to = "Marker", values_to = "Alleles") %>% 
  filter(!(Marker %in% c("AMEL", "Yindel", "DYS391"))) %>% 
  mutate(Allele1 = str_split_i(Alleles, ",", 1),
         Allele2 = str_split_i(Alleles, ",", 2)) %>% 
  select(-Alleles)


write.csv(genotypes00,
          paste0(output_folder, "genotypes.csv"), row.names=F, quote = F)


# Save also as single files for DNAStatistX

for(contr in unique(genotypes00$ID)){
  
}
















