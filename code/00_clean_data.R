# Code for cleaning 1p and mixture data for paper "Non-parametric simulation of DNA mixture profiles from one-person trace profiles"
# We shall turn the data into a format needed for DNAmixtures and remove replicate alleles (for some profiles, there are some alleles
# present which have 2 or more peak height data points, we only choose the highest peak).
# Author: Kai Budrikas (kaib@itu.dk)

# TO DO: adjust other files according to this



# -----------------------------------------------------------------------------------------------

# Load packages and set data locations

library(tidyverse)
library(readxl)


setwd("..") #Project location
input_provedit <- "data/data_provedit/PROVEDIt_1-5-Person CSVs UnFiltered_3500_GF29cycles/" # Provedit mixtures folder
output_folder <- "data/data_provedit_cleaned/" # Folder for saving the results, must end with "/"
source("code/helping_functions.R") # Some functions for transforming the profile data from long to wide format and vice versa




# -----------------------------------------------------------------------------------------------




########## GENOTYPES ################################################

genotypes00 <- read_xlsx(paste0(input_provedit, "PROVEDIt_RD14-0003 GF Known Genotypes.xlsx")) %>% 
  select(-"Research ID") %>% 
  rename("SampleName" = "Sample ID") %>% 
  mutate(SampleName = paste0("c", SampleName)) %>%
  rename("Sample Name" = SampleName) %>% 
  pivot_longer(cols=2:25, names_to = "Marker", values_to = "Alleles") %>% 
  filter(!(Marker %in% c("AMEL", "Yindel", "DYS391"))) %>% 
  mutate("Allele 1" = str_split_i(Alleles, ",", 1),
         "Allele 2" = str_split_i(Alleles, ",", 2)) %>% 
  select(-Alleles) %>% 
  arrange(Marker)


write.csv(genotypes00,
          paste0(output_folder, "genotypes.csv"), row.names=F, quote = F)


# Save also as single files for DNAStatistX

for(contr in unique(genotypes00$`Sample Name`)){
  subset <- genotypes00 %>% 
    filter(`Sample Name`== contr)
  
  write.csv(subset,
            paste0(output_folder, "genotypes/", contr, ".csv"), row.names=F, quote = F)
}












########## EVIDENCE PROFILES ##################################

all_traces00 <- list.files(input_provedit, full.names = T, recursive = T)
all_traces01 <- all_traces00[str_detect(all_traces00, "/15 sec/") & !str_detect(all_traces00, "/5-Person/")] #Choose only 15 sec ones and remove 5p

all_traces02 <- all_traces01 %>% 
  lapply(function(file) {
    read.csv(file) %>%
      mutate(across(everything(), as.character))
  }) %>%
  bind_rows() %>% 
  relocate(Sample.File, Marker, starts_with("Allele"), starts_with("Height")) %>% 
  
  filter(!(Marker %in% c("AMEL", "Yindel", "DYS391"))) %>% 
  select(contains(c("Sample.File", "Marker", "Allele", "Height"))) %>% 
  select_if(function(x) !(all(is.na(x)) | all(x==""))) %>%  #Drop empty columns
  rename(SampleName = Sample.File) %>%
  mutate(across(contains(c("Allele", "Height")), as.character)) %>% 
  
  #Bring to long format
  pivot_longer(
    starts_with(c("Allele", "Height")),
    cols_vary = "fastest",
    names_to = c(".value"),
    names_pattern = "(.)"
  ) %>% 
  rename(Allele = A,
         Height = H) %>% 
  
  #Drop some rows
  filter(!is.na(Allele) & Allele != "" & Allele != "OL") %>% 
  mutate(Allele = as.character(as.numeric(Allele)),
         Height = as.numeric(Height)) %>% 
  
  #Fix duplicates -- some alleles have two height values, we take only the highest peak
  group_by(SampleName, Allele, Marker) %>% 
  summarise(Height = max(Height)) %>% 
  ungroup() %>% 
  
  arrange(SampleName, Marker, Allele)

write.csv(all_traces02,
          paste0(output_folder, "traces.csv"), row.names=F, quote = F)








