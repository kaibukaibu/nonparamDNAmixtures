# Code for applying stutter filter and analytical thresholds in paper "Non-parametric simulation of DNA mixture profiles from one-person trace profiles"
# Author: Kai Budrikas (kaib@itu.dk)




# -----------------------------------------------------------------------------------------------

library(tidyverse)
library(openxlsx)

setwd("..")
source("code/helping_functions.R") # Some functions to help
genotypes_folder <- "data/data_provedit/genotypes/" #Folder for genotypes, must have "/" at the end. Only used for making sure which markers to look at.
results_folder <- "data/data_replicated/" #Where to save the results, must have "/" at the end



# -----------------------------------------------------------------------------------------------




# Make folders for stutter filter + AT info
input_folders <- list.files(results_folder,
                            full.names = T)[!str_detect(list.files(results_folder), "sFRiman") &
                                              !str_detect(list.files(results_folder), "log") &
                                              !str_detect(list.files(results_folder), "realmix")]
for(file in input_folders){
  if (dir.exists(paste0(file, "_sFRiman"))) {
    unlink(paste0(file, "_sFRiman"), recursive = TRUE, force = T)  # Deletes the folder and all contents
  }
  dir.create(paste0(file, "_sFRiman"))
}
all_files <- lapply(input_folders, function(folder) list.files(folder, full.names=T)) %>% unlist()



# Add expected stutter
stutter_string <- "CSF1PO 8.77
D10S1248 11.46
D12S391 13.66
D13S317 9.19
D16S539 9.48
D18S51 12.42
D19S433 9.97
D1S1656 12.21
D21S11 10.45
D22S1045 16.26
D2S1338 11.73
D2S441 8.10
D3S1358 10.98
D5S818 9.16
D7S820 8.32
D8S1179 9.60
DYS391 7.43
FGA 11.55
SE33 14.49
TH01 4.45
TPOX 5.55
vWA 10.73"

stutter_percentages <- str_split(stutter_string, "\n") %>% 
  as.data.frame() %>% 
  dplyr::rename(thing = 1) %>% 
  group_by(thing) %>% 
  mutate(Marker = str_split(thing, " ")[[1]][1],
         p = as.numeric(str_split(thing, " ")[[1]][2])) %>% 
  ungroup() %>% 
  select(!thing)


# ATs from Riman et al.-s paper
AT_Riman <- "D3S1358=35 vWA=35 D16S539=35 CSF1PO=35 TPOX=35 D8S1179=65 D21S11=65 D18S51=65 D2S441=45 D19S433=45 TH01=45 FGA=45 D22S1045=50 D5S818=50 D13S317=50 D7S820=50 SE33=50 D10S1248=60 D1S1656=60 D12S391=60 D2S1338=60 AMEL=10 Yindel=10 DYS391=10"
AT_Riman_data <- str_split(AT_Riman, " ") %>% 
  as.data.frame() %>% 
  dplyr::rename(thing = 1) %>% 
  group_by(thing) %>% 
  mutate(Marker = str_split(thing, "=")[[1]][1],
         AT = as.numeric(str_split(thing, "=")[[1]][2])) %>% 
  ungroup() %>% 
  select(!thing)




################## Read in all trace files and adjust them ------------------------------

cl <- makeCluster(5)
registerDoParallel(cl)
all_markers <- read.csv(paste0(genotypes_folder, "c1.csv")) %>% #load all necessary markers to make sure nothing disappears later
  pull(Marker) #all 21 of them


loop <- foreach(trace=all_files,
                .packages = c("tidyverse")) %dopar% {
                  tracedata <- read.csv(trace)
                  trace_long <- wide_to_long(tracedata) %>% 
                    mutate(Allele = as.numeric(Allele),
                           Allele_plus1 = as.numeric(Allele) +1)
                  
                  t <- trace_long %>%
                    select(Marker, Allele, Height) %>%
                    rename("Height_plus1"=Height,
                           "Allele_plus1" = Allele) %>% 
                    left_join(stutter_percentages, by=join_by("Marker")) %>% #add all expected stutters
                    mutate(stutter_threshold = p/100 * as.numeric(Height_plus1),
                           Allele_plus1 = as.numeric(Allele_plus1))
                  
                  trace_long01 <- trace_long %>% 
                    left_join(t,
                              by=join_by("Allele_plus1", "Marker")) %>% #add +1 alleles to compare them to
                    filter(is.na(stutter_threshold) | as.numeric(Height) > stutter_threshold) %>% # stutter filter
                    
                    left_join(AT_Riman_data, by=join_by("Marker")) %>% 
                    filter(as.numeric(Height) > AT) %>%   # Apply AT
                    select(SampleName, Marker, Allele, Height) %>% 
                    arrange(Marker, as.numeric(Allele)) %>% 
                    mutate(Allele = as.character(as.numeric(Allele))) %>% 
                    mutate(SampleName = gsub(";", "-", gsub("_", "-", SampleName)))
                  
                  missing_markers <- setdiff(all_markers,
                                             unique(trace_long01$Marker))
                  
                  if(length(missing_markers)!=0){
                    trace_wide <- trace_long01 %>%
                      rbind(data.frame(
                        SampleName = unique(trace_long01$SampleName),
                        Marker = missing_markers,
                        Allele = "",
                        Height = ""
                      )) %>% 
                      long_to_wide(samplelabel = unique(trace_long01$SampleName)) %>% 
                      arrange(Marker)
                    
                    print(paste0("Missing markers: ", trace))
                  } else{
                    trace_wide <- long_to_wide(trace_long01, unique(trace_long01$SampleName)) %>% 
                      arrange(Marker)
                  }
                  
                  to_location <- paste0(str_split(trace, "/")[[1]][1], "/",
                                        str_split(trace, "/")[[1]][2],
                                        "_sFRiman")
                  write.csv(trace_wide, paste0(to_location, "/", unique(trace_wide$SampleName), ".csv"),
                            row.names=F, quote = F)
                }
stopCluster(cl)






# Read in key file
keys <- read.csv(paste0(results_folder, "logfile.csv")) %>% 
  mutate(joinby = gsub(";", "-", real)) %>% 
  rename(nomod = unmod)



# Make traces - references log file (writing up all calculations to go through in DNAStatistX)
input_folders_afterfilter <- list.files(results_folder, full.names = T)[str_detect(list.files(results_folder), "sFRiman") &
                                                          str_detect(list.files(results_folder), "real") &
                                                          !str_detect(list.files(results_folder), "batches")]
all_files_afterfilter <- lapply(input_folders_afterfilter, function(folder) list.files(folder, full.names=T)) %>% unlist()
traces <- all_files_afterfilter %>% 
  as.data.frame() %>% 
  rename("trace" = ".") %>% 
  filter(!str_detect(trace, "traces_references")) %>% 
  rowwise() %>% 
  mutate(NOC = as.numeric(str_split(str_split_i(trace, "/", 2), "p")[[1]][1]),
         
         contr1 = ifelse(str_detect(trace, "real"),
                         paste0("c",str_split(trace, "-")[[1]][4]),
                         str_split(str_split(trace, "-")[[1]][1], "/")[[1]][2]),
         contr2 = ifelse(str_detect(trace, "real"),
                         paste0("c",str_split(trace, "-")[[1]][5]),
                         str_split(trace, "-")[[1]][2]),
         
         contr3 = ifelse(NOC==2, NA,
                         ifelse(str_detect(trace, "real"),
                                paste0("c",str_split(trace, "-")[[1]][6]),
                                str_split(trace, "-")[[1]][3])),
         
         contr4 = ifelse(NOC==4,
                         ifelse(str_detect(trace, "real"),
                                paste0("c",str_split(trace, "-")[[1]][7]),
                                str_split(trace, "-")[[1]][4]),
                         NA)) %>%
  
  mutate(joinby = gsub(".csv", "", str_split(trace, "/")[[1]][3])) %>% 
  ungroup() %>% 
  full_join(keys,
            by=join_by("joinby", "NOC")) %>% 
  select(-real) %>% 
  rename(real=joinby) %>% 
  
  group_by(nomod) %>% 
  mutate(n_nomod = 1:n(),
         n_nomod_total = n()) %>% 
  ungroup()



traces2 <- traces %>% 
  pivot_longer(cols=real:mod) %>% 
  select(-trace) %>% 
  mutate(trace = paste0(NOC, "p_", name, "_sFRiman/", value, ".csv")) %>% 
  select(trace, NOC, ends_with(as.character(1:4))) %>% 
  unique() %>% 
  pivot_longer(cols = contr1:contr4) %>% 
  filter(!is.na(value)) %>% 
  select(trace, value, NOC) %>% 
  rename(reference = value) %>% 
  arrange(trace, reference) %>% 
  mutate(trace = gsub(";", "-", trace))




# Save combinations in each file
for(NOCp in unique(traces2$NOC)){
  NOCdata_nomod <- traces2 %>% 
    filter(NOC == NOCp) %>% 
    rowwise() %>%
    filter(str_detect(trace, "nomod")) %>% 
    mutate(trace = str_split_i(trace, "/", 2)) %>% 
    ungroup() %>% 
    select(trace, reference) %>% 
    arrange(trace, reference) %>% 
    unique()
  
  NOCdata_mod <- traces2 %>% 
    filter(NOC == NOCp) %>% 
    rowwise() %>%
    filter(str_detect(trace, "_mod")) %>% 
    mutate(trace = str_split_i(trace, "/", 2)) %>% 
    ungroup() %>% 
    select(trace, reference) %>% 
    arrange(trace, reference) %>% 
    unique()
  
  NOCdata_real <- traces2 %>% 
    filter(NOC == NOCp) %>% 
    rowwise() %>%
    filter(str_detect(trace, "real")) %>% 
    mutate(trace = str_split_i(trace, "/", 2)) %>% 
    ungroup() %>% 
    select(trace, reference) %>% 
    arrange(trace, reference) %>% 
    unique()
  
  # Save to mod and nomod and real
  write.xlsx(NOCdata_mod, paste0(results_folder, NOCp, "p_mod_sFRiman/traces_references.xlsx"))
  write.xlsx(NOCdata_nomod, paste0(results_folder, NOCp, "p_nomod_sFRiman/traces_references.xlsx"))
  write.xlsx(NOCdata_real, paste0(results_folder, NOCp, "p_real_sFRiman/traces_references.xlsx"))
}














