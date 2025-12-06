# Code for generating replicate mixtures in paper "Non-parametric simulation of DNA mixture profiles from one-person trace profiles"
# Author: Kai Budrikas (kaib@itu.dk)




# -----------------------------------------------------------------------------------------------

# Load packages and set data locations

library(tidyverse)
library(doParallel)


#setwd("..") #Project location
input_provedit <- "data/data_provedit_cleaned/traces.csv" # Provedit trace profiles
output_folder <- "data/data_replicated/" # Folder for saving the results, must end with "/"
GTs_provedit <- "data/data_provedit_cleaned/genotypes/" # Folder for the genotypes of the contributors, must end with "/"
source("code/helping_functions.R") # Some functions for transforming the profile data from long to wide format and vice versa
set.seed(66)



# -----------------------------------------------------------------------------------------------




# make results subfolders
if (dir.exists(output_folder)) {
  unlink(output_folder, recursive = TRUE, force = T)  # Deletes the folder and all contents
  dir.create(output_folder)
}
for(NOC in 2:4){
  for(type in c("real", "nomod", "nomodQ", "mod")){
    dir.create(paste0(output_folder, NOC, "p_", type),
               recursive = F, showWarnings = T)
  }
}


# expected stutter
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
         p = as.numeric(str_split(thing, " ")[[1]][2]),
         p23 = p*(2/3) / 100) %>% 
  ungroup() %>% 
  select(!thing)






# Read and transform each file, then bind rows
GTs00 <- do.call(rbind, lapply(list.files(GTs_provedit, full.names = T), function(file) {
  read.csv(file) %>%
    filter(!(Marker %in% c("AMEL", "Yindel", "DYS391")))
}))


# Creating indicators for which peaks in the profile are allelic and which are -1 stutter
GTs01 <- rbind(
  # Allelic peaks
  GTs00 %>% 
    select(Sample.Name, Marker, Allele.1) %>% 
    rename(Allele=Allele.1),
  GTs00 %>% 
    select(Sample.Name, Marker, Allele.2) %>% 
    rename(Allele=Allele.2),
  
  # -1 stutters
  GTs00 %>% 
    select(Sample.Name, Marker, Allele.1) %>% 
    rename(Allele=Allele.1) %>% 
    mutate(Allele = as.character(as.numeric(Allele)-1)),
  GTs00 %>% 
    select(Sample.Name, Marker, Allele.2) %>% 
    rename(Allele=Allele.2) %>% 
    mutate(Allele = as.character(as.numeric(Allele)-1))
) %>% 
  mutate(Allele = as.character(as.numeric(Allele)),
         Sample.Name = gsub("c", "", Sample.Name)) %>% 
  rename(cont = Sample.Name) %>% 
  unique()







# Add some info for the mixtures coming from the naming convention
#all_mixes <- list.files(input_proveditmixes, full.names = T)
mixtures00 <- read.csv2(input_provedit, sep = ",") %>% 
  
  #Add NOC and take only mixtures
  mutate(NOC = str_count(SampleName, ";") + 1) %>% #produces NAs but that's fine
  filter(NOC != 1) %>% 
  
  #Add contributors
  rowwise() %>% 
  mutate(contr1 = str_split(str_split(SampleName, "-")[[1]][3], "_")[[1]][1],
         contr2 = str_split(str_split(SampleName, "-")[[1]][3], "_")[[1]][2],
         contr3 = str_split(str_split(SampleName, "-")[[1]][3], "_")[[1]][3],
         contr4 = str_split(str_split(SampleName, "-")[[1]][3], "_")[[1]][4],
         
         prop1 = as.numeric(str_split(str_split(SampleName, "-")[[1]][4], ";")[[1]][1]),
         prop2 = as.numeric(str_split(str_split(SampleName, "-")[[1]][4], ";")[[1]][2]),
         prop3 = as.numeric(str_split(str_split(SampleName, "-")[[1]][4], ";")[[1]][3]),
         prop4 = as.numeric(str_split(str_split(SampleName, "-")[[1]][4], ";")[[1]][4])) %>% 
  ungroup() %>% 
  
  group_by(SampleName) %>% 
  
  #Add treatment stuff -- this throws some warnings but it seems fine
  mutate(
    diltreat = str_split(SampleName, "-")[[1]][5],
    contdil = substr(diltreat, 1, 2),
    treat = gsub(contdil, "", diltreat),
    treattype = substr(treat, 1,1),
    treattype2 = case_when(
      treattype=="a" ~ "Untreated",
      treattype %in% c("b", "c", "d", "e") ~ "DNase",
      treattype=="U" ~ "UV",
      treattype=="S" ~ "Sonication",
      treattype=="I" ~ "Humic Acid"
    )
  ) %>% 
  
  
  # Add allelic and -1 stut peak info
  left_join(
    GTs01 %>% 
      mutate(contr1_allele = 1),
    by = join_by("contr1"=="cont", "Marker"=="Marker", "Allele"=="Allele")
  ) %>% 
  left_join(
    GTs01 %>% 
      mutate(contr2_allele = 1),
    by = join_by("contr2"=="cont", "Marker"=="Marker", "Allele"=="Allele")
  ) %>% 
  left_join(
    GTs01 %>% 
      mutate(contr3_allele = 1),
    by = join_by("contr3"=="cont", "Marker"=="Marker", "Allele"=="Allele")
  ) %>% 
  left_join(
    GTs01 %>% 
      mutate(contr4_allele = 1),
    by = join_by("contr4"=="cont", "Marker"=="Marker", "Allele"=="Allele")
  ) %>% 
  
  mutate(
    sumsum = ifelse(ifelse(is.na(contr1_allele), 0, contr1_allele) +
                      ifelse(is.na(contr2_allele), 0, contr2_allele) +
                      ifelse(is.na(contr3_allele), 0, contr3_allele) +
                      ifelse(is.na(contr4_allele), 0, contr4_allele) == 0,
                    0, 1)
  ) %>% 
  
  group_by(SampleName) %>% 
  
  #Add rfu info (sum of peak heights, which have been multiplied by the allelic/stut 1-0 indicator)
  mutate(rfu_total = sum(Height * sumsum)) %>% 
  ungroup() %>% 
  
  rowwise() %>% 
  mutate(rfu1 = rfu_total * prop1 / sum(c(prop1, prop2, prop3, prop4), na.rm = T),
         rfu2 = rfu_total * prop2 / sum(c(prop1, prop2, prop3, prop4), na.rm = T),
         rfu3 = rfu_total * prop3 / sum(c(prop1, prop2, prop3, prop4), na.rm = T),
         rfu4 = rfu_total * prop4 / sum(c(prop1, prop2, prop3, prop4), na.rm = T)) %>% 
  ungroup()



#Save rfus for later interpretation
k <- mixtures00 %>% 
  select(SampleName, rfu_total) %>% 
  unique() %>% 
  mutate(rfu_total = rfu_total,
         SampleName = gsub("\\.", "", SampleName),
         SampleName = gsub(";", "-", SampleName),
         SampleName = gsub("_", "-", SampleName),
         SampleName = paste0(SampleName, ".csv"))
write.csv(k,
          paste0(output_folder, "realmix_rfus.csv"), row.names=F, quote = F)



# For looping later
mixtures01 <- mixtures00 %>% 
  select(SampleName, treattype2, contr1, contr2, contr3, contr4,
         rfu1, rfu2, rfu3, rfu4, NOC) %>% 
  mutate(contr1_mod = NA,
         contr2_mod = NA,
         contr3_mod = NA,
         contr4_mod = NA,
         
         contr1_nomodQ = NA,
         contr2_nomodQ = NA,
         contr3_nomodQ = NA,
         contr4_nomodQ = NA,
         
         contr1_nomod = NA,
         contr2_nomod = NA,
         contr3_nomod = NA,
         contr4_nomod = NA) %>% 
  unique()





# How many mixtures do we have?
# should be 176-162-176
mixtures00 %>% 
  select(NOC, SampleName) %>%
  unique() %>% 
  group_by(NOC) %>% 
  summarise(n=n())





# 1p traces
provedit1p <- read.csv2(input_provedit, sep = ",") %>% 
  
  #Add NOC and take only 1ps
  mutate(NOC = str_count(SampleName, ";") + 1) %>% 
  filter(NOC == 1) %>%  
  rename(trace = SampleName) %>% 
  
  group_by(trace) %>% 
  mutate(contr = as.character(as.numeric(substr(str_split(trace, "-")[[1]][3],1,2))),
         DNA = as.numeric(gsub("H", "", gsub("G", "", gsub("F", "", rev(str_split(trace, "-")[[1]])[2]))))
  ) %>% 
  mutate(diltreat = ifelse(NOC==1,
                           str_split(trace, "-")[[1]][3],
                           str_split(trace, "-")[[1]][5]),
         
         treattype2 = case_when(
           substr(diltreat, nchar(diltreat), nchar(diltreat))=="a" ~ "Untreated",
           substr(diltreat, nchar(diltreat), nchar(diltreat)) %in% c("b", "c", "d", "e") ~ "DNase",
           
           str_detect(diltreat, "U") ~ "UV",
           
           str_detect(diltreat, "S") ~ "Sonication",
           
           str_detect(diltreat, "I") ~ "Humic Acid",
           
           TRUE ~ "Fragmentase"
         )) %>% 
  ungroup() %>% 
  
  # Filter out everything but allelic and -1 stutter peaks after adding that info
  left_join(
    GTs01 %>% 
      mutate(contr1_allele = 1),
    by = join_by("contr"=="cont", "Marker"=="Marker", "Allele"=="Allele")
  ) %>% 
  
  mutate(
    sumsum = ifelse(is.na(contr1_allele), 0, contr1_allele)
  ) %>% 
  
  group_by(trace) %>% 
  mutate(rfu = sum(as.numeric(Height) * sumsum)) %>% 
  ungroup()





provedit1p_noheights <- provedit1p %>% 
  select(trace, contr, NOC, treattype2, rfu) %>% 
  unique()





diff_function <- function(diff) abs(diff)

# Find most similar 1p traces for contributors
# If the loop breaks unexpectedly, then there likely isn't a 1p trace with the wished properties
for(row in 1:nrow(mixtures01)){ 
  
  ###################### Mixtures with unmodified peaks - nomod (SET1) -------------------------------
  trace1_sametreatsamecontr <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2,
           contr == mixtures01[row,]$contr1) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu1)) %>% 
    arrange(rfu_diff) %>% 
    head(1) %>% 
    pull(trace)
  
  if(length(trace1_sametreatsamecontr)==0){
    trace1_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(treattype2 == "Untreated", #if real treatment type is missing then take untreated
             contr == mixtures01[row,]$contr1) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu1)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  trace2_sametreatsamecontr <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2,
           contr == mixtures01[row,]$contr2) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu2)) %>% 
    arrange(rfu_diff) %>% 
    head(1) %>% 
    pull(trace)
  
  if(length(trace2_sametreatsamecontr)==0){
    trace2_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(treattype2 == "Untreated",
             contr == mixtures01[row,]$contr2) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu2)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  trace3_sametreatsamecontr <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2,
           contr == mixtures01[row,]$contr3) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu3)) %>% 
    arrange(rfu_diff) %>% 
    head(1) %>% 
    pull(trace)
  
  if(length(trace3_sametreatsamecontr)==0){
    trace3_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(treattype2 == "Untreated",
             contr == mixtures01[row,]$contr3) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu3)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  if(length(trace3_sametreatsamecontr)==0){
    trace3_sametreatsamecontr <- NA #if still missing then because NOC=2
  }
  
  trace4_sametreatsamecontr <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2,
           contr == mixtures01[row,]$contr4) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu4)) %>% 
    arrange(rfu_diff) %>% 
    head(1) %>% 
    pull(trace)
  
  if(length(trace4_sametreatsamecontr)==0){
    trace4_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(treattype2 == "Untreated",
             contr == mixtures01[row,]$contr4) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu4)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  if(length(trace4_sametreatsamecontr)==0){
    trace4_sametreatsamecontr <- NA #if still missing then because NOC=2 or 3
  }
  
  mixtures01[row, "contr1_nomod"] <- trace1_sametreatsamecontr
  mixtures01[row, "contr2_nomod"] <- trace2_sametreatsamecontr
  mixtures01[row, "contr3_nomod"] <- trace3_sametreatsamecontr
  mixtures01[row, "contr4_nomod"] <- trace4_sametreatsamecontr
  
  
  
  
  
  
  
  ########################## Mixes with modified peaks - mod (SET2) ----------------------------------
  trace1_sametreat <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2 & 
             contr != mixtures01[row,]$contr1) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu1)) %>% 
    arrange(rfu_diff) %>% 
    head(1) %>% 
    pull(trace)
  
  trace2_sametreat <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2,
           trace != trace1_sametreat & 
             contr != mixtures01[row,]$contr2) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu2)) %>% 
    arrange(rfu_diff) %>% 
    head(1) %>% 
    pull(trace)
  
  trace3_sametreat <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2,
           trace != trace1_sametreat,
           trace != trace2_sametreat & 
             contr != mixtures01[row,]$contr3) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu3)) %>% 
    arrange(rfu_diff) %>% 
    filter(!is.na(rfu_diff)) %>% 
    head(1) %>% 
    pull(trace)
  
  trace3_sametreat <- ifelse(length(trace3_sametreat)==0, NA, trace3_sametreat)
  
  trace4_sametreat <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2,
           trace != trace1_sametreat,
           trace != trace2_sametreat,
           trace != trace3_sametreat & 
             contr != mixtures01[row,]$contr4) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu4)) %>% 
    arrange(rfu_diff) %>% 
    filter(!is.na(rfu_diff)) %>% 
    head(1) %>% 
    pull(trace)
  
  trace4_sametreat <- ifelse(length(trace4_sametreat)==0, NA, trace4_sametreat)
  
  
  mixtures01[row, "contr1_mod"] <- trace1_sametreat
  mixtures01[row, "contr2_mod"] <- trace2_sametreat
  mixtures01[row, "contr3_mod"] <- trace3_sametreat
  mixtures01[row, "contr4_mod"] <- trace4_sametreat
  
  
  print(round(row/nrow(mixtures01)*100, 2))
  
  
  ###################### Mixtures with unmodified peaks, relaxing the treatment criteria - nomodQ (SET3) -------------------------------
  trace1_sametreatsamecontr <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2,
           contr == mixtures01[row,]$contr1) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu1)) %>% 
    arrange(rfu_diff) %>% 
    head(1) %>% 
    pull(trace)
  
  if(length(trace1_sametreatsamecontr)==0){
    trace1_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(treattype2 != "Untreated", #if real treatment type is missing then take any treated
             contr == mixtures01[row,]$contr1) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu1)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  if(length(trace1_sametreatsamecontr)==0){
    trace1_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(#treattype2 == "Untreated", #if real treatment type is missing then take any
        contr == mixtures01[row,]$contr1) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu1)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  
  
  
  
  
  
  
  trace2_sametreatsamecontr <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2,
           contr == mixtures01[row,]$contr2) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu2)) %>% 
    arrange(rfu_diff) %>% 
    head(1) %>% 
    pull(trace)
  
  if(length(trace2_sametreatsamecontr)==0){
    trace2_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(treattype2 != "Untreated",
             contr == mixtures01[row,]$contr2) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu2)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  if(length(trace2_sametreatsamecontr)==0){
    trace2_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(#treattype2 == "Untreated",
        contr == mixtures01[row,]$contr2) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu2)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  
  
  
  
  
  
  
  trace3_sametreatsamecontr <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2,
           contr == mixtures01[row,]$contr3) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu3)) %>% 
    arrange(rfu_diff) %>% 
    head(1) %>% 
    pull(trace)
  
  
  if(length(trace3_sametreatsamecontr)==0){
    trace3_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(treattype2 != "Untreated",
             contr == mixtures01[row,]$contr3) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu3)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  if(length(trace3_sametreatsamecontr)==0){
    trace3_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(#treattype2 == "Untreated",
        contr == mixtures01[row,]$contr3) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu3)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  if(length(trace3_sametreatsamecontr)==0){
    trace3_sametreatsamecontr <- NA #if still missing then because NOC=2
  }
  
  trace4_sametreatsamecontr <- provedit1p_noheights %>% 
    filter(treattype2 == mixtures01[row,]$treattype2,
           contr == mixtures01[row,]$contr4) %>% 
    mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu4)) %>% 
    arrange(rfu_diff) %>% 
    head(1) %>% 
    pull(trace)
  
  if(length(trace4_sametreatsamecontr)==0){
    trace4_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(treattype2 != "Untreated",
             contr == mixtures01[row,]$contr4) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu4)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  if(length(trace4_sametreatsamecontr)==0){
    trace4_sametreatsamecontr <- provedit1p_noheights %>% 
      filter(#treattype2 == "Untreated",
        contr == mixtures01[row,]$contr4) %>% 
      mutate(rfu_diff = diff_function(rfu-mixtures01[row,]$rfu4)) %>% 
      arrange(rfu_diff) %>% 
      head(1) %>% 
      pull(trace)
  }
  
  if(length(trace4_sametreatsamecontr)==0){
    trace4_sametreatsamecontr <- NA #if still missing then because NOC=2 or 3
  }
  
  mixtures01[row, "contr1_nomodQ"] <- trace1_sametreatsamecontr
  mixtures01[row, "contr2_nomodQ"] <- trace2_sametreatsamecontr
  mixtures01[row, "contr3_nomodQ"] <- trace3_sametreatsamecontr
  mixtures01[row, "contr4_nomodQ"] <- trace4_sametreatsamecontr
}











cl <- makeCluster(5)
registerDoParallel(cl)

# Pasting the mixtures together and saving
log_file <- foreach(row = 1:nrow(mixtures01),
                    .combine=rbind,
                    .packages = c("tidyverse")) %dopar% {
                      # for(i in 1:nrow(mixtures01)){
                      rowdata <- mixtures01[row,]
                      
                      NOC_mix <- rowdata$NOC
                      
                      contributors <- rowdata[, c("contr1", "contr2", "contr3", "contr4")] %>% as.matrix() %>% c()
                      contributors <- contributors[!is.na(contributors)]
                      
                      
                      
                      # Save real data as single files
                      real_subset <- mixtures00 %>% 
                        filter(SampleName == rowdata$SampleName) %>% 
                        select(SampleName, Marker, Allele, Height) %>% 
                        mutate(Allele = as.character(as.numeric(Allele)))
                      
                      missing_markers <- setdiff(unique(real_subset$Marker),
                                                 unique(mixtures00$Marker))
                      
                      if(length(missing_markers)!=0){
                        real_subset2 <- real_subset %>%
                          rbind(data.frame(
                            Marker = missing_markers,
                            Allele = "",
                            Height = ""
                          )) %>% 
                          long_to_wide(samplelabel = rowdata$SampleName) %>% 
                          arrange(Marker)
                      } else{
                        real_subset2 <- real_subset %>%
                          long_to_wide(samplelabel = rowdata$SampleName) %>% 
                          arrange(Marker)
                      }
                      
                      # Save results
                      write.csv(real_subset2,
                                paste0(output_folder, NOC_mix, "p_real/", unique(real_subset2$SampleName), ".csv"), row.names=F, quote = F)
                      
                      
                      
                      
                      
                      ###################### No modifications first (SET1) ------------------------------------------------------------
                      traces_1p <- rowdata[, c("contr1_nomod", "contr2_nomod", "contr3_nomod", "contr4_nomod")] %>% as.matrix() %>% c()
                      traces_1p <- traces_1p[!is.na(traces_1p)]
                      
                      target_contrs00 <- do.call(rbind, lapply(contributors,
                                                               function(contr){
                                                                 read.csv(paste0(GTs_provedit, "c", contr, ".csv"))
                                                               }))
                      
                      # Paste GTs together
                      target_contr_GTs01 <- target_contrs00 %>% 
                        select(Sample.Name, Marker, Allele.1) %>% 
                        rename(SampleName = Sample.Name,
                               Allele = Allele.1) %>% 
                        rbind(
                          target_contrs00 %>% 
                            select(Sample.Name, Marker, Allele.2) %>% 
                            rename(SampleName = Sample.Name,
                                   Allele = Allele.2)
                        ) %>% 
                        mutate(Allele = as.character(as.numeric(Allele))) %>% 
                        unique()
                      
                      
                      
                      
                      
                      # Paste 1p traces into one data set
                      target_contr_traces00 <- do.call(rbind, lapply(traces_1p,
                                                                     function(con){
                                                                       provedit1p %>% 
                                                                         filter(trace==con)
                                                                     })) %>% 
                        mutate(Allele = as.character(as.numeric(Allele)),
                               Height = as.numeric(Height),
                               contr = paste0("c", contr)) %>% 
                        rename(SampleName = trace)
                      
                      
                      # Add indicators for what is allelic and what is stutter
                      target_contr_traces01 <- target_contr_traces00 %>% 
                        left_join(target_contr_GTs01 %>% 
                                    mutate(amallelic = 1),
                                  by=join_by("contr"=="SampleName", "Marker"=="Marker", "Allele"=="Allele")) %>% 
                        left_join(target_contr_GTs01 %>% 
                                    mutate(Allele = as.character(as.numeric(Allele)-1),
                                           amstutter = 1),
                                  by=join_by("contr"=="SampleName", "Marker"=="Marker", "Allele"=="Allele")) %>% 
                        mutate(amallelic = ifelse(is.na(amallelic), 0, amallelic),
                               amstutter = ifelse(is.na(amstutter), 0, amstutter),
                               
                               Height2 = Height * ifelse(amallelic==1 | amstutter==1, 1, 0)) %>% 
                        
                        group_by(contr) %>% 
                        mutate(total_rfu = sum(Height2)) %>% 
                        ungroup()
                      
                      
                      
                      
                      
                      treatments <- target_contr_traces01 %>% 
                        select(treattype2, contr) %>% 
                        unique() %>% 
                        pull(treattype2)
                      
                      
                      # If there are only untreated traces then choose noise from largest untreated contributor
                      if(length(treatments)==sum(str_detect(treatments, "Untreated"))){
                        
                        noisefromthis <- target_contr_traces01 %>% 
                          select(contr, total_rfu) %>% 
                          unique() %>% 
                          arrange(total_rfu) %>% 
                          tail(1) %>% 
                          pull(contr)
                        
                      } else{ #If there are other treatments besides untreated then take noise from largest treated contributor
                        
                        noisefromthis <- target_contr_traces01 %>% 
                          filter(treattype2!="Untreated") %>% 
                          select(contr, total_rfu) %>% 
                          unique() %>% 
                          arrange(total_rfu) %>% 
                          tail(1) %>% 
                          pull(contr)
                      }
                      
                      
                      
                      
                      
                      ending <- target_contr_traces01 %>% 
                        group_by(SampleName) %>% 
                        mutate(origin = str_split(SampleName, "-")[[1]][3],
                               temp = sub("F","",sub("H","",sub("G","",str_split(SampleName, "-")[[1]][4])))) %>% 
                        ungroup() %>% 
                        select(contr, total_rfu, origin, temp) %>% 
                        unique() %>% 
                        arrange(total_rfu)
                      
                      rfus <- ending$total_rfu
                      contrs <- ending$contr
                      origins <- ending$origin
                      templates <- ending$temp
                      
                      # Choose noise only from biggest contributors and add all peaks up
                      target_nomod <- target_contr_traces01 %>% 
                        filter(contr == noisefromthis | (amallelic==1  | amstutter==1)) %>% 
                        group_by(Marker, Allele) %>% 
                        summarise(Height = sum(Height)) %>% 
                        ungroup() 
                      
                      missing_markers <- setdiff(unique(target_contr_GTs01$Marker),
                                                 unique(target_nomod$Marker))
                      
                      if(length(missing_markers)!=0){
                        target_nomod2 <- target_nomod %>%
                          rbind(data.frame(
                            Marker = missing_markers,
                            Allele = "",
                            Height = ""
                          )) %>% 
                          long_to_wide(samplelabel = paste0(c(contrs, origins, templates, rfus
                          ), collapse = "-")) %>% 
                          arrange(Marker)
                        print("Missing markers")
                      } else{
                        target_nomod2 <- target_nomod %>%
                          long_to_wide(samplelabel = paste0(c(contrs, origins, templates, rfus
                          ), collapse = "-")) %>% 
                          arrange(Marker)
                      }
                      
                      
                      
                      # Save results
                      write.csv(target_nomod2,
                                paste0(output_folder, NOC_mix, "p_nomod/", unique(target_nomod2$SampleName), ".csv"), row.names=F, quote = F)
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      ###################### Now no modifications but using all treatments (SET3)  ------------------------------------------------------------
                      traces_1p <- rowdata[, c("contr1_nomodQ", "contr2_nomodQ", "contr3_nomodQ", "contr4_nomodQ")] %>% as.matrix() %>% c()
                      traces_1p <- traces_1p[!is.na(traces_1p)]
                      
                      target_contrs00 <- do.call(rbind, lapply(contributors,
                                                               function(contr){
                                                                 read.csv(paste0(GTs_provedit, "c", contr, ".csv"))
                                                               }))
                      
                      # Paste GTs together
                      target_contr_GTs01 <- target_contrs00 %>% 
                        select(Sample.Name, Marker, Allele.1) %>% 
                        rename(SampleName = Sample.Name,
                               Allele = Allele.1) %>% 
                        rbind(
                          target_contrs00 %>% 
                            select(Sample.Name, Marker, Allele.2) %>% 
                            rename(SampleName = Sample.Name,
                                   Allele = Allele.2)
                        ) %>% 
                        mutate(Allele = as.character(as.numeric(Allele))) %>% 
                        unique()
                      
                      
                      
                      
                      
                      # Paste 1p traces into one data set
                      target_contr_traces00 <- do.call(rbind, lapply(traces_1p,
                                                                     function(con){
                                                                       provedit1p %>% 
                                                                         filter(trace==con)
                                                                     })) %>% 
                        mutate(Allele = as.character(as.numeric(Allele)),
                               Height = as.numeric(Height),
                               contr = paste0("c", contr)) %>% 
                        rename(SampleName = trace)
                      
                      
                      # Add indicators for what is allelic and what is stutter
                      target_contr_traces01 <- target_contr_traces00 %>% 
                        left_join(target_contr_GTs01 %>% 
                                    mutate(amallelic = 1),
                                  by=join_by("contr"=="SampleName", "Marker"=="Marker", "Allele"=="Allele")) %>% 
                        left_join(target_contr_GTs01 %>% 
                                    mutate(Allele = as.character(as.numeric(Allele)-1),
                                           amstutter = 1),
                                  by=join_by("contr"=="SampleName", "Marker"=="Marker", "Allele"=="Allele")) %>% 
                        mutate(amallelic = ifelse(is.na(amallelic), 0, amallelic),
                               amstutter = ifelse(is.na(amstutter), 0, amstutter),
                               
                               Height2 = Height * ifelse(amallelic==1 | amstutter==1, 1, 0)) %>% 
                        
                        group_by(contr) %>% 
                        mutate(total_rfu = sum(Height2)) %>% 
                        ungroup()
                      
                      
                      
                      
                      
                      treatments <- target_contr_traces01 %>% 
                        select(treattype2, contr) %>% 
                        unique() %>% 
                        pull(treattype2)
                      
                      
                      # If there are only untreated traces then choose noise from largest untreated contributor
                      if(length(treatments)==sum(str_detect(treatments, "Untreated"))){
                        
                        noisefromthis <- target_contr_traces01 %>% 
                          select(contr, total_rfu) %>% 
                          unique() %>% 
                          arrange(total_rfu) %>% 
                          tail(1) %>% 
                          pull(contr)
                        
                      } else{ #If there are other treatments besides untreated then take noise from largest treated contributor
                        
                        noisefromthis <- target_contr_traces01 %>% 
                          filter(treattype2!="Untreated") %>% 
                          select(contr, total_rfu) %>% 
                          unique() %>% 
                          arrange(total_rfu) %>% 
                          tail(1) %>% 
                          pull(contr)
                      }
                      
                      
                      
                      
                      
                      ending <- target_contr_traces01 %>% 
                        group_by(SampleName) %>% 
                        mutate(origin = str_split(SampleName, "-")[[1]][3],
                               temp = sub("F","",sub("H","",sub("G","",str_split(SampleName, "-")[[1]][4])))) %>% 
                        ungroup() %>% 
                        select(contr, total_rfu, origin, temp) %>% 
                        unique() %>% 
                        arrange(total_rfu)
                      
                      rfus <- ending$total_rfu
                      contrs <- ending$contr
                      origins <- ending$origin
                      templates <- ending$temp
                      
                      # Choose noise only from biggest contributors and add all peaks up
                      target_nomod <- target_contr_traces01 %>% 
                        filter(contr == noisefromthis | (amallelic==1  | amstutter==1)) %>% 
                        group_by(Marker, Allele) %>% 
                        summarise(Height = sum(Height)) %>% 
                        ungroup() 
                      
                      missing_markers <- setdiff(unique(target_contr_GTs01$Marker),
                                                 unique(target_nomod$Marker))
                      
                      if(length(missing_markers)!=0){
                        target_nomod2Q <- target_nomod %>%
                          rbind(data.frame(
                            Marker = missing_markers,
                            Allele = "",
                            Height = ""
                          )) %>% 
                          long_to_wide(samplelabel = paste0(c(contrs, origins, templates, rfus
                          ), collapse = "-")) %>% 
                          arrange(Marker)
                        print("Missing markers")
                      } else{
                        target_nomod2Q <- target_nomod %>%
                          long_to_wide(samplelabel = paste0(c(contrs, origins, templates, rfus
                          ), collapse = "-")) %>% 
                          arrange(Marker)
                      }
                      
                      
                      
                      # Save results
                      write.csv(target_nomod2Q,
                                paste0(output_folder, NOC_mix, "p_nomodQ/", unique(target_nomod2Q$SampleName), ".csv"), row.names=F, quote = F)
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      ###################### Modified 1p mixtures -------------------------------------------------------------
                      traces_1p <- rowdata[, c("contr1_mod", "contr2_mod", "contr3_mod", "contr4_mod")] %>% as.matrix() %>% c()
                      traces_1p <- traces_1p[!is.na(traces_1p)]
                      
                      origin_contributors <- paste0(str_split_i(traces_1p, "-", 3), "_", gsub("GF", "", str_split_i(traces_1p, "-", 4)))
                      origin_contributors_ident <- as.character(as.numeric(substr(str_split_i(traces_1p, "-", 3), 1, 2)))
                      
                      target_contributors <- rowdata[, c("contr1", "contr2", "contr3", "contr4")] %>% as.matrix() %>% c()
                      target_contributors <- target_contributors[!is.na(target_contributors)]
                      
                      settings <- data.frame(
                        origin = origin_contributors,
                        target = paste0("c", target_contributors),
                        trace = traces_1p
                      )
                      
                      origin_contrs00 <- do.call(rbind, lapply(settings$origin,
                                                               function(contr){
                                                                 read.csv(paste0(GTs_provedit, "c",
                                                                                 as.character(as.numeric(substr(contr, 1, 2))),
                                                                                 ".csv")) %>% 
                                                                   mutate(Sample.Name = contr)
                                                               })) %>% 
                        mutate(Allele.1 = as.character(as.numeric(Allele.1)),
                               Allele.2 = as.character(as.numeric(Allele.2))
                        )
                      
                      # Paste GTs together
                      origin_contr_GTs01 <- origin_contrs00 %>%
                        select(Sample.Name, Marker, Allele.1) %>%
                        rename(SampleName = Sample.Name,
                               Allele = Allele.1) %>% 
                        
                        # Add Allele1 and 2 indicatos
                        mutate(Allele1 = 1,
                               Allele2 = 0,
                               Stutter1 = 0,
                               Stutter2 = 0) %>%
                        rbind(
                          origin_contrs00 %>%
                            select(Sample.Name, Marker, Allele.2) %>%
                            rename(SampleName = Sample.Name,
                                   Allele = Allele.2) %>% 
                            mutate(Allele1 = 0,
                                   Allele2 = 1,
                                   Stutter1 = 0,
                                   Stutter2 = 0)
                        ) %>%
                        # Add -1 stutters
                        rbind(
                          origin_contrs00 %>%
                            select(Sample.Name, Marker, Allele.1) %>%
                            rename(SampleName = Sample.Name,
                                   Allele = Allele.1) %>% 
                            mutate(Allele = as.character(as.numeric(Allele)-1),
                                   Allele1 = 0,
                                   Allele2 = 0,
                                   Stutter1 = 1,
                                   Stutter2 = 0)
                        ) %>% 
                        rbind(
                          origin_contrs00 %>%
                            select(Sample.Name, Marker, Allele.2) %>%
                            rename(SampleName = Sample.Name,
                                   Allele = Allele.2) %>% 
                            mutate(Allele = as.character(as.numeric(Allele)-1),
                                   Allele1 = 0,
                                   Allele2 = 0,
                                   Stutter1 = 0,
                                   Stutter2 = 1)
                        ) %>% 
                        arrange(SampleName, Marker, as.numeric(Allele)) %>% 
                        left_join(settings,
                                  by=join_by("SampleName"=="origin"))
                      
                      
                      # Paste 1p traces into one data set
                      origin_contr_traces00 <- do.call(rbind, lapply(traces_1p,
                                                                     function(con){
                                                                       provedit1p %>% 
                                                                         filter(trace==con)
                                                                     })) %>% 
                        mutate(Allele = as.character(as.numeric(Allele)),
                               Height = as.numeric(Height)) %>% 
                        rename(SampleName = trace) %>% 
                        mutate(contr = paste0(str_split_i(SampleName, "-", 3), "_", gsub("GF", "", str_split_i(SampleName, "-", 4))))
                      
                      
                      origin_contr_traces01 <- origin_contr_traces00 %>% 
                        full_join(origin_contr_GTs01,
                                  by = join_by("contr"=="SampleName", "Marker"=="Marker", "Allele"=="Allele")) %>% 
                        arrange(contr, Marker, as.numeric(Allele))
                      
                      
                      
                      origin_contr_traces02 <- origin_contr_traces01 %>% 
                        select(-SampleName) %>% 
                        relocate(contr) %>% 
                        mutate(Height = ifelse(is.na(Height), 0, Height))
                      
                      
                      treatments <- origin_contr_traces02 %>% 
                        filter(!is.na(Allele1)) %>% 
                        filter(Height != 0) %>% 
                        select(treattype2, contr) %>% 
                        unique() %>% 
                        pull(treattype2)
                      
                      
                      #Take noise from largest contributor
                      noisefromthis <- origin_contr_traces02 %>% 
                        filter(!is.na(Allele1)) %>% 
                        filter(Height != 0) %>% 
                        select(contr, Marker, Allele, Height) %>% 
                        unique() %>% 
                        group_by(contr) %>%
                        summarise(heightsum = sum(Height)) %>%
                        ungroup() %>%
                        arrange(heightsum) %>%
                        tail(1) %>%
                        pull(contr)
                      
                      
                      
                      noise <- origin_contr_traces02 %>% 
                        filter(contr==noisefromthis & is.na(Allele1)) #remove new allelic peak locations later before pasting together
                      
                      
                      
                      # Remove noise and add the 4 components which will then be added up into peaks later
                      origin_contr_traces03 <- origin_contr_traces02 %>% 
                        
                        # Add +1 height info for (a, b=a+1) hetzyg
                        left_join(origin_contr_traces00 %>% 
                                    mutate(Allele_minus1 = as.character(as.numeric(Allele)-1)) %>% 
                                    select(contr, Allele_minus1, Height, Marker) %>% 
                                    rename(Height_plus1 = Height),
                                  by=join_by("contr"=="contr", "Marker"=="Marker", "Allele"=="Allele_minus1")) %>% 
                        mutate(Height_plus1 = ifelse(is.na(Height_plus1), 0, Height_plus1)) %>% 
                        
                        filter(!is.na(Allele1)) %>% 
                        
                        #Add how many unique peaks there are for figuring out which setting is homoz, hetz1dif, hetzno1dif
                        group_by(contr, Marker) %>% 
                        mutate(n_unique_alleles = n_distinct(Allele)) %>% 
                        
                        #Add random components and stutter stuff
                        mutate(
                          q = runif(1, min=0.4, max=0.6)
                        ) %>% 
                        ungroup() %>% 
                        left_join(stutter_percentages %>% 
                                    select(Marker, p23),
                                  by = join_by("Marker"=="Marker")) %>% 
                        
                        #Add the 4 components
                        # Stutter 1
                        mutate(
                          S1 = case_when(
                            n_unique_alleles == 2 & Stutter1 == 1 ~ round(q * Height),
                            n_unique_alleles == 3 & Stutter1 == 1 ~ Height,
                            n_unique_alleles == 4 & Stutter1 == 1 ~ Height,
                            TRUE ~ 0
                          )
                        ) %>% 
                        group_by(contr, Marker) %>% # We need to spread the value downwards for later calculations
                        mutate(S1 = sum(S1)) %>% 
                        ungroup() %>% 
                        
                        # Stutter 2
                        mutate(
                          S2 = case_when(
                            n_unique_alleles == 2 & Stutter2 == 1 ~ Height - S1,
                            n_unique_alleles == 3 & Stutter2 == 1 ~ case_when(Height == 0 ~ 0,
                                                                              Height != 0 & Height < round(p23 * Height_plus1) ~ Height, #If expected S2 is larger than the actual A1 peak, we take the A1 height itself
                                                                              TRUE ~ round(p23 * Height_plus1)),
                            n_unique_alleles == 4 & Stutter2 == 1 ~ Height,
                            TRUE ~ 0
                          )
                        ) %>% 
                        group_by(contr, Marker) %>% 
                        mutate(S2 = sum(S2)) %>% 
                        ungroup() %>% 
                        
                        # Allele 1
                        mutate(
                          A1 = case_when(
                            n_unique_alleles == 2 & Allele1 == 1 ~ round(q * Height),
                            n_unique_alleles == 3 & Allele1 == 1 ~ Height - S2,
                            n_unique_alleles == 4 & Allele1 == 1 ~ Height,
                            TRUE ~ 0
                          )
                        ) %>% 
                        group_by(contr, Marker) %>% 
                        mutate(A1 = sum(A1)) %>% 
                        ungroup() %>% 
                        
                        # Allele 2
                        mutate(
                          A2 = case_when(
                            n_unique_alleles == 2 & Allele2 == 1 ~ Height - A1,
                            n_unique_alleles == 3 & Allele2 == 1 ~ Height,
                            n_unique_alleles == 4 & Allele2 == 1 ~ Height,
                            TRUE ~ 0
                          )
                        ) %>% 
                        group_by(contr, Marker) %>% 
                        mutate(A2 = sum(A2)) %>% 
                        ungroup() 
                      
                      
                      # Test that they sum up correctly
                      test <- origin_contr_traces03 %>% 
                        select(contr, Marker, Allele, Height, S1, S2, A1, A2) %>% 
                        unique() %>% 
                        mutate(testsum1 = S1+S2+A1+A2) %>% 
                        group_by(contr, Marker) %>% 
                        mutate(testsum2 = sum(Height)) %>% 
                        ungroup() %>% 
                        mutate(t = testsum1-testsum2) %>% 
                        filter(t != 0)
                      
                      if(nrow(test) != 0){
                        warning(paste0("Error: selecting 4 components didn't go correctly! Rownr: ", row))
                        break
                      }
                      
                      if(origin_contr_traces03 %>%
                         filter(A1 < 0) %>%
                         nrow() != 0){
                        warning(paste0("Error: A1/S2 split didn't go correctly! Rownr: ", row))
                        break
                      }
                      
                      
                      
                      
                      
                      
                      
                      
                      # Target GTs
                      target_contrs00 <- do.call(rbind, lapply(contributors,
                                                               function(contr){
                                                                 read.csv(paste0(GTs_provedit, "c", contr, ".csv"))
                                                               })) %>% 
                        mutate(Allele.1 = as.character(as.numeric(Allele.1)),
                               Allele.2 = as.character(as.numeric(Allele.2))
                        )
                      
                      
                      # Paste GTs together
                      target_contr_GTs01 <- target_contrs00 %>% 
                        select(Sample.Name, Marker, Allele.1) %>%
                        rename(SampleName = Sample.Name,
                               Allele = Allele.1) %>% 
                        
                        # Add Allele1 and 2 indicators
                        mutate(Allele1 = 1,
                               Allele2 = 0,
                               Stutter1 = 0,
                               Stutter2 = 0) %>%
                        rbind(
                          target_contrs00 %>%
                            select(Sample.Name, Marker, Allele.2) %>%
                            rename(SampleName = Sample.Name,
                                   Allele = Allele.2) %>% 
                            mutate(Allele1 = 0,
                                   Allele2 = 1,
                                   Stutter1 = 0,
                                   Stutter2 = 0)
                        ) %>%
                        # Add -1 stutters
                        rbind(
                          target_contrs00 %>%
                            select(Sample.Name, Marker, Allele.1) %>%
                            rename(SampleName = Sample.Name,
                                   Allele = Allele.1) %>% 
                            mutate(Allele = as.character(as.numeric(Allele)-1),
                                   Allele1 = 0,
                                   Allele2 = 0,
                                   Stutter1 = 1,
                                   Stutter2 = 0)
                        ) %>% 
                        rbind(
                          target_contrs00 %>%
                            select(Sample.Name, Marker, Allele.2) %>%
                            rename(SampleName = Sample.Name,
                                   Allele = Allele.2) %>% 
                            mutate(Allele = as.character(as.numeric(Allele)-1),
                                   Allele1 = 0,
                                   Allele2 = 0,
                                   Stutter1 = 0,
                                   Stutter2 = 1)
                        ) %>% 
                        arrange(SampleName, Marker, as.numeric(Allele))
                      
                      
                      
                      
                      # Add two sets of genotypes and peak heights together
                      data_complete00 <- origin_contr_traces03 %>% 
                        full_join(target_contr_GTs01 %>% 
                                    rename(target_Allele = Allele),
                                  by=join_by("target"=="SampleName", "Marker"=="Marker",
                                             "Allele1"=="Allele1",
                                             "Allele2"=="Allele2",
                                             "Stutter1"=="Stutter1",
                                             "Stutter2"=="Stutter2")) %>% 
                        
                        # Turn some to 0 to sum over later
                        mutate(S1_new = ifelse(Stutter1==1, S1, 0),
                               S2_new = ifelse(Stutter2==1, S2, 0),
                               A1_new = ifelse(Allele1==1, A1, 0),
                               A2_new = ifelse(Allele2==1, A2, 0))
                      
                      settings2 <- data_complete00 %>% 
                        select(target, trace, Marker, S1, S2, A1, A2) %>% 
                        unique() %>% 
                        mutate(rfu = S1+S2+A1+A2) %>% 
                        group_by(target, trace) %>% 
                        summarise(total_rfu = sum(rfu)) %>% 
                        ungroup()
                      
                      
                      data_complete01 <- data_complete00 %>% 
                        select(target, Marker, target_Allele, S1_new, S2_new, A1_new, A2_new) %>% 
                        
                        #Pivot to longer format
                        pivot_longer(cols=S1_new:A2_new) %>% 
                        
                        # Add up the peaks
                        group_by(Marker, target_Allele) %>% 
                        summarise(Height = sum(value)) %>% 
                        ungroup() %>% 
                        rename(Allele = target_Allele)
                      
                      
                      # Add noise
                      data_complete02 <- data_complete01 %>% 
                        rbind(noise %>% 
                                #filter(!paste0(Marker, "_", Allele) %in% paste0(data_complete01$Marker, "_", data_complete01$Allele)) %>% 
                                select(Marker, Allele, Height)) %>% 
                        
                        #Add up the peaks
                        group_by(Marker, Allele) %>% 
                        summarise(Height = sum(Height)) %>% 
                        ungroup() %>% 
                        
                        # Check which markers dont have any peaks, we need to keep them in our dataset even if they lose the peaks
                        group_by(Marker) %>% 
                        mutate(nopeaks = ifelse(sum(Height)==0, 1, 0)) %>% 
                        ungroup()
                      
                      data_complete03 <- data_complete02 %>% 
                        
                        #First we take markers which have at least one peak
                        filter(nopeaks==0) %>% 
                        select(-nopeaks) %>% 
                        filter(Height != 0) %>% 
                        
                        #And now we add the markers that have missing peaks
                        rbind(data_complete02 %>% 
                                filter(nopeaks==1) %>% 
                                mutate(Allele = "",
                                       Height = "") %>% 
                                select(-nopeaks) %>% 
                                unique())
                      
                      if(data_complete02 %>% 
                         filter(nopeaks==1) %>% 
                         select(-nopeaks) %>% 
                         unique() %>% 
                         nrow()!=0) print("Missing peaks")
                      
                      
                      
                      data_complete04 <- data_complete03 %>% 
                        long_to_wide(samplelabel = paste0(c(settings2$target,
                                                            str_split_i(settings2$trace, "-", 3),
                                                            gsub("GF", "", str_split_i(settings2$trace, "-", 4)),
                                                            settings2$total_rfu
                        ), collapse = "-")) %>% 
                        arrange(Marker)
                      
                      
                      # Save results
                      write.csv(data_complete04,
                                paste0(output_folder, NOC_mix, "p_mod/", unique(data_complete04$SampleName), ".csv"), row.names=F, quote = F)
                      
                      
                      
                      data.frame(
                        real = unique(real_subset2$SampleName),
                        nomod = unique(target_nomod2$SampleName),
                        nomodQ = unique(target_nomod2Q$SampleName),
                        mod = unique(data_complete04$SampleName),
                        NOC = NOC_mix
                      )
                      #print(round(row/nrow(mixtures01)*100, 2))
                    }


stopCluster(cl)


write.csv(log_file,
          paste0(output_folder, "logfile.csv"), row.names=F, quote = F)

openxlsx::write.xlsx(log_file,
                     paste0(output_folder, "logfile.xlsx"))






