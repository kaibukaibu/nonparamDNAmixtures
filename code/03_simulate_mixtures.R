# Code for simulating mixtures in paper "Non-parametric simulation of DNA mixture profiles from one-person trace profiles"
# Author: Kai Budrikas (kaib@itu.dk)

# The code consists of 3 parts:
# 1. First we manipulate 1p traces we wish to sample from
# 2. Then we create a dataset that specifies all combinations of 1p traces + contributors we wish to create mixtures of
# 3. Then we update the peaks and paste 1p traces together (SF+AT can be applied later if you wish by adjusting the 02 code)

# In this code we shall generate 10 2-person mixtures for each combination of contributors in the "GTs_target" folder.
# Let us assume we are only interested in untreated samples with 1p profile rfus being between 10000 and 10000.
# Code is easily adjustable to also condition on target mixture proportions
# One could also adjust other parameters, fx if you want 1000 mixtures of the same set of contributors

# TO DO: smoothen the duplicate situation / add data cleaning files separately



# -----------------------------------------------------------------------------------------------

# Load packages and set data locations

library(tidyverse)
library(doParallel)


#setwd("..") #Project location
#setwd("repDNAmixtures/")
input_proved1p <- "data/data_provedit/1p_unfiltered/"  # 1p traces folder, must end with "/"
output_folder <- "data/data_simulated/" # Folder for saving the results, must end with "/"
GTs_provedit <- "data/data_provedit/genotypes/" # Folder for the genotypes of the contributors, must end with "/"
GTs_target <- "data/target_genotypes/" # Folder for the genotypes of contributors for whom we wish to create mixtures, must end with "/"
source("code/helping_functions.R") # Some functions for transforming the profile data from long to wide format and vice versa
set.seed(66)


# Set parameters
NOC <- 2
n_contrs <- length(list.files(GTs_target))
n_mixturereps <- 10 # specifying how many mixtures of the same contributor combination we want
rfu_min <- 10000
rfu_max <- 100000
treatment_type <- "Untreated" # takes values "DNase", "Humic acid", "Sonication", "Untreated", "UV"

if(NOC > n_contrs){
  stop("NOC cannot be larger than total number of target contributors!")
}


# -----------------------------------------------------------------------------------------------




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






#################### PART 1 -- selecting a subset of 1p traces ###########################################

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









# 1p traces
provedit1p00 <- list.files(input_proved1p) %>% 
  as.data.frame() %>% 
  rename("trace" = ".") %>% 
  group_by(trace) %>% 
  mutate(contr = as.character(as.numeric(substr(str_split(trace, "-")[[1]][3],1,2))),
         DNA = as.numeric(gsub("G", "", gsub("F", "", rev(str_split(trace, "-")[[1]])[2]))),
         diltreat = str_split(trace, "-")[[1]][3],
         contdil = substr(diltreat, 1, 4),
         treat = gsub(contdil, "", diltreat),
         treattype = substr(treat, 1,1),
         treattype2 = case_when(
           treattype=="a" ~ "Untreated",
           treattype %in% c("b", "c", "d", "e") ~ "DNase",
           treattype=="U" ~ "UV",
           treattype=="S" ~ "Sonication",
           treattype=="I" ~ "Humic acid"
         )) %>% 
  filter(treat!="") %>% 
  ungroup() %>% 
  mutate(rfu = NA) %>% 
  
  ungroup()


# add total rfu
for(i in 1:nrow(provedit1p00)){
  trace_i <- read.csv(paste0(input_proved1p, provedit1p00$trace[i])) %>% 
    filter(!Marker %in% c("AMEL", "Yindel", "DYS391")) %>% 
    wide_to_long2() %>% 
    
    #Fix duplicates -- some alleles have two height values, we take only the highest peak
    group_by(SampleName, Allele, Marker) %>% 
    summarise(Height = max(Height)) %>% 
    ungroup() %>% 
    
    mutate(contr = as.character(as.numeric(substr(str_split(SampleName, "-")[[1]][3],1,2)))) %>% 
    
    # Filter out everything but allelic and -1 stutter peaks after adding that info
    left_join(
      GTs01 %>% 
        mutate(contr1_allele = 1),
      by = join_by("contr"=="cont", "Marker"=="Marker", "Allele"=="Allele")
    ) %>% 
    
    mutate(
      sumsum = ifelse(is.na(contr1_allele), 0, contr1_allele)
    ) %>% 
    
    mutate(sum_perrow = sum(as.numeric(Height) * sumsum))
  
  provedit1p00[i, "rfu"] <- unique(trace_i$sum_perrow)
}


# Filter here
provedit1p_subset <- provedit1p00 %>% 
  filter(treattype2 == treatment_type &
           rfu_min <= rfu &
           rfu_max >= rfu)












######### PART 2 -- make dataset that specifies all combinations of mixtures we wish to make #################
target_contrs00 <- do.call(rbind, lapply(list.files(GTs_target, full.names = T), function(file) {
  read.csv(file) %>%
    filter(!(Marker %in% c("AMEL", "Yindel", "DYS391")))
}))

target_contrs01 <- target_contrs00 %>% 
  left_join(
    target_contrs00 %>% 
      select(SampleName) %>% 
      unique() %>% 
      mutate(contrnr = 1:n()),
    by=join_by("SampleName")
  )

combs00 <- gtools::combinations(n = n_contrs,
                                r = NOC) %>% 
  as.data.frame() %>% 
  slice(rep(1:n(), each = n_mixturereps))


combs01 <- combs00 %>% 
  
  # add mixture identifiers
  left_join(
    combs00 %>% 
      unique() %>% 
      mutate(repnr = 1:n())
  ) %>% 
  
  group_by(repnr) %>% 
  mutate(mixnr = 1:n()) %>% 
  ungroup() %>% 
  
  pivot_longer(
    cols=starts_with("V"), values_to = "contrnr"
  ) %>% 
  select(-name) %>% 
  
  left_join(
    target_contrs01 %>% 
      select(SampleName, contrnr) %>% 
      unique() %>% 
      rename(contr_target = SampleName)
  ) %>% 
  
  # Sample one 1p for each at random -- this can be adjusted if one wants other certain parameters
  # Currently it just samples two origins for each mixture without replacement (so that the two origins would not be the same in one mixture)
  # from the set of 1ps
  group_by(repnr, mixnr) %>% 
  mutate(origin_trace = sample(provedit1p_subset$trace, size = NOC, replace = F)) %>% 
  ungroup() %>% 
  
  #Add all other information
  left_join(
    provedit1p_subset, by=join_by("origin_trace"=="trace")
  ) %>% 
  
  mutate(rep_mix = paste0(repnr, "_", mixnr))















cl <- makeCluster(5)
registerDoParallel(cl)

# Pasting the mixtures together and saving
foreach(mixident = unique(combs01$rep_mix),
        .combine=rbind,
        .packages = c("tidyverse")) %dopar% {
          
          # for(mixident in unique(combs01$rep_mix)){
          
          rowdata <- combs01 %>% 
            filter(rep_mix == mixident)
          
          
          
          
          
          
          
          # # Modified 1p mixtures -------------------------------------------------------------
          
          # 
          # origin_contributors <- paste0(str_split_i(traces_1p, "-", 3), "_", gsub("GF", "", str_split_i(traces_1p, "-", 4)))
          # origin_contributors_ident <- as.character(as.numeric(substr(str_split_i(traces_1p, "-", 3), 1, 2)))
          # 
          # target_contributors <- rowdata[, c("contr1", "contr2", "contr3", "contr4")] %>% as.matrix() %>% c()
          # target_contributors <- target_contributors[!is.na(target_contributors)]
          
          settings <- data.frame(
            origin = paste0("c",rowdata$contr),  #paste0(origin_contributors),
            target = rowdata$contr_target,
            trace = rowdata$origin_trace
          )
          
          traces_1p <- settings$trace
          
          origin_contrs00 <- GTs00 %>% 
            filter(Sample.Name %in% settings$origin) %>% 
            mutate(Allele.1 = as.character(as.numeric(Allele.1)),
                   Allele.2 = as.character(as.numeric(Allele.2))
            )
          
          # Paste GTs together
          origin_contr_GTs01 <- origin_contrs00 %>%
            select(Sample.Name, Marker, Allele.1) %>%
            rename(SampleName = Sample.Name,
                   Allele = Allele.1) %>% 
            
            # Add Allele1 and 2 indicators
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
                                                         function(contr){
                                                           read.csv(paste0(input_proved1p, contr)) %>% 
                                                             filter(!(Marker %in% c("AMEL", "Yindel", "DYS391"))) %>% 
                                                             wide_to_long()
                                                         })) %>% 
            mutate(Allele = as.character(as.numeric(Allele)),
                   Height = as.numeric(Height)) %>% 
            group_by(SampleName) %>% 
            mutate(contr = paste0("c",as.character(as.numeric(substr(str_split_i(SampleName, "-", 3),1,2))))) %>% 
            ungroup() %>% 
            left_join(settings %>% 
                        mutate(trace = gsub(".csv", "", trace)) %>% 
                        select(target, trace), by=join_by("SampleName"=="trace")) %>% 
            ungroup() #%>% 
          # group_by(SampleName) %>% 
          # mutate(diltreat = str_split(SampleName, "-")[[1]][3],
          #        contdil = substr(diltreat, 1, 4),
          #        treat = gsub(contdil, "", diltreat),
          #        treattype = substr(treat, 1,1),
          #        treattype2 = case_when(
          #          treattype=="a" ~ "Untreated",
          #          treattype %in% c("b", "c", "d", "e") ~ "DNase",
          #          treattype=="U" ~ "UV",
          #          treattype=="S" ~ "Sonication",
          #          treattype=="I" ~ "Humic acid"
          #        )) %>% 
          # ungroup()
          
          
          origin_contr_traces01 <- origin_contr_traces00 %>% 
            full_join(origin_contr_GTs01 %>% 
                        #select(-target) %>% 
                        mutate(trace = gsub(".csv", "", trace)),
                      by = join_by("SampleName"=="trace", "contr"=="SampleName", "Marker"=="Marker", "Allele"=="Allele", "target"=="target")) %>% 
            arrange(contr, Marker, as.numeric(Allele))
          
          
          
          origin_contr_traces02 <- origin_contr_traces01 %>% 
            select(-SampleName) %>% 
            relocate(contr) %>% 
            mutate(Height = ifelse(is.na(Height), 0, Height))
          
          
          # treatments <- origin_contr_traces02 %>% 
          #   filter(!is.na(Allele1)) %>% 
          #   filter(Height != 0) %>% 
          #   select(treattype, contr) %>% 
          #   unique() %>% 
          #   pull(treattype)
          
          
          
          noisefromthis <- origin_contr_traces02 %>% 
            filter(!is.na(Allele1)) %>% 
            #filter(Height != 0) %>% 
            select(contr, target, Marker, Allele, Height) %>% 
            unique() %>% 
            group_by(contr, target) %>%
            summarise(heightsum = sum(Height)) %>%
            ungroup() %>%
            arrange(heightsum) %>%
            tail(1) %>%
            pull(target)
          
          
          
          noise <- origin_contr_traces02 %>% 
            filter(target==noisefromthis & is.na(Allele1)) #remove new allelic peak locations later before pasting together
          
          
          
          # Remove noise and add the 4 components which will then be added up into peaks later
          origin_contr_traces03 <- origin_contr_traces02 %>% 
            filter(!is.na(Allele1)) %>% 
            
            #Add how many unique peaks there are for figuring out which setting is homoz, hetz1dif, hetzno1dif
            group_by(contr, Marker, target) %>% 
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
            group_by(contr, Marker, target) %>% # We need to spread the value downwards for later calculations
            mutate(S1 = sum(S1)) %>% 
            ungroup() %>% 
            
            # Stutter 2
            mutate(
              S2 = case_when(
                n_unique_alleles == 2 & Stutter2 == 1 ~ Height - S1,
                n_unique_alleles == 3 & Stutter2 == 1 ~ round(p23 * Height),
                n_unique_alleles == 4 & Stutter2 == 1 ~ Height,
                TRUE ~ 0
              )
            ) %>% 
            group_by(contr, Marker, target) %>% 
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
            group_by(contr, Marker, target) %>% 
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
            group_by(contr, Marker, target) %>% 
            mutate(A2 = sum(A2)) %>% 
            ungroup() 
          
          
          # Test that they sum up correctly
          test <- origin_contr_traces03 %>% 
            select(contr, target, Marker, Allele, Height, S1, S2, A1, A2) %>% 
            unique() %>% 
            mutate(testsum1 = S1+S2+A1+A2) %>% 
            group_by(contr, target, Marker) %>% 
            mutate(testsum2 = sum(Height)) %>% 
            ungroup() %>% 
            mutate(t = testsum1-testsum2) %>% 
            filter(t != 0)
          
          if(nrow(test) != 0){
            warning("Error: selecting 4 components didn't go correctly!!!!")
            break
          }
          
          
          
          
          
          
          
          
          # Target GTs
          target_contrs05 <- target_contrs01 %>% 
            filter(SampleName %in% settings$target) %>% 
            rename(Allele.1 = Allele1,
                   Allele.2 = Allele2,
                   Sample.Name = SampleName) %>% 
            mutate(Allele.1 = as.character(as.numeric(Allele.1)),
                   Allele.2 = as.character(as.numeric(Allele.2))
            )
          
          
          # Paste GTs together
          target_contr_GTs01 <- target_contrs05 %>% 
            select(Sample.Name, Marker, Allele.1) %>%
            rename(SampleName = Sample.Name,
                   Allele = Allele.1) %>% 
            
            # Add Allele1 and 2 indicators
            mutate(Allele1 = 1,
                   Allele2 = 0,
                   Stutter1 = 0,
                   Stutter2 = 0) %>%
            rbind(
              target_contrs05 %>%
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
              target_contrs05 %>%
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
              target_contrs05 %>%
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
            select(target, Marker, S1, S2, A1, A2) %>% 
            unique() %>% 
            mutate(rfu = S1+S2+A1+A2) %>% 
            group_by(target) %>% 
            summarise(total_rfu = sum(rfu)) %>% 
            ungroup() %>% 
            left_join(
              settings %>% 
                select(target, trace)
            )
          
          ### Checking that the rfus are equal to the origin traces' ones
          # t = origin_contr_traces01 %>% 
          #   filter(!(is.na(Allele1) | is.na(SampleName))) %>% 
          #   select(SampleName, Marker, Allele, Height) %>% 
          #   unique() %>% 
          #   group_by(SampleName) %>% 
          #   summarise(Height = sum(Height)) %>% 
          #   ungroup()
          
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
                    filter(!paste0(Marker, "_", Allele) %in% paste0(data_complete01$Marker, "_", data_complete01$Allele)) %>% 
                    select(Marker, Allele, Height)) %>% 
            
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
             nrow()!=0) print(paste0("Missing peaks: ",
                                     paste0(c(paste0("rep", unique(rowdata$repnr)),
                                              settings2$target,
                                              str_split_i(settings2$trace, "-", 3),
                                              gsub("G", "", gsub("F", "", str_split_i(settings2$trace, "-", 4))),
                                              settings2$total_rfu), collapse = "-")))
          
          
          
          data_complete04 <- data_complete03 %>% 
            long_to_wide(samplelabel = paste0(c(paste0("rep", unique(rowdata$repnr)),
                                                settings2$target,
                                                str_split_i(settings2$trace, "-", 3),
                                                gsub("G", "", gsub("F", "", str_split_i(settings2$trace, "-", 4))),
                                                settings2$total_rfu), collapse = "-")) %>% 
            arrange(Marker)
          
          
          # Save results
          write.csv(data_complete04,
                    paste0(output_folder, unique(data_complete04$SampleName), ".csv"), row.names=F, quote = F)
          
          
          
          # data.frame(
          #   real = unique(real_subset2$SampleName),
          #   unmod = unique(target_nomod2$SampleName),
          #   mod = unique(data_complete04$SampleName),
          #   NOC = NOC_mix
          # )
          #print(round(row/nrow(mixtures01)*100, 2))
        }


stopCluster(cl)


write.csv(combs01,
          paste0(output_folder, "logfile.csv"), row.names=F, quote = F)
