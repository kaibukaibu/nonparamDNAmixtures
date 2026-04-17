# Code for generating replicate mixtures with EFM in paper "Non-parametric simulation of DNA mixture profiles from one-person trace profiles"
# Author: Kai Budrikas (kaib@itu.dk)




# -----------------------------------------------------------------------------------------------

# Load packages and set data locations

#setwd("..")
library(tidyverse)
#remotes::install_github("oyvble/euroformix")
library(euroformix)
library(rlist)
source("code/helping_functions.R")
set.seed(66)



input_parameters_fileloc <- "LR_results/replicates_unfiltered001riman.csv"
output <- "data/data_replicated/"
pop <- "data/NIST1036_Cauc.csv"
input_1p_fileloc <- "data/data_provedit_cleaned/genotypes.csv"




# -----------------------------------------------------------------------------------------------




# make results subfolders (and delete old results)
for(NOC in 2:4){
  for(type in c("EFM")){
    if(dir.exists(paste0(output, NOC, "p_", type))) unlink(paste0(output, NOC, "p_", type), recursive = TRUE, force = T)
    dir.create(paste0(output, NOC, "p_", type),
               recursive = F, showWarnings = T)
  }
}



# Read in LR calculation results
input_parameters <- read.csv2(input_parameters_fileloc) %>%  #should be of size 176+176+162=514
  filter(label=="real") %>% 
  select(trace, NoC, starts_with("Hu")) %>% 
  group_by(trace) %>% 
  filter(Hu_expectedPeakHeight == sample(Hu_expectedPeakHeight,1))



# mixture contributors
input_1p <- read.csv(input_1p_fileloc) %>%
  pivot_longer(cols=c("Allele.1", "Allele.2")) %>%
  rename(Allele = value) %>%
  select(!name) %>%
  mutate(Sample.Name = gsub("c", "", Sample.Name),
         Allele = as.character(as.numeric(Allele)),
         Marker = toupper(Marker)) #capitalize all the markers because otherwise EFM won't work



# dataset of popfreqs
popfreqs <- read.csv(pop) %>% 
  pivot_longer(cols=2:25) %>%
  rename(freq = value,
         Marker = name) %>%
  filter(!is.na(freq)) %>% 
  filter(toupper(Marker) %in% unique(input_1p$Marker)) %>% 
  mutate(Marker = toupper(Marker))

popfreqs2 <- lapply(split(popfreqs, popfreqs$Marker), function(x) {
  freqs <- x$freq
  names(freqs) <- x$Allele
  freqs })






for(i in 1:nrow(input_parameters)){
  tracerow <- input_parameters[i,]
  
  #Get contribution props from naming convention
  prop1 <- case_when(
    tracerow$NoC==2 ~ str_split_i(tracerow$trace, "-", 6),
    tracerow$NoC==3 ~ str_split_i(tracerow$trace, "-", 7),
    tracerow$NoC==4 ~ str_split_i(tracerow$trace, "-", 8)
  )
  prop2 <- case_when(
    tracerow$NoC==2 ~ str_split_i(tracerow$trace, "-", 7),
    tracerow$NoC==3 ~ str_split_i(tracerow$trace, "-", 8),
    tracerow$NoC==4 ~ str_split_i(tracerow$trace, "-", 9)
  )
  prop3 <- case_when(
    tracerow$NoC==2 ~ NA,
    tracerow$NoC==3 ~ str_split_i(tracerow$trace, "-", 9),
    tracerow$NoC==4 ~ str_split_i(tracerow$trace, "-", 10)
  )
  prop4 <- case_when(
    tracerow$NoC==2 ~ NA,
    tracerow$NoC==3 ~ NA,
    tracerow$NoC==4 ~ str_split_i(tracerow$trace, "-", 11)
  )
  
  props00 <- as.numeric(na.exclude(c(prop1, prop2, prop3, prop4)))
  props <- props00 / sum(props00)
  
  
  # Get contributors from naming convention
  contr1 <- str_split_i(tracerow$trace, "-", 4)
  contr2 <- str_split_i(tracerow$trace, "-", 5)
  contr3 <- case_when(
    tracerow$NoC==2 ~ NA,
    tracerow$NoC %in% c(3,4) ~ str_split_i(tracerow$trace, "-", 6)
  )
  contr4 <- case_when(
    tracerow$NoC %in% c(2,3) ~ NA,
    tracerow$NoC==4 ~ str_split_i(tracerow$trace, "-", 7)
  )
  
  contrs <- as.character(na.exclude(c(contr1, contr2, contr3, contr4)))
  
  # Loop over contrs, read in their genotypes and save them in a format necessary for EFM
  locs <- unique(input_1p$Marker)
  refdat <- setNames(vector("list", length(locs)), locs)
  
  for (j in seq_along(contrs)) {
    cont <- contrs[j]
    sub <- input_1p %>% filter(Sample.Name == cont)
    
    byloc <- split(sub, sub$Marker)
    
    for (loc in names(byloc)) {
      alleles <- byloc[[loc]]$Allele
      refdat[[loc]][[j]] <- alleles
      names(refdat[[loc]])[j] <- paste0("genref", j)
    }
  }
  
  
  
  
  # Sample one mixture per PROVEDIt mixture
  sampled <- genDataset(
    nC = tracerow$NoC, #true NOC
    popFreq = popfreqs2, #our population frequencies
    mu = tracerow$Hu_expectedPeakHeight, #mu from fitted DNAStatistX model
    sigma = tracerow$Hu_peakHeightVariance, #coef of variance from fitted DNAStatistX model
    threshT = c(
      "D3S1358"=35,
      "vWA"=35,
      "D16S539"=35,
      "CSF1PO"=35, 
      "TPOX"=35, 
      "D8S1179"=65, 
      "D21S11"=65, 
      "D18S51"=65, 
      "D2S441"=45, 
      "D19S433"=45, 
      "TH01"=45, 
      "FGA"=45, 
      "D22S1045"=50, 
      "D5S818"=50, 
      "D13S317"=50, 
      "D7S820"=50, 
      "SE33"=50, 
      "D10S1248"=60, 
      "D1S1656"=60,
      "D12S391"=60,
      "D2S1338"=60
    ), #Required allele peak height in mixture
    refData = refdat, #contributor genotypes
    mx = props, #mixing proportions from PROVEDIT file names
    nrep = 1, #1 replicate
    stutt = 0, #-1 stutter proportion
    stuttFW = 0, #+1 stutter proportion
    prC = 0.05, #numerical dropin probability, default in DNAStatistX
    lambda = 0.01, #The rate parameter in the exponential distribution for simulating drop-in peak heights
    beta = tracerow$Hu_degradationSlope, # Degradation slope parameter from fitted DNAStatistX model
    kit = "GlobalFiler"
  )
  
  #plotEPG2(sampled$samples,kit = "GlobalFiler",AT=0) # to visualize
  
  mixture <- map_df(names(sampled$samples), function(samp) {
    map_df(names(sampled$samples[[samp]]), function(loc) {
      dat <- sampled$samples[[samp]][[loc]]
      tibble(
        Marker = loc,
        Allele = dat$adata,
        Height = dat$hdata
      )
    })
  }) %>% 
    long_to_wide(samplelabel = paste0("EFM-", gsub(".csv", "", tracerow$trace)))
  
  if(nrow(mixture)!=21){ #then we have missing markers
    missingmarkers00 <- unique(input_1p$Marker)
    missingmarkers01 <- missingmarkers00[!missingmarkers00 %in% unique(mixture$Marker)]
    for(mark in missingmarkers01){
      mixture <- mixture %>% 
        bind_rows(data.frame(
          SampleName = unique(mixture$SampleName),
          Marker = mark
        ))
    }
  }
  
  mixture <- mixture %>% 
    replace(., is.na(.), "") %>% 
    arrange(Marker) %>% 
    mutate(Marker = ifelse(Marker=="VWA", "vWA", Marker))
  
  write.csv(mixture, paste0(output, tracerow$NoC, "p_", type, "/", unique(mixture$SampleName), ".csv"),
            row.names=F, quote = F)
}










