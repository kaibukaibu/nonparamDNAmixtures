library(tidyverse)
library(euroformix)
set.seed(66)

#setwd("..")


input_parameters_fileloc <- "LR_results/replicates_unfiltered001riman.csv"
output <- "data/data_replicated/"
pop <- "data/NIST1036_Cauc.csv"



# # make results subfolders
# if (dir.exists(output)) {
#   unlink(output, recursive = TRUE, force = T)  # Deletes the folder and all contents
#   dir.create(output)
# }
for(NOC in 2:4){
  for(type in c("EFM")){
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


# dataset of popfreqs
popfreqs <- read.csv(pop) %>% 
  pivot_longer(cols=2:25) %>%
  rename(freq = value,
         Marker = name) %>%
  filter(!is.na(freq))

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
  
  sampled <- genDataset(
    nC = tracerow$NoC,
    popFreq = popfreqs2,
    mu = tracerow$Hu_expectedPeakHeight,
    sigma = tracerow$Hu_peakHeightVariance,
    threshT = 0,
    #refData = #somehow add contributors,
    mx = props,
    nrep = 1,
    stutt = 0,
    stuttFW = 0,
    prC = 0.05, #default of DNAStatistX
    lambda = 0.01,
    beta = tracerow$Hu_degradationSlope,
    kit = "GlobalFiler"
  )
  # plotEPG2(sampled$samples,
  #          kit = "GlobalFiler",
  #          refData=sampled$refData,
  #          AT=50)
}

















