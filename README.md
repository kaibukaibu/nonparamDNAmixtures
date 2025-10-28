# repDNAmixtures
Code for non-parametric simulation/replication of DNA mixtures from one-person trace profiles by Kai Budrikas (kaib@itu.dk). This code is part of paper "Non-Parametric Simulation of DNA Mixture Profiles from One-Person Trace Profiles" (2026) by Kai Budrikas and Klaas Slooten. In this code we replicate 2-4-person DNA mixtures from PROVEDIt by pasting together several 1-person DNA profiles. 

# Necessary software and packages:
- **R** (tested in version 4.4.1)
- **tidyverse**
- **doParallel**
- **openxlsx**

# Explanation of files
The code consists of three R files:
- **01_replicate_mixtures.R**: code for replicating PROVEDIt mixtures. Outputs raw mixture profiles without applying stutter filter and analytical thresholds.
- **02_apply_stutterfilter.R**: code for applying analytical thresholds from Riman et al. (2021) and a stutter filter. Also generates collections of trace-POI combinations to go through in LR calculations (by whatever software).
- **helping_functions.R**: collection of functions for transforming profile data from wide to long format and vice versa.

In addition, there are three folders for data:
- **data_provedit**: mixture profiles, 1-person profiles and genotype profiles from PROVEDIt.
- **data_simulated**: replicate mixture profiles that the code generates and that are further examined in the paper. "real" represents PROVEDIt data, "nomod" represents replicate mixtures pasted from unmodified 1-person profiles (set 1 in paper) and "mod" represents mixtures pasted from modified 1-person profiles (set 2 in paper).
- **LR_results**: collection of LR results obtained in DNAStatistX for the replicated mixtures.

