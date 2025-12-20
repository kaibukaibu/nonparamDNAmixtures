The LR calculation results file in this folder contains of LR calculations of 3084 mixture-PoI combinations. The variables in the dataset represent the following:

- trace: name of mixture (either replicated or from PROVEDIt); for replicated mixtures the trace name includes the indicators for two contributors, indicators for used 1p origin profiles (contr, dilution number, treatment), template masses (delimiter removed, i.e. 0125 indicates 0.125 ng) and contribution rfus -- all in the same order as the contributors
- poi: PoI
- NoC
- label: mixture type -- "mod" indicates replicated mixtures and "real" indicates PROVEDIt mixtures
- contr1-contr4: contributors
- temp1-temp4: template mass of contributors (for replicates it is the template masses of the 1p origin traces; for PROVEDIt mixtures it is the template mass divided by given mixing proportions)
- rfu_allestut1-rfu_allestut4: rfu contribution of contributors for replicate mixtures
- rfu_allestut_total: sum of all rfu_allestut contributions
- rfu_allestut_total_real: rfu contribution of the PROVEDIt mixture that the replicate mixture was replicated for
- poi_iswho: who is the PoI to the true contributor? in this case they are only the true contributors
- whoserel: who is the PoI related to? in this case the true contributors
- threshold: which analytic thresholds were applied, in this case the ones from Riman et al. paper
- filt: was filtered or unfiltered 1p data used for simulation?
- Fst
- real_trace: corresponding PROVEDIt trace
- Hc_prop1-Hc_prop4: estimated contribution proportions assuming Hc, the first prop is for PoI and rest is ordered decreasingly
- Hu_prop1-Hu_prop4: estimated contribution proportions assuming Hu, here all contributor proportions are ordered decreasingly
- Hc_... and Hu_...: other estimated model parameters and likelihoods obtained from DNAStatistX for both hypotheses
