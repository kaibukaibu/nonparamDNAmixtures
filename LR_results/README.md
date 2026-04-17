## Likelihood ratio results

The LR calculation results file in this folder consists of LR calculations of 4626 mixture-PoI (who were true contributors) combinations. Note that in this dataset, a few non-parametrically simulated traces are given as duplicates, as they have two corresponding PROVEDIt mixtures. This is due to limitations in the sampling space of the 1p trace profiles. For example, for PROVEDIt mixtures "A02-RD14-0003-44-45-1-1-M4I35-0062GF-QLAND-0115sechid.csv" and "F01-RD14-0003-44-45-1-1-M4I22-0062GF-Q20-0615sechid.csv", there is only one non-parametric replicate trace profile because of these two mixture profiles being very similar (total contribution rfus 10059 and 10031; same contributors, treatment, mixing proportions, and template mass). Our simulation algorithm chose the exact same 1p profiles for both of these mixtures due to the high similarity.

The variables in the dataset represent the following:

- trace: name of mixture (either replicated or from PROVEDIt); for non-parametrically replicated mixtures the trace name includes the indicators for two contributors, indicators for used 1p origin profiles (contr, dilution number, treatment), template masses (delimiter removed, i.e. 0125 indicates 0.125 ng) and contribution rfus -- all in the same order as the contributors. For the EFM-replicated mixtures the trace name is "EFM-" pasted together with the target PROVEDIt mixture.
- poi: PoI
- NoC
- label: mixture type -- "mod" indicates replicated (modified) mixtures, "real" indicates PROVEDIt mixtures, and "EFM" indicates EuroForMix mixtures
- contr1-contr4: contributors
- temp1-temp4: template mass of contributors (for nonparam replicates it is the template masses of the 1p origin traces; for PROVEDIt mixtures it is the total template mass split according to given mixing proportions; for EFM mixtures it is the same as their targeted PROVEDIt mixtures)
- rfu_allestut1-rfu_allestut4: rfu contribution (allelic and stutter) of contributors for nonparam replicate mixtures
- rfu_allestut_total: sum of all peaks in allelic or -1 stutter positions (over all contributors)
- rfu_allestut_total_real: rfu_allestut_total of the PROVEDIt mixture that the replicate mixture was replicated for
- poi_iswho: who is the PoI? in this case they are only the true contributors
- whoserel: who is the PoI related to? in this case they are only the true contributors
- threshold: which analytic thresholds were applied, in this case we used the same dye-specific ATs as in the article by [Riman, Iyer and Vallone](https://doi.org/10.1371/journal.pone.0256714)
- filt: was filtered or unfiltered 1p data used for simulation? here we only used unfiltered
- Fst
- real_trace: corresponding PROVEDIt trace
- Hc_prop1-Hc_prop4: estimated contribution proportions assuming Hc, the first prop is for PoI and rest is ordered decreasingly
- Hu_prop1-Hu_prop4: estimated contribution proportions assuming Hu, here all contributor proportions are ordered decreasingly
- Hc_... and Hu_...: other estimated model parameters and likelihoods obtained from DNAStatistX for both hypotheses
