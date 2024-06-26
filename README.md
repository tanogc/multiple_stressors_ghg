# Code to reproduce the analyses of "Multiple stressors alter greenhouse gas concentrations in streams through local and distal processes"

Using data from an extensive survey of 50 sites across multiple stressor gradients in North Portugal, we investigated the combined effects of stressors on stream metabolism, pCO2 and pCH4 and the role of local and distal processes in driving these responses. To do this, we conducted simultaneous measurements of stream metabolism (gross primary production, ecosystem respiration and net ecosystem production) from diel changes in dissolved oxygen and algal production, and stream water concentrations of CO2 and CH4. In our models, we considered the effects of four stressors (nutrient enrichment, oxygen depletion, thermal stress and riparian degradation) and hydrology, which is a crucial factor for stream biogeochemistry.   

This code re-creates analysis of the field and experimental data

## Original article:

Please, use this citation to reference the database:
```
Gutiérrez-Cánovas, C., von Schiller, D., Pace, G., Gómez-Gener, L., & Pascoal, C. Multiple stressors alter greenhouse gas
concentrations in streams through local and distal processes. Global Change Biology (accepted)
```

# R files description:

* **0_calculating DO_def & GHG ratios.R:** R code to calculate DO deficit and GHG ratios
* **0_part_r2.R:** R function to perform variance partition for multiple models
* **1_preparing&exploring_data.R:** R code to prepare data and produce exploratory plots
* **2_multimodel_inference.R:** R script to explore the effects and importance of stressors and hydrology on stream metabolism and GHG concentrations
* **3_SEMs.R:** R script to assess the scales at which stressors influence the local stream concentrations of CO2 and CH4
* **4_scenarios.R:** R script to generate the CO2-equivalent concentrations for different combinations of dissolved inorganic nitrogen concentrations an dissolved oxygen deficit.

# Database description:
* **dat.txt:** data used to produce the results presented in the study, including stressors, discharge, metabolic rates and concentrations and emissions of greenhouse gases. It also includes information for the calculation of the reaereation rate K600 and the metabolic rates
* **read_me.xlsx:** Dataset (sheet dat) dataset used to reproduce the results and description and units of the variables (sheet variable codes).
* **folder dissolved_oxygen_data:** measures of dissolved oxygen and water temperature used to calculate the metabolic rates of each site (site01.txt to site50.txt)


```
Please, for questions related with the use of this code, send an email to cayetano.gutierrez@urjc.es.

