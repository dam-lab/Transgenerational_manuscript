This file is meant as a descriptive README for the scripts and data files included for "Rapid, but limited, zooplankton adaptation to simultaneous warming and acidification".


### Notes for scripts and Data ####

To use these scripts, first download R (https://www.r-project.org/) and R Studio (https://www.rstudio.com/) for the most effective implementation of the scripts.

The Script titled: "Phenotypic_analysis_all_traits.R" is used to analyze all life-history traits discussed in the manuscript (survival, egg production rate, egg hatching success, development time, sex ratio, and population fitness). Survival, development time, and sex ratio analyses are all calculated from the available survival data  "Survival_data_total.txt" file. Egg production and egg hatching success analyses are all performed using the "EPR_HF_data_total.txt" file. Resulting population fitness data, along with the corresponding values of survival, egg production, hatching success, development time and sex ratio used to calculate fitness can be accessed directly from the "lambda_results_devtime_surv_epr_hf_sex_standardized_relative.txt" file. This also includes relative fitness values and standardized fitness values.

The Script titled: "A_tonsa_physical_data.R" is used to analyze all temperature, pH, pCO2, alkalinity, and accompanying carbonate chemistry data collected during the transgenerational experiment. The resulting data is summarized in the MS Excel workbook titled "A_tonsa_physical_data_complete_MS.xls" file. Included, you will find all the collected temperature and pH measurements, as well as all calculated values of pCO2, alkalinity, omega CA, omega AR, fCO2, and DIC. There is also all the statisically evaluated temperature and pH contrast results.

The script used to evaluate genetic (nucleotide) diversity is titled "genetic_diversity.md". It includes sections of bash script for bioinformatic analysis and subsequent R code for analysis of nucleotide diversity.

Thank you for expressing interest in our research! For an easy to understand video that explains more about our research, check out: https://youtu.be/YrI2188-ejM
