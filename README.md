## Overview
R code for analyses and figures used as part of this project/manuscript. Project examines the additional feature variables and there utility in predictive models for disease severity in influenza infected ferret animal models. 

See Manuscripts and Abstract below for further details. 

## Data
Data available on data.cdc.gov: https://data.cdc.gov/National-Center-for-Immunization-and-Respiratory-D/An-aggregated-dataset-of-serially-collected-influe/cr56-k9wj/about_data

Kieran TJ, Sun X, Creager HM, Tumpey TM, Maines TR, Belser JA. An aggregated dataset of serial morbidity and titer measurements from influenza A virus-infected ferrets. Sci Data 11, 510 (2024). https://doi.org/10.1038/s41597-024-03256-6

## Scripts

1_morbidity2_data_cleaning.R - R code for processing and merging the data input. Intial vairable creation and calculations.

2_morbidity2_analysis.R - R code for the primary analyses conducting, including machine learning, elastic net, and linear mixed effects models.

3_morbidity2_figures.R - R code for figures used in the manuscript. 

## Manuscript
Kieran TJ, Maines TR, Belser JA. (Under Revision). Limited ‘heft’ of weight-based outcomes in predicting influenza A virus disease severity in ferrets. PLoS Computational Biology

### Abstract
Studies evaluating viral pathogenicity in small mammalian models often quantify disease severity using the magnitudes of temperature rise and weight loss post-challenge. However, no rigorous assessment on the transformation of serially collected data into features suitable for predictive models has been conducted. Using data aggregated from ferrets inoculated with a diverse panel of influenza A viruses (IAV) spanning a broad range of clinical outcomes, we assessed statistical correlations and predictive performance of temperature and weight loss, summarized by conventional and novel approaches. Conventional summary metrics (peak values or area under the curve) were weak and inconsistent correlates of overall disease severity and viral titers. Novel dynamic weight metrics capturing onset, duration, slope, and volatility over 14 days showed lower coefficients of variation than conventional summary approaches. However, inclusion of novel metrics did not meaningfully improve the predictive performance of machine learning models for disease severity outcomes in IAV-inoculated ferrets. Mixed-effects models indicated that weight loss post-IAV infection is driven by time and viral burden, with temperature contributing little additional information. Collectively, these findings support that derived metrics are at least comparable, if not enhanced, to conventional summaries for data science analyses of serially generated clinical data from in vivo pathogen studies. However, because pathogen disease severity in mammals is multifactorial, models that rely solely on weight and temperature metrics without additional quantitative measures of clinical perturbation within-host are unlikely to achieve strong predictive performance. 
