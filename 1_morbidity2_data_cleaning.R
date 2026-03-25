########################################################################################
########################################################################################
### Morbidity2 Data Cleaning
### Updated to include new features & stats
### 26 November 2025
### Troy J. Kieran
########################################################################################
########################################################################################
### Load packages

library(tidyverse)
library(tidylog)
library(janitor)
library(DescTools) ## for AUC

########################################################################################
### Import Data
fullData <- ## data available on data.cdc.gov:
            ## https://data.cdc.gov/National-Center-for-Immunization-and-Respiratory-D/An-aggregated-dataset-of-serially-collected-influe/cr56-k9wj/about_data
            ## Cite Data Descriptor: https://www.nature.com/articles/s41597-024-03256-6
morbidity2 <- ## data available from corresponding authors

########################################################################################

### Data Clean morbidity2 File
## impute means for NAs with values on either side, except euth columns

## identify tp_ and wt_ columns
tp_cols <- names(morbidity2)[str_detect(names(morbidity2), "^tp_")]
wt_cols <- names(morbidity2)[str_detect(names(morbidity2), "^wt_")]
all_cols <- c(tp_cols, wt_cols)

## create euth mask BEFORE converting anything
euth_mask <- morbidity2 %>%
  select(all_of(all_cols)) %>%
  mutate(across(everything(), ~ .x == "euth"))

## interpolation mean helper function
interpolate_mean <- function(values, euth_flags) {
  euth_flags <- ifelse(is.na(euth_flags), FALSE, euth_flags)
  for (i in seq_along(values)) {
    if (is.na(values[i]) && !euth_flags[i]) {
      left <- if (i > 1) values[i - 1] else NA
      right <- if (i < length(values)) values[i + 1] else NA
      if (!is.na(left) && !is.na(right)) {
        values[i] <- mean(c(left, right), na.rm = TRUE)}}}
  return(values)}

## interpolation function
interpolate_block <- function(df, cols, euth_mask) {
  numeric_data <- df[cols] %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(replace(.x, .x == "euth", NA)))))
  
  filled <- map2_dfr(
    split(numeric_data, seq(nrow(numeric_data))),
    split(euth_mask[cols], seq(nrow(euth_mask))),
    ~ interpolate_mean(unlist(.x), unlist(.y)) %>% as.list())
  
  names(filled) <- cols
  df[cols] <- filled
  return(df)}

## apply interpolation to tp_ and wt_ columns
morbidity2 <- interpolate_block(morbidity2, tp_cols, euth_mask)
morbidity2 <- interpolate_block(morbidity2, wt_cols, euth_mask)

## add euth_flag and euth_day columns
morbidity2 <- morbidity2 %>%
  mutate(
    euth_flag = apply(euth_mask, 1, function(row) any(row, na.rm = TRUE)),
    euth_day = apply(euth_mask, 1, function(row) {
      euth_indices <- which(row)
      if (length(euth_indices) > 0) {
        as.numeric(str_extract(names(row)[euth_indices[1]], "\\d+"))} 
      else {NA}}))

## rename column for convenience
morbidity2 <- morbidity2 %>% rename(euth = euth_flag)

## normalize tp_ columns (change from tp_0)
tp_norm <- morbidity2 %>%
  select(all_of(tp_cols)) %>%
  mutate(across(everything(), ~ . - tp_0, .names = "{.col}_n"))

## normalize wt_ columns (percent change from wt_0)
wt_norm <- morbidity2 %>%
  select(all_of(wt_cols)) %>%
  mutate(across(everything(), ~ (. / wt_0) * 100, .names = "{.col}_n"))

## actual percent change
wt_norm_pct <- wt_norm %>%
  mutate(across(everything(), ~ . - 100))

## extract the normalized columns only
tp_norm <- tp_norm %>% select(ends_with("_n"))
wt_norm <- wt_norm %>% select(ends_with("_n"))
wt_norm_pct <- wt_norm_pct %>% select(ends_with("_n"))

## bind normalized columns to original data
morbidity2 <- bind_cols(morbidity2, tp_norm, wt_norm_pct)

#str(morbidity2)
#View(morbidity2)

########################################################################################

### merge fullData and morbidity2

fullMorbidity <- left_join(fullData, morbidity2, by = c("Ferret")) %>%
  subset(., select = -c(Virus.y, CDCID.y)) %>%
  rename(c(Virus = Virus.x, CDCID = CDCID.x))

#str(fullMorbidity)
#View(fullMorbidity)

########################################################################################

## recreate temp and temp_5 parameters using morbidity dataset

fullMorbidity$tp_max <- fullMorbidity %>%
  dplyr::select(tp_0_n:tp_14_n) %>%
  apply(., 1, function(x) max(x, na.rm = TRUE))

fullMorbidity$tp_max_5 <- fullMorbidity %>%
  dplyr::select(tp_0_n:tp_5_n) %>%
  apply(., 1, function(x) max(x, na.rm = TRUE))

## convert all -Inf and Inf to NA
fullMorbidity[fullMorbidity == "-Inf"] <- NA
fullMorbidity[fullMorbidity == "Inf"] <- NA

########################################################################################

### generate per-ferret AUC of weight loss & temp using different gatings

## Define a function to calculate AUC per ferret
## This function calculates the AUC per ferret for a given set of days of infection.
## The function takes two arguments:
##   - df: The dataframe containing the infection data (here it is fullMorbidity)
##   - day_inoc_cols: A vector of the column names that contain the day of infection data
## The function returns a dataframe with the AUC per ferret.

## Create object with column names for AUC function
day_inoc_cols <- c("wt_0_n", "wt_1_n", "wt_2_n", "wt_3_n", "wt_4_n", 
                   "wt_5_n", "wt_6_n", "wt_7_n", "wt_8_n", "wt_9_n",
                   "wt_10_n", "wt_11_n", "wt_12_n", "wt_13_n", "wt_14_n")

calculate_AUC_per_ferret <- function(fullMorbidity, day_inoc_cols, method = "trapezoid") { 
  ## Create a dataframe with the day of infection and the wt loss
  df_day_inoc <- fullMorbidity %>% gather(all_of(day_inoc_cols), key = day_inoc, value = inoc_wt) %>%
    mutate(day = case_when(
      day_inoc == "wt_0_n" ~ 0,
      day_inoc == "wt_1_n" ~ 1,
      day_inoc == "wt_2_n" ~ 2,
      day_inoc == "wt_3_n" ~ 3,
      day_inoc == "wt_4_n" ~ 4,
      day_inoc == "wt_5_n" ~ 5,
      day_inoc == "wt_6_n" ~ 6,
      day_inoc == "wt_7_n" ~ 7,
      day_inoc == "wt_8_n" ~ 8,
      day_inoc == "wt_9_n" ~ 9,
      day_inoc == "wt_10_n" ~ 10,
      day_inoc == "wt_11_n" ~ 11,
      day_inoc == "wt_12_n" ~ 12,
      day_inoc == "wt_13_n" ~ 13,
      day_inoc == "wt_14_n" ~ 14)) 
  
  ## Calculate the AUC for each ferret
  auc_per_ferret <- sapply(unique(df_day_inoc$Ferret), function(x) {
    AUC(df_day_inoc[df_day_inoc$Ferret == x, ]$day, 
        df_day_inoc[df_day_inoc$Ferret == x, ]$inoc_wt, 
        method = method,
        na.rm = TRUE)})
  
  ## Create a dataframe with the AUC per ferret
  df_auc_per_ferret <- data.frame(Ferret = unique(df_day_inoc$Ferret), AUC_f = auc_per_ferret)
  return(df_auc_per_ferret)} 

## Calculate AUC per ferret for 14 days of infection
df_auc_14_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols) %>%
  rename(AUC_14_wt = AUC_f)
## Calculate AUC per ferret for 13 days of infection
df_auc_13_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:14]) %>%
  rename(AUC_13_wt = AUC_f)
## Calculate AUC per ferret for 12 days of infection
df_auc_12_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:13]) %>%
  rename(AUC_12_wt = AUC_f)
## Calculate AUC per ferret for 11 days of infection
df_auc_11_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:12]) %>%
  rename(AUC_11_wt = AUC_f)
## Calculate AUC per ferret for 10 days of infection
df_auc_10_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:11]) %>%
  rename(AUC_10_wt = AUC_f)
## Calculate AUC per ferret for 9 days of infection
df_auc_9_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:10]) %>%
  rename(AUC_9_wt = AUC_f)
## Calculate AUC per ferret for 8 days of infection
df_auc_8_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:9]) %>%
  rename(AUC_8_wt = AUC_f)
## Calculate AUC per ferret for 7 days of infection
df_auc_7_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:8]) %>%
  rename(AUC_7_wt = AUC_f)
## Calculate AUC per ferret for 6 days of infection
df_auc_6_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:7]) %>%
  rename(AUC_6_wt = AUC_f)
## Calculate AUC per ferret for 5 days of infection
df_auc_5_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:6]) %>%
  rename(AUC_5_wt = AUC_f)
## Calculate AUC per ferret for 4 days of infection
df_auc_4_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:5]) %>%
  rename(AUC_4_wt = AUC_f)
## Calculate AUC per ferret for 3 days of infection
df_auc_3_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:4]) %>%
  rename(AUC_3_wt = AUC_f)
## Calculate AUC per ferret for 2 days of infection
df_auc_2_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:3]) %>%
  rename(AUC_2_wt = AUC_f)
## Calculate AUC per ferret for 1 day of infection
df_auc_1_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:2]) %>%
  rename(AUC_1_wt = AUC_f)

## join the AUC per ferret dataframes to the fullData dataframe
fullMorbidity <- fullMorbidity %>% 
  left_join(df_auc_1_f, by = "Ferret") %>%
  left_join(df_auc_2_f, by = "Ferret") %>%
  left_join(df_auc_3_f, by = "Ferret") %>%
  left_join(df_auc_4_f, by = "Ferret") %>%
  left_join(df_auc_5_f, by = "Ferret") %>%
  left_join(df_auc_6_f, by = "Ferret") %>%
  left_join(df_auc_7_f, by = "Ferret") %>%
  left_join(df_auc_8_f, by = "Ferret") %>%
  left_join(df_auc_9_f, by = "Ferret") %>%
  left_join(df_auc_10_f, by = "Ferret") %>%
  left_join(df_auc_11_f, by = "Ferret") %>%
  left_join(df_auc_12_f, by = "Ferret") %>%
  left_join(df_auc_13_f, by = "Ferret") %>%
  left_join(df_auc_14_f, by = "Ferret")

## add NA's back into data to reflect ferrets euth during observation period
fullMorbidity <- 
  fullMorbidity %>%
  mutate(AUC_14_wt = ifelse(euth == TRUE, NA, AUC_14_wt),
         AUC_13_wt = ifelse((euth == TRUE & euth_day < 14), NA, AUC_13_wt),
         AUC_12_wt = ifelse((euth == TRUE & euth_day < 13), NA, AUC_12_wt),
         AUC_11_wt = ifelse((euth == TRUE & euth_day < 12), NA, AUC_11_wt),
         AUC_10_wt = ifelse((euth == TRUE & euth_day < 11), NA, AUC_10_wt),
         AUC_9_wt = ifelse((euth == TRUE & euth_day < 10), NA, AUC_9_wt),
         AUC_8_wt = ifelse((euth == TRUE & euth_day < 9), NA, AUC_8_wt),
         AUC_7_wt = ifelse((euth == TRUE & euth_day < 8), NA, AUC_7_wt),
         AUC_6_wt = ifelse((euth == TRUE & euth_day < 7), NA, AUC_6_wt),
         AUC_5_wt = ifelse((euth == TRUE & euth_day < 6), NA, AUC_5_wt),
         AUC_4_wt = ifelse((euth == TRUE & euth_day < 5), NA, AUC_4_wt),
         AUC_3_wt = ifelse((euth == TRUE & euth_day < 4), NA, AUC_3_wt),
         AUC_2_wt = ifelse((euth == TRUE & euth_day < 3), NA, AUC_2_wt),
         AUC_1_wt = ifelse((euth == TRUE & euth_day < 2), NA, AUC_1_wt))

## add NA's back into data to reflect missing weight values not collected
fullMorbidity <- fullMorbidity %>%
  rowwise() %>%
  mutate(across(
    .cols = paste0("AUC_", 1:14, "_wt"),
    .fns = ~ if (any(is.na(c_across(paste0("wt_", 0:as.numeric(str_extract(cur_column(), "\\d+")), "_n"))))) NA else .,
    .names = "{.col}")) %>%
  ungroup()

###
### repeat for temp
###

## Create object with column names for AUC function
day_inoc_cols <- c("tp_0_n", "tp_1_n", "tp_2_n", "tp_3_n", "tp_4_n", 
                   "tp_5_n", "tp_6_n", "tp_7_n", "tp_8_n", "tp_9_n",
                   "tp_10_n", "tp_11_n", "tp_12_n", "tp_13_n", "tp_14_n")

###run this function in one go otherwise R says NO
calculate_AUC_per_ferret <- function(fullMorbidity, day_inoc_cols, method = "trapezoid") { 
  ## Create a dataframe with the day of infection and the wt loss
  df_day_inoc <- fullMorbidity %>% gather(all_of(day_inoc_cols), key = day_inoc, value = inoc_tp) %>%
    mutate(day = case_when(
      day_inoc == "tp_0_n" ~ 0,
      day_inoc == "tp_1_n" ~ 1,
      day_inoc == "tp_2_n" ~ 2,
      day_inoc == "tp_3_n" ~ 3,
      day_inoc == "tp_4_n" ~ 4,
      day_inoc == "tp_5_n" ~ 5,
      day_inoc == "tp_6_n" ~ 6,
      day_inoc == "tp_7_n" ~ 7,
      day_inoc == "tp_8_n" ~ 8,
      day_inoc == "tp_9_n" ~ 9,
      day_inoc == "tp_10_n" ~ 10,
      day_inoc == "tp_11_n" ~ 11,
      day_inoc == "tp_12_n" ~ 12,
      day_inoc == "tp_13_n" ~ 13,
      day_inoc == "tp_14_n" ~ 14))
  
  ## Calculate the AUC for each ferret
  auc_per_ferret <- sapply(unique(df_day_inoc$Ferret), function(x) {
    AUC(df_day_inoc[df_day_inoc$Ferret == x, ]$day,
        df_day_inoc[df_day_inoc$Ferret == x, ]$inoc_tp,
        na.rm = TRUE)}) ### right before na.rm I deleted method = method,
  
  ## Create a dataframe with the AUC per ferret
  df_auc_per_ferret <- data.frame(Ferret = unique(df_day_inoc$Ferret), AUC_f = auc_per_ferret)
  return(df_auc_per_ferret)}

## Calculate AUC per ferret for 14 days of infection
df_auc_14_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols) %>%
  rename(AUC_14_tp = AUC_f)
## Calculate AUC per ferret for 13 days of infection
df_auc_13_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:14]) %>%
  rename(AUC_13_tp = AUC_f)
## Calculate AUC per ferret for 12 days of infection
df_auc_12_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:13]) %>%
  rename(AUC_12_tp = AUC_f)
## Calculate AUC per ferret for 11 days of infection
df_auc_11_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:12]) %>%
  rename(AUC_11_tp = AUC_f)
## Calculate AUC per ferret for 10 days of infection
df_auc_10_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:11]) %>%
  rename(AUC_10_tp = AUC_f)
## Calculate AUC per ferret for 9 days of infection
df_auc_9_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:10]) %>%
  rename(AUC_9_tp = AUC_f)
## Calculate AUC per ferret for 8 days of infection
df_auc_8_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:9]) %>%
  rename(AUC_8_tp = AUC_f)
## Calculate AUC per ferret for 7 days of infection
df_auc_7_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:8]) %>%
  rename(AUC_7_tp = AUC_f)
## Calculate AUC per ferret for 6 days of infection
df_auc_6_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:7]) %>%
  rename(AUC_6_tp = AUC_f)
## Calculate AUC per ferret for 5 days of infection
df_auc_5_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:6]) %>%
  rename(AUC_5_tp = AUC_f)
## Calculate AUC per ferret for 4 days of infection
df_auc_4_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:5]) %>%
  rename(AUC_4_tp = AUC_f)
## Calculate AUC per ferret for 3 days of infection
df_auc_3_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:4]) %>%
  rename(AUC_3_tp = AUC_f)
## Calculate AUC per ferret for 2 days of infection
df_auc_2_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:3]) %>%
  rename(AUC_2_tp = AUC_f)
## Calculate AUC per ferret for 1 day of infection
df_auc_1_f <- calculate_AUC_per_ferret(fullMorbidity = fullMorbidity, day_inoc_cols = day_inoc_cols[0:2]) %>%
  rename(AUC_1_tp = AUC_f)

## Join the AUC per ferret dataframes to the fullData dataframe
fullMorbidity <- fullMorbidity %>%
  left_join(df_auc_1_f, by = "Ferret") %>%
  left_join(df_auc_2_f, by = "Ferret") %>%
  left_join(df_auc_3_f, by = "Ferret") %>%
  left_join(df_auc_4_f, by = "Ferret") %>%
  left_join(df_auc_5_f, by = "Ferret") %>%
  left_join(df_auc_6_f, by = "Ferret") %>%
  left_join(df_auc_7_f, by = "Ferret") %>%
  left_join(df_auc_8_f, by = "Ferret") %>%
  left_join(df_auc_9_f, by = "Ferret") %>%
  left_join(df_auc_10_f, by = "Ferret") %>%
  left_join(df_auc_11_f, by = "Ferret") %>%
  left_join(df_auc_12_f, by = "Ferret") %>%
  left_join(df_auc_13_f, by = "Ferret") %>%
  left_join(df_auc_14_f, by = "Ferret")

## Add NA's back into data to reflect ferrets euth during observation period
fullMorbidity <-
  fullMorbidity %>%
  mutate(AUC_14_tp = ifelse(euth == TRUE, NA, AUC_14_wt),
         AUC_13_tp = ifelse((euth == TRUE & euth_day < 14), NA, AUC_13_tp),
         AUC_12_tp = ifelse((euth == TRUE & euth_day < 13), NA, AUC_12_tp),
         AUC_11_tp = ifelse((euth == TRUE & euth_day < 12), NA, AUC_11_tp),
         AUC_10_tp = ifelse((euth == TRUE & euth_day < 11), NA, AUC_10_tp),
         AUC_9_tp = ifelse((euth == TRUE & euth_day < 10), NA, AUC_9_tp),
         AUC_8_tp = ifelse((euth == TRUE & euth_day < 9), NA, AUC_8_tp),
         AUC_7_tp = ifelse((euth == TRUE & euth_day < 8), NA, AUC_7_tp),
         AUC_6_tp = ifelse((euth == TRUE & euth_day < 7), NA, AUC_6_tp),
         AUC_5_tp = ifelse((euth == TRUE & euth_day < 6), NA, AUC_5_tp),
         AUC_4_tp = ifelse((euth == TRUE & euth_day < 5), NA, AUC_4_tp),
         AUC_3_tp = ifelse((euth == TRUE & euth_day < 4), NA, AUC_3_tp),
         AUC_2_tp = ifelse((euth == TRUE & euth_day < 3), NA, AUC_2_tp),
         AUC_1_tp = ifelse((euth == TRUE & euth_day < 2), NA, AUC_1_tp))

## Add NA's back into data to reflect missing temp values not collected
fullMorbidity <- fullMorbidity %>%
  rowwise() %>%
  mutate(across(
    .cols = paste0("AUC_", 1:14, "_tp"),
    .fns = ~ if (any(is.na(c_across(paste0("tp_", 0:as.numeric(str_extract(cur_column(), "\\d+")), "_n"))))) NA else .,
    .names = "{.col}")) %>%
  ungroup()

########################################################################################

## longer for only wt_#_n columns
long_wt <- fullMorbidity %>%
  pivot_longer(cols = matches("^wt_\\d+_n$"),
               names_to = "day_label",
               values_to = "pct_loss") %>%
  mutate(day = as.numeric(str_extract(day_label, "\\d+")))

## find day of max weight loss (most negative pct_loss)
max_loss_day <- long_wt %>%
  group_by(Ferret) %>%
  ## filtering rows where pct_loss equals the min value for that ferret
  filter(pct_loss == min(pct_loss, na.rm = TRUE)) %>%
  ## then 'select' that row w/ slice
  slice(1) %>%
  ungroup() %>%
  select(Ferret, max_wt_loss_day = day, max_wt_loss_pct = pct_loss)

## find first day reaching ≤ –5% weight loss
first_5pct_day <- long_wt %>%
  group_by(Ferret) %>%
  ## 5% threshold
  filter(pct_loss <= -5) %>%
  ## selects the earliest day where that threshold was met
  slice_min(day) %>%
  ungroup() %>%
  select(Ferret, day_5pct_wt_loss = day, value_5pct_wt_loss = pct_loss)

## join back into fullMorbidity and calculate slopes
fullMorbidity <- fullMorbidity %>%
  left_join(max_loss_day, by = "Ferret") %>%
  left_join(first_5pct_day, by = "Ferret") %>%
  ## slope of max weight loss & day of max weight loss
  mutate(wt_loss_slope_max = (max_wt_loss_pct / max_wt_loss_day)) %>%
  ## slope between difference in max weight loss & 5% weight loss
  mutate(wt_loss_slope_max5 = (max_wt_loss_pct - value_5pct_wt_loss) / 
           (max_wt_loss_day - day_5pct_wt_loss))

## replace NaN with NA in wt_loss_slope_max
fullMorbidity$wt_loss_slope_max[is.nan(fullMorbidity$wt_loss_slope_max)] <- NA

## replace NaN with 0 in wt_loss_slope_max5
fullMorbidity$wt_loss_slope_max5[is.nan(fullMorbidity$wt_loss_slope_max5)] <- 0

#view(fullMorbidity)

########################################################################################
########################################################################################

### export file

## save/export as CSV file w/ current date in the format YYYY-MM-DD
## this will prevent overwriting mistakes, but create 'duplicates'
# fullMorbidity %>% write.csv(paste("inputs/fullMorbidity_", format(Sys.Date(), "%Y-%m-%d"), 
#                              ".csv", sep = ""), row.names = FALSE)


########################################################################################
########################################################################################

### create new features

## helpful IDs / factors
data <- fullMorbidity %>%
  mutate(
    Ferret = as.factor(Ferret),
    Virus = as.factor(Virus),
    expt  = factor(expt, levels = c("path","DC","RD")),
    lethal = factor(lethal, levels = c("no","yes")),
    Subtype = as.factor(Subtype),
    Origin  = as.factor(Origin),
    PB2_627 = as.factor(PB2_627))

skimr::skim(data)

########################################################################################

### Reshape to long format

## temperature (raw)
tp_long <- data %>%
  pivot_longer(cols = matches("^tp_\\d+$"), names_to = "tp_day", values_to = "tp") %>%
  mutate(day = as.integer(str_extract(tp_day, "\\d+"))) %>%
  select(-tp_day)

## temperature (normalized)
tpn_long <- data %>%
  pivot_longer(cols = matches("^tp_\\d+_n$"), names_to = "tpn_day", values_to = "tp_n") %>%
  mutate(day = as.integer(str_extract(tpn_day, "\\d+"))) %>%
  select(-tpn_day)

## weight (raw)
wt_long <- data %>%
  pivot_longer(cols = matches("^wt_\\d+$"), names_to = "wt_day", values_to = "wt") %>%
  mutate(day = as.integer(str_extract(wt_day, "\\d+"))) %>%
  select(-wt_day)

## merge long tables + keep metadata
id_cols <- intersect(names(tp_long), names(tpn_long))
id_cols2 <- intersect(names(tpn_long), names(wt_long))

long <- tp_long %>%
  left_join(tpn_long, by = id_cols) %>%
  left_join(wt_long,  by = id_cols2)

## create derived % change in weight from baseline (redundant even if wt_n cols exist)
long <- long %>%
  group_by(Ferret) %>%
  mutate(wt0 = first(wt[day == 0], default = NA_real_),
         wt_pct = 100 * (wt / wt0 - 1)) %>%
  ungroup()

########################################################################################

### Create new variables for time to event/burden

## calculation functions

## AUC
trapz <- function(x, y) {    # trapezoidal rule
  ok <- which(!is.na(x) & !is.na(y))
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  sum(diff(x) * (head(y,-1) + tail(y,-1)) / 2)}

## root mean square of successive differences (volatility)
rmssd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  sqrt(mean(diff(x)^2))}

## mean absolute successive difference (robust volatility)
masd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  mean(abs(diff(x)))}

## rate normalized root mean square of successive differences (volatility)
## for time gaps (due to NAs)
rmssd_rate <- function(t, x) {
  t <- as.numeric(t); x <- as.numeric(x)
  d <- x[-1] - x[-length(x)]
  dt <- t[-1] - t[-length(t)]
  keep <- !is.na(d) & !is.na(dt) & dt > 0
  d <- d[keep]; dt <- dt[keep]
  if (length(d) == 0) return(NA_real_)
  sqrt(mean((d/dt)^2))}

## rate normalized mean absolute successive difference (robust volatility)
## for time gaps (due to NAs)
masd_rate <- function(t, x) {
  t <- as.numeric(t); x <- as.numeric(x)
  d <- x[-1] - x[-length(x)]
  dt <- t[-1] - t[-length(t)]
  keep <- !is.na(d) & !is.na(dt) & dt > 0
  d <- d[keep]; dt <- dt[keep]
  if (length(d) == 0) return(NA_real_)
  mean(abs(d/dt))}

## onset finder with run length
first_run <- function(x, day, thresh, k = 1, dir = ">=") {
  if (all(is.na(x))) return(NA_integer_)
  comp <- switch(dir, ">=" = x >= thresh, ">" = x > thresh, "<=" = x <= thresh, "<" = x < thresh)
  comp[is.na(comp)] <- FALSE
  r <- rle(comp)
  pos <- which(r$values & r$lengths >= k)[1]
  if (is.na(pos)) return(NA_integer_)
  start_idx <- sum(r$lengths[seq_len(pos-1)]) + 1
  day[start_idx]}

## lag-1 autocorrelation with guards
acf1 <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  x1 <- head(x, -1); x2 <- tail(x,  -1)
  ok <- !is.na(x1) & !is.na(x2)
  if (!any(ok)) return(NA_real_)
  cor(x1[ok], x2[ok])}

########################################################################################

## function for new feature computation per ferret
compute_features <- function(df,
                             tp_thr = 1.25,          # threshold degree °C from baseline = fever onset
                             temp_recover_tol = 0.5, # threshold degree °C within baseline = recovered
                             wt_thr_onset = -5,      # threshold % change in wt = weight loss onset
                             wt_thr_recov = -5,      # threshold % change in wt = recovered
                             k_consec = 2,           # number of consecutive days needed to confirm onset or recovery
                             early_thr = 5) {        # day that separates early vs late disease periods (days 0-5 vs 6+)
  df    <- df %>% arrange(day)
  d     <- df$day
  tpn   <- df$tp_n
  tpabs <- df$tp
  wtpc  <- df$wt_pct
  
  ## safety helpers
  day_of_max <- function(x) if (all(is.na(x))) NA_integer_ else d[which.max(x)]
  day_of_min <- function(x) if (all(is.na(x))) NA_integer_ else d[which.min(x)]
  
  ### A. Onset and Duration Metrics ###
  
  ## First day where temperature (Δ from baseline) exceeds a threshold (e.g., ≥1 °C) for at least 1 day.
  tp_onset_1d     <- first_run(tpn, d, tp_thr, k = 1, dir = ">=")
  ## Same, but requires 2 consecutive days above threshold (set above, stricter definition of onset). = 1st day of 2 day tp streak.
  tp_onset_2d     <- first_run(tpn, d, tp_thr, k = k_consec, dir = ">=")
  ## First day where weight loss exceeds a threshold (e.g., ≤−5% set above) for ≥1 day.
  wt_onset_1d     <- first_run(wtpc, d, wt_thr_onset, k = 1, dir = "<=")
  ## Same, but requires 2 consecutive days above threshold (set above, stricter definition of onset). = 1st day of 2 day wt streak.
  wt_onset_2d     <- first_run(wtpc, d, wt_thr_onset, k = k_consec, dir = "<=")
  ## Total number of days temperature increase was ≥1 or 1.5 °C.
  tp_days_1C      <- sum(tpn >= 1, na.rm = TRUE)
  tp_days_1p5C    <- sum(tpn >= 1.5, na.rm = TRUE)
  ## overall tp_0 baseline mean = 38.50971 & median 38.5
  ## max baseline = 40.1, max temp = 42.2
  ## Total number of days absolute temperature ≥39, 39.5, or 40 °C.
  tp_days_39     <- sum(tpabs >= 39, na.rm = TRUE)
  tp_days_39p5   <- sum(tpabs >= 39.5, na.rm = TRUE)
  tp_days_40     <- sum(tpabs >= 40, na.rm = TRUE)
  ## Number of days weight loss met/exceeded thresholds from −2.5% to −15%.
  #wt_days_2p5    <- sum(wtpc <= -2.5, na.rm = TRUE)
  wt_days_5_d5     <- sum(df$wt_pct[df$day %in% 1:5] <= -5, na.rm = TRUE)
  wt_days_5_d7     <- sum(df$wt_pct[df$day %in% 1:7] <= -5, na.rm = TRUE)
  wt_days_7p5_d5  <- sum(df$wt_pct[df$day %in% 1:5] <= -7.5, na.rm = TRUE)
  wt_days_7p5_d7 <- sum(df$wt_pct[df$day %in% 1:7] <= -7.5, na.rm = TRUE)
  wt_days_10_d5     <- sum(df$wt_pct[df$day %in% 1:5] <= -10, na.rm = TRUE)
  wt_days_10_d7     <- sum(df$wt_pct[df$day %in% 1:7] <= -10, na.rm = TRUE)
  wt_days_12p5_d5  <- sum(df$wt_pct[df$day %in% 1:5] <= -12.5, na.rm = TRUE)
  wt_days_12p5_d7 <- sum(df$wt_pct[df$day %in% 1:7] <= -12.5, na.rm = TRUE)
  wt_days_15_d5     <- sum(df$wt_pct[df$day %in% 1:5] <= -15, na.rm = TRUE)
  wt_days_15_d7     <- sum(df$wt_pct[df$day %in% 1:7] <= -15, na.rm = TRUE)
  
  
  ### B. Peak and Recovery Metrics ###
  
  ## peak days
  ## Day of peak fever (maximum Δ temperature).
  tpk_day <- day_of_max(tpn)
  ## Day of lowest (peak) weight (most negative % change).
  wpk_day <- day_of_min(wtpc)
  
  ## recovery
  ## Day when fever returned close to baseline (e.g., within 0.5 °C for 2+ days).
  temp_recovery_day <- {
    if (is.na(tpk_day)) NA_integer_ else {
      post <- d >= tpk_day
      cond <- df$tp_n <= temp_recover_tol ## threshold set in arguments
      cond[!post] <- FALSE
      r <- rle(cond)
      starts <- cumsum(c(1, head(r$lengths, -1)))
      pos <- which(r$values & r$lengths >= k_consec)[1]
      if (is.na(pos)) NA_integer_ else d[starts[pos]]}}
  
  ## Day when weight recovered above a set threshold (e.g., better than −2.5%).
  wt_recovery_day <- {
    if (is.na(wpk_day)) NA_integer_ else {
      post <- d >= wpk_day
      cond <- wtpc >= wt_thr_recov
      cond[!post] <- FALSE
      r <- rle(cond)
      starts <- cumsum(c(1, head(r$lengths, -1)))
      pos <- which(r$values & r$lengths >= k_consec)[1]
      if (is.na(pos)) NA_integer_ else d[starts[pos]]}}
  
  ### C. Burden / Area-Under-Curve (AUC) Metrics ###
  
  ## early days (i.e. 0-5), threshold values set in arguments
  early <- d <= early_thr
  ## Total area under the curve of fever (all days with positive ΔT).
  tp_auc_pos_all   <- trapz(d, pmax(tpn, 0))
  ## Same, restricted to early days 0–5.
  tp_auc_pos_0to5  <- trapz(d[early], pmax(tpn[early], 0))
  ## AUC of fever above threshold (e.g., 1°C, set in arguments).
  tp_auc_above_thr <- trapz(d, pmax(tpn - tp_thr, 0))
  ## Quadratic AUC (penalizes high spikes more strongly).
  tp_auc_sq        <- trapz(d, pmax(tpn, 0)^2)
  
  ## AUC of weight loss below baseline (all days).
  wt_auc_below_all    <- trapz(d, pmax(-wtpc, 0))
  ## Same, days 0–5 only.
  wt_auc_below_0to5   <- trapz(d[early], pmax(-wtpc[early], 0))
  ## AUC of weight loss below −2.5% or -5% (burden beyond a mild clinical threshold).
  wt_auc_excess_2p5   <- trapz(d, pmax(-(wtpc + 2.5), 0))
  wt_auc_excess_5     <- trapz(d, pmax(-(wtpc + 5), 0))
  ## Quadratic AUC of weight loss (heavier penalty for deeper losses).
  wt_auc_sq           <- trapz(d, (pmax(-wtpc, 0))^2)
  
  ### D. Early vs Later Days Ratios ###
  
  ## Calculate AUC for early and late days for temp and weight
  wt_early_auc    <- trapz(d[early], pmax(-wtpc[early], 0))
  wt_late_auc     <- trapz(d[!early], pmax(-wtpc[!early], 0))
  temp_early_auc  <- trapz(d[early], pmax(tpn[early], 0))
  temp_late_auc   <- trapz(d[!early], pmax(tpn[!early], 0))
  
  ## Early (0–5 days) vs late (6–14 days) AUC of weight loss. Threshold set in arguments.
  wt_early_late_ratio   <- wt_early_auc  / wt_late_auc
  ## Same for temp/fever.
  temp_early_late_ratio <- temp_early_auc / temp_late_auc
  ## High ratio (>1) = disease hits early and strong, then eases.
  ## Low ratio (<1) = slower, more prolonged burden.
  
  ### E. Variability / Dynamics ###
  
  ## rmssd: Root mean square of successive differences (day-to-day variability).
  temp_rmssd <- rmssd(tpn)
  wt_rmssd   <- rmssd(wtpc)
  temp_rmssd_rate <- rmssd_rate(d, tpn)
  wt_rmssd_rate   <- rmssd_rate(d, wtpc)
  ## masd: Mean absolute successive difference (simpler volatility measure).
  temp_masd  <- masd(tpn)
  wt_masd    <- masd(wtpc)
  temp_masd_rate  <- masd_rate(d, tpn)
  wt_masd_rate    <- masd_rate(d, wtpc)
  ## acf1: Lag-1 autocorrelation (how similar today is to yesterday; high = smoother trajectory).
  temp_acf1  <- acf1(tpn)
  wt_acf1    <- acf1(wtpc)
  ## High values = fluctuating illness course (unstable fever or weight changes).
  ## Low values =  smooth, monotonic course.
  
  ### F. Rise and Recovery Dynamics ###
  
  ## Speed of temperature increase from onset to peak.
  ##   High value = fever escalates very quickly — severe acute onset.
  ##   Low value = fever slow to reach peak
  rise_rate_temp <- {
    if (is.na(tpk_day) || is.na(tp_onset_1d)) NA_real_ else {
      (max(tpn, na.rm = TRUE) - 0) / max(1, tpk_day - tp_onset_1d)}}
  ## Speed of recovery from peak back to baseline.
  ##   High value = fever resolves quickly.
  ##   Low value = fever prolonged, slow recovery
  recov_rate_temp <- {
    if (is.na(tpk_day) || is.na(temp_recovery_day)) NA_real_ else {
      (max(tpn, na.rm = TRUE) - 0) / max(1, temp_recovery_day - tpk_day)}}
  ## Compares rise vs recovery (e.g., sharp rise, slow recovery).
  ##   >1: Rapid rise but slow recovery = more damaging illness.
  ##   ~1: Symmetrical rise and fall.
  ##   <1: Gradual rise, fast resolution.
  rise_recover_ratio_temp <- rise_rate_temp / recov_rate_temp
  
  ### G. Gather all variables into a dataframe ###
  
  tibble::tibble(
    tp_onset_1d, tp_onset_2d, wt_onset_1d, wt_onset_2d,
    tp_days_1C, tp_days_1p5C, tp_days_39, tp_days_39p5, tp_days_40,
    wt_days_5_d5, wt_days_5_d7, wt_days_7p5_d5, wt_days_7p5_d7, wt_days_10_d5, 
    wt_days_10_d7, wt_days_12p5_d5, wt_days_12p5_d7, wt_days_15_d5, wt_days_15_d7,
    tpk_day, wpk_day, temp_recovery_day, wt_recovery_day,
    tp_auc_pos_all, tp_auc_pos_0to5, tp_auc_above_thr, tp_auc_sq,
    wt_auc_below_all, wt_auc_below_0to5, wt_auc_excess_2p5, wt_auc_excess_5, wt_auc_sq,
    wt_early_late_ratio, temp_early_late_ratio,
    temp_rmssd, temp_rmssd_rate, temp_masd, temp_masd_rate, temp_acf1, 
    wt_rmssd, wt_rmssd_rate, wt_masd, wt_masd_rate, wt_acf1,
    rise_rate_temp, recov_rate_temp, rise_recover_ratio_temp
  )}

########################################################################################

### run function to calculate new variables/features

features <- long %>%
  group_by(Ferret) %>%
  group_modify(~ compute_features(.x, tp_thr = 1.25, temp_recover_tol = 0.5,
                                  wt_thr_onset = -5, wt_thr_recov = -5, 
                                  k_consec = 2, early_thr = 5)) %>%
  ungroup()

features1 <- long %>%
  group_by(Ferret) %>%
  group_modify(~ compute_features(.x, tp_thr = 1.25, temp_recover_tol = 0.5,
                                  wt_thr_onset = -7.5, wt_thr_recov = -5, 
                                  k_consec = 2, early_thr = 5)) %>%
  ungroup() %>%
  select(Ferret, wt_onset_1d, wt_onset_2d) %>%
  rename(wt_onset_1d_7p5 = wt_onset_1d,
         wt_onset_2d_7p5 = wt_onset_2d)

features2 <- long %>%
  group_by(Ferret) %>%
  group_modify(~ compute_features(.x, tp_thr = 1.25, temp_recover_tol = 0.5,
                                  wt_thr_onset = -10, wt_thr_recov = -5, 
                                  k_consec = 2, early_thr = 5)) %>%
  ungroup() %>%
  select(Ferret, wt_onset_1d, wt_onset_2d) %>%
  rename(wt_onset_1d_10 = wt_onset_1d,
         wt_onset_2d_10 = wt_onset_2d)

features3 <- long %>%
  group_by(Ferret) %>%
  group_modify(~ compute_features(.x, tp_thr = 1.25, temp_recover_tol = 0.5,
                                  wt_thr_onset = -15, wt_thr_recov = -5, 
                                  k_consec = 2, early_thr = 5)) %>%
  ungroup() %>%
  select(Ferret, wt_onset_1d, wt_onset_2d) %>%
  rename(wt_onset_1d_15 = wt_onset_1d,
         wt_onset_2d_15 = wt_onset_2d)

## check data for Inf & NaN values
sapply(features, function(x) any(is.infinite(x)))
sapply(features, function(x) any(is.nan(x)))
## found in '_ratio' columns

## set Inf & NaN values to NA
features[sapply(features, is.infinite)] <- NA
features[sapply(features, is.nan)] <- NA

features1[sapply(features1, is.infinite)] <- NA
features1[sapply(features1, is.nan)] <- NA

features2[sapply(features2, is.infinite)] <- NA
features2[sapply(features2, is.nan)] <- NA

features3[sapply(features3, is.infinite)] <- NA
features3[sapply(features3, is.nan)] <- NA

#View(features)

## join wt_onset_1d_7p5 and wt_onset_1d_10 to features also wt_onset_1d_15

features <- left_join(features, features1, by = "Ferret")
features <- left_join(features, features2, by = "Ferret")
features <- left_join(features, features3, by = "Ferret")

########################################################################################

# features %>% write.csv(paste("inputs/features_new_", format(Sys.Date(), "%Y-%m-%d"),
#                              ".csv", sep = ""), row.names = FALSE)

########################################################################################

## then join back the full data
## rename or merge with fullMorbidity when finalized
fullMorbidity <- left_join(data, features, by = 'Ferret')

########################################################################################
########################################################################################

### export file

fullMorbidity %>% write.csv(paste("inputs/fullMorbidity_features_new_", 
                                  format(Sys.Date(), "%Y-%m-%d"),
                                  ".csv", sep = ""), row.names = FALSE)

########################################################################################
########################################################################################

### fullMorbidityStats ###

fullMorbidityStats <- fullMorbidity %>%
  group_by(Virus, CDCID, units, Origin, Origin_orig, HA, RBS, PBS, PB2_627, HPAI,
           BnOB_rep, Bn_rep, Int_rep, Lg_rep, Tr_rep) %>%
  dplyr::summarize(AUC_9 = mean(AUC_9_v, na.rm = TRUE),
                   AUC_8 = mean(AUC_8_v, na.rm = TRUE),
                   AUC_6 = mean(AUC_6_v, na.rm = TRUE),
                   AUC_4 = mean(AUC_4_v, na.rm = TRUE),
                   NW_avg = mean(NW_avg, na.rm = TRUE),
                   Lg_avg = mean(Lg_avg, na.rm = TRUE),
                   NT_avg = mean(NT_avg, na.rm = TRUE),
                   Tr_avg = mean(Tr_avg, na.rm = TRUE),
                   BnOB_avg = mean(BnOB_avg, na.rm = TRUE),
                   Bn_avg = mean(Bn_avg, na.rm = TRUE),
                   Int_avg = mean(Int_avg, na.rm = TRUE),
                   d1_avg = mean(d1_inoc, na.rm = TRUE),
                   d3_avg = mean(d3_inoc, na.rm = TRUE),
                   peak_inoc = mean(peak_inoc, na.rm = TRUE),
                   slope13 = mean(slope13_f, na.rm = TRUE),
                   slope35 = mean(slope35_f, na.rm = TRUE),
                   slope57 = mean(slope57_f, na.rm = TRUE),
                   lethal_yes_p = mean(lethal_yes_p, na.rm = TRUE),
                   max_wt_loss_pct = mean(max_wt_loss_pct, na.rm = TRUE),
                   temp_5_avg = mean(temp_5, na.rm = TRUE),
                   wt_loss_slope_max = mean(wt_loss_slope_max, na.rm = TRUE),
                   wt_loss_slope_max5 = mean(wt_loss_slope_max5, na.rm = TRUE),
                   RD_trans_yes_p = mean(RD_trans_yes_p, na.rm = TRUE),
                   DC_trans_yes_p = mean(DC_trans_yes_p, na.rm = TRUE)) %>%
  ungroup() %>% 
  ## keep unique rows only
  distinct() %>%
  mutate(RD_trans_cat_p = ifelse(DC_trans_yes_p == 0 & 
                                   is.na(RD_trans_yes_p), 0, RD_trans_yes_p),
         DC_trans_cat_p = ifelse(RD_trans_yes_p == 1 & 
                                   is.na(DC_trans_yes_p), 1, DC_trans_yes_p)) %>%
  ## create categorical columns
  mutate(wt_loss_cat = case_when(max_wt_loss_pct < -5 ~ 'none',
                                 max_wt_loss_pct < -9.5 ~ 'low',
                                 max_wt_loss_pct < -14.5 ~ 'med',
                                 max_wt_loss_pct < -27.6 ~ 'high'),
         lethal_cat = ifelse(lethal_yes_p < 0.5, "low", "high"),
         RD_trans_cat = case_when(RD_trans_cat_p < 0.3 ~ 'low',
                                  RD_trans_cat_p < 0.67 ~ 'med',
                                  RD_trans_cat_p < 1.1 ~ 'high',
                                  is.na(RD_trans_cat_p) ~ NA_character_),
         DC_trans_cat = case_when(DC_trans_cat_p < 0.3 ~ 'low',
                                  DC_trans_cat_p < 0.67 ~ 'med',
                                  DC_trans_cat_p < 1.1 ~ 'high',
                                  is.na(DC_trans_cat_p) ~ NA_character_)) %>%
  ## Add transmission threshold categories
  mutate(RD_33 = ifelse((RD_trans_yes_p < 0.33 | 
                           is.na(RD_trans_yes_p) & DC_trans_yes_p < 0.1), "no",
                        ifelse(RD_trans_yes_p > 0.33, "yes", "no")),
         RD_50 = ifelse((RD_trans_yes_p < 0.5 | 
                           is.na(RD_trans_yes_p) & DC_trans_yes_p < 0.1), "no",
                        ifelse(RD_trans_yes_p > 0.5, "yes", "no")),
         RD_67 = ifelse((RD_trans_yes_p < 0.67 | 
                           is.na(RD_trans_yes_p) & DC_trans_yes_p < 0.1), "no",
                        ifelse(RD_trans_yes_p > 0.67, "yes", "no")),
         DC_33 = ifelse(DC_trans_yes_p < 0.33 & !is.na(DC_trans_yes_p), "no",
                        ifelse(DC_trans_yes_p > 0.33 | 
                                 (is.na(DC_trans_yes_p) & RD_trans_yes_p == 1), "yes", "no")),
         DC_50 = ifelse(DC_trans_yes_p < 0.5 & !is.na(DC_trans_yes_p), "no",
                        ifelse(DC_trans_yes_p > 0.5 | 
                                 (is.na(DC_trans_yes_p) & RD_trans_yes_p == 1), "yes", "no")),
         DC_67 = ifelse(DC_trans_yes_p < 0.67 & !is.na(DC_trans_yes_p), "no",
                        ifelse(DC_trans_yes_p > 0.67 | 
                                 (is.na(DC_trans_yes_p) & RD_trans_yes_p == 1), "yes", "no"))) %>%
  ## drop redundant columns
  dplyr::select(-c(RD_trans_yes_p, DC_trans_yes_p)) 

## Convert NaN values to NA in the entire dataframe
fullMorbidityStats <- apply(fullMorbidityStats, 2, function(col) ifelse(is.nan(col), NA, col))

########################################################################################
########################################################################################

## Save/export as CSV file w/ current date in the format YYYY-MM-DD
## This will prevent overwriting mistakes, but create 'duplicates'
fullMorbidityStats %>% write.csv(paste("inputs/fullMorbidityStats_clean_", 
                                       format(Sys.Date(), "%Y-%m-%d"), 
                                       ".csv", sep = ""), row.names = FALSE)

########################################################################################
########################################################################################
### End of Code
########################################################################################
########################################################################################
