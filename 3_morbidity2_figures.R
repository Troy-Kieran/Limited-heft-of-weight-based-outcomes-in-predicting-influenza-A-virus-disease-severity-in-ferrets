########################################################################################
########################################################################################
### Morbidity2 Manuscript Figures
### 1 August 2025 - 25 March 2026
### Jessica A. Belser & Troy J. Kieran 
########################################################################################
########################################################################################
### Load packages

library(tidyverse)
library(tidylog)
library(ggstatsplot)
library(janitor)
library(GGally)
library(ggpubr)
library(patchwork)
library(correlation)
library(GGally)
library(FSA)
library(ragg)

########################################################################################
### Import Data
fullMorbidity <- ## data generated from 1_morbidity2_data_cleaning.R
fullMorbidityStats <- ## data generated from 1_morbidity2_data_cleaning.R

########################################################################################

### segregate out lethal and lethal_day records that don't have serial morbidity parameters

fullMorbidity <- fullMorbidity %>%
  mutate(lethal_day = ifelse(is.na(tpk_day) & is.na(wpk_day), NA, lethal_day))  

fullMorbidity <- fullMorbidity %>%
  mutate(lethal = ifelse(is.na(tpk_day) & is.na(wpk_day), NA, lethal)) 

### split datasets by titration matrix when comparing titer v non-titer measures

fullMorbidity_eggs <- fullMorbidity %>%
  filter(units == "EID")

fullMorbidity_cells <- fullMorbidity %>%
  filter(units == "London")

fullMorbidityStats_eggs <- fullMorbidityStats %>%
  filter(units == "EID")

fullMorbidityStats_cells <- fullMorbidityStats %>%
  filter(units == "London")

### filter split datasets for origin skew

fullMorbidity_avian <- fullMorbidity_eggs %>%
  filter(Origin == "avian")
###this removes 134 rows (20%), 547 rows remaining

fullMorbidity_mammal <- fullMorbidity_cells %>%
  filter(Origin == "mammal")
###this removes 132 rows (25%), 389 rows remaining

fullMorbidityStats_avian <- fullMorbidityStats %>%
  filter(Origin == "avian")
###this removes 58 viruses (41%), 82 viruses remaining

fullMorbidityStats_mammal <- fullMorbidityStats %>%
  filter(Origin == "mammal")
###this removes 82 viruses (59%), 58 viruses remaining

###cor tests to support data shown in figure 2b-c
fullMorbidity_notlethal <- fullMorbidity %>%
  filter(lethal == "no")

fullMorbidity_lethal <- fullMorbidity %>%
  filter(lethal == "yes")

### ok run code until here

### content for first results section

########################################################################################

### figure 1 plot
Fig_1A <- fullMorbidity %>%
  ggplot(data = ., aes(tpk_day, fill = Origin))+
  geom_histogram(color = "black", alpha = 0.8, bins = 15) +
  geom_density(alpha = 0.5) +
  ## Opuyama - poisonfrogs - Panama
  scale_fill_manual(values = c("#0086A1", "#451A0E")) +
  scale_x_continuous(breaks = seq(0, 14, 1)) +
  xlab("Day of peak temperature") +
  ylab("Number of ferrets")+
  theme_bw() +
  theme(legend.position = 'bottom')

Fig_1B <- fullMorbidity %>%
  ggplot(data = ., aes(wpk_day, fill = Origin))+
  geom_histogram(color = "black", alpha = 0.8, bins = 15) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#0086A1", "#451A0E")) +
  scale_x_continuous(breaks = seq(0, 14, 1)) +
  xlab("Day of maximum weight loss") +
  ylab("Number of ferrets")+
  theme_bw() +
  theme(legend.position = 'bottom')

Fig_1C <- fullMorbidity %>%
   filter(lethal_day != '0') %>%
  ggplot(data = ., aes(lethal_day, fill = Origin))+
  geom_histogram(color = "black", alpha = 0.8, bins = 12) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#0086A1", "#451A0E")) +
  scale_x_continuous(breaks = seq(2, 15, 1)) +
  xlab("Day of lethality") +
  ylab("Number of ferrets")+
  theme_bw() +
  theme(legend.position = 'none')

Fig_1D <- fullMorbidity %>%
  ggplot(aes(x = abs(max_wt_loss_pct), y = tp_max, fill = Origin)) +
  geom_jitter(aes(color = Origin), alpha = 0.5,
              width = 0.00, height = 0.00, size = 1.5) +
  geom_smooth(method = 'lm', se = TRUE, aes(color = Origin),
              linetype = "longdash", alpha = 0.3) +
  scale_fill_manual(values = c("#0086A1", "#451A0E")) +
  scale_color_manual(values = c("#0086A1", "#451A0E")) +
  xlab('Maximum body weight loss (%)') +
  ylab(expression(paste("Peak temperature rise (", degree, "C)"))) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(color = "Category")
###for this figure legend add cor values, avian is 0.3 [0.21-0.4], p=7.18e-9
###and mammal is 0.052 [-0.065-0.17], p=0.3833
 
Fig_1E <- fullMorbidity_avian %>%
  filter(units == "EID") %>%
  ggplot(aes(x = abs(max_wt_loss_pct), y = Lg_avg)) +
  geom_jitter(color = "#0086A1", alpha = 0.5,
             width = 0.00, height = 0.00, size = 1.5) +
  geom_smooth(method = 'lm', se = TRUE, color = "#0086A1", fill = "#0086A1",
             linetype = "longdash", alpha = 0.3) +
  xlab('Maximum body weight loss (%)') +
  ylab(expression(mean~Lg~titer~(EID[50]/ml))) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(color = "Category")
 
Fig_1F <- fullMorbidity_cells %>%
  filter(units == "London") %>%
  ggplot(aes(x = abs(max_wt_loss_pct), y = NT_avg, fill = Origin)) +
  geom_jitter(aes(color = Origin), alpha = 0.5,
             width = 0.00, height = 0.00, size = 1.5) +
  geom_smooth(method = 'lm', se = TRUE, aes(color = Origin),
             linetype = "longdash", alpha = 0.3) +
  scale_fill_manual(values = c("#0086A1", "#451A0E")) +
  scale_color_manual(values = c("#0086A1", "#451A0E")) +
  xlab('Maximum body weight loss (%)') +
  ylab('mean NT titer (PFU/ml)') +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(color = "Category")
 ### for figure legend, avian cor is r=0.196 [0.087-0.299], p=0.0004769
 ###and mammal cor is -0.032 [-0.16-0.096], p=0.6287
 
 
Fig_1A + Fig_1B + Fig_1C + Fig_1D + Fig_1E + Fig_1F + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = 'A') -> Figure_1


## export TIFF
ragg::agg_tiff(filename = "Figure_1.tiff",
               width = 6, height = 6.5, units = 'in', res = 600,
               scaling = 0.8, compression = "lzw")
plot(Figure_1)
invisible(dev.off())
 
########################################################################################
 
### content for results section 2

###redoing the density plot from DMM for weight onset day
 
Fig_2A <- fullMorbidity %>%
   ggplot() +
   geom_density(aes(x = wt_onset_1d, fill = "Below 5%"), alpha = 0.5) +
   geom_density(aes(x = wt_onset_1d_7p5, fill = "Below 7.5%"), alpha = 0.5) +
   geom_density(aes(x = wt_onset_1d_10, fill = "Below 10%"), alpha = 0.5) +
   geom_density(aes(x = wt_onset_1d_15, fill = "Below 15%"), alpha = 0.5) +
   scale_fill_manual(breaks = c("Below 5%", "Below 7.5%", "Below 10%", "Below 15%"),
                     ## Daratus - poisonfrogs - Colombia
                     values = c("#1B445B", "#4AA3B8", "#2AB8A2", "#96D02B")) +
   scale_color_manual(values = c("#1B445B", "#4AA3B8", "#2AB8A2", "#96D02B")) +
   xlab("Day of first weight loss below cutoff") +
   ylab("Density") +
   scale_x_continuous(breaks = seq(0, 14, 1)) +
   theme_bw() +
   theme(legend.position = c(0.8, 0.8),
         legend.title = element_blank())

Fig_2B <- fullMorbidity %>%
  filter(!is.na(lethal)) %>%
  filter(expt == "DC" | expt == "RD" | expt == "path") %>%
  ggplot(aes(x = abs(max_wt_loss_pct), y = wt_onset_1d_7p5, fill = lethal)) +
  geom_jitter(aes(color = lethal, shape = lethal), alpha = 0.5,
              width = 0.00, height = 0.08, size = 1.5) +
  geom_smooth(method = 'lm', se = TRUE, aes(color = lethal),
              linetype = "longdash", alpha = 0.4, show.legend = FALSE) +
  ## Pvaillantii - poisonfrogs - Colombia
  scale_fill_manual(values = c("#979136", "#5E2421"), guide = "none") +
  scale_color_manual(values = c(
      colorspace::lighten("#979136", 0.35),
      colorspace::darken("#5E2421",  0.9)), labels = c("not lethal", "lethal")) +
  scale_shape_manual(labels = c("not lethal", "lethal"), values = c(16, 17)) +
  scale_x_continuous(limits = c(7.5, 27), breaks = seq(7, 27, 2)) +
  scale_y_continuous(limits = c(1, 14), breaks = seq(1, 14, 2)) +
  xlab('Maximum body weight loss (%)') +
  ylab('Day of first weight loss >7.5%') +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_blank())

Fig_2C <- fullMorbidity %>%
  filter(!is.na(lethal)) %>%
  filter(expt == "DC" | expt == "RD" | expt == "path") %>%
  ggplot(aes(x = abs(max_wt_loss_pct), y = wt_days_7p5_d7, fill = lethal)) +
  geom_jitter(aes(color = lethal, shape = lethal), alpha = 0.5,
              width = 0.00, height = 0.08, size = 1.5) +
  geom_smooth(method = 'lm', se = TRUE, aes(color = lethal),
              linetype = "longdash", alpha = 0.4) +
  scale_fill_manual(values = c("#979136", "#5E2421")) +
  scale_color_manual(values = c(
    colorspace::lighten("#979136", 0.35),
    colorspace::darken("#5E2421",  0.9))) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_continuous(limits = c(7.5, 27), breaks = seq(7, 27, 2)) +
  scale_y_continuous(limits = c(1, 6.5), breaks = seq(1, 7, 1)) +
  xlab('Maximum body weight loss (%)') +
  ylab('Number of days weight loss >7.5% 1-7 p.i.') +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(color = "Category")
  

Fig_2A + (Fig_2B + Fig_2C) + 
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'A') -> Figure_2

## export TIFF
ragg::agg_tiff(filename = "Figure_2.tiff",
               width = 6, height = 6, units = 'in', res = 600,
               scaling = 0.84, compression = "lzw")
plot(Figure_2)
invisible(dev.off())

########################################################################################

###time for 3rd results section

Fig_3A <- fullMorbidity %>%
  filter(expt == "DC" | expt == "RD" | expt == "path") %>%
  filter(!is.na(lethal)) %>%
  ggplot(aes(x = lethal, y = max_wt_loss_pct, color = lethal, fill = lethal)) +
  geom_violin(alpha = 0.5) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), 
              aes(group = lethal, shape = lethal), size = 1.5, alpha = 0.5) + 
  ggpubr::geom_pwc(method = 'tukey_hsd', label = 'p.adj.signif', 
                   step.increase = 0.07, hide.ns = TRUE, vjust = 0.5) +
  labs(x = "Treatment Group", y = "Weight loss (% maximum)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  scale_color_manual(values = c(
    colorspace::lighten("#979136", 0.35),
    colorspace::darken("#5E2421",  0.9))) +
  scale_fill_manual(values = c("#979136", "#5E2421")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_discrete(labels = c("no" = "not lethal", "yes" = "lethal"))
 
Fig_3B <- fullMorbidity %>%
  filter(expt == "DC" | expt == "RD" | expt == "path") %>%
  filter(!is.na(lethal)) %>%
  ggplot(aes(x = lethal, y = wt_loss_slope_max, color = lethal, fill = lethal)) +
  geom_violin(alpha = 0.5) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), 
              aes(group = lethal, shape = lethal), size = 1.5, alpha = 0.5) + 
  ggpubr::geom_pwc(method = 'tukey_hsd', label = 'p.adj.signif', 
                   step.increase = 0.07, hide.ns = TRUE, vjust = 0.5) +
  labs(x = "Treatment Group", y = "Weight loss (slope)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  scale_color_manual(values = c(
    colorspace::lighten("#979136", 0.35),
    colorspace::darken("#5E2421",  0.9))) +
  scale_fill_manual(values = c("#979136", "#5E2421")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_discrete(labels = c("no" = "not lethal", "yes" = "lethal"))
 
Fig_3C <- fullMorbidity %>%
   filter(!is.na(lethal)) %>%
   filter(expt == "DC" | expt == "RD" | expt == "path") %>%
   ggplot(aes(x = abs(max_wt_loss_pct), y = wt_loss_slope_max, fill = lethal)) +
   geom_jitter(aes(color = lethal, shape = lethal), alpha = 0.5,
               width = 0.00, height = 0.00, size = 1.5) +
   geom_smooth(method = 'lm', se = TRUE, aes(color = lethal),
               linetype = "longdash", alpha = 0.5) +
   scale_fill_manual(values = c("#979136", "#5E2421")) +
   scale_color_manual(values = c(
     colorspace::lighten("#979136", 0.35),
     colorspace::darken("#5E2421",  0.9))) +
   scale_shape_manual(values = c(16, 17)) +
   xlab('Maximum body weight loss (%)') +
   ylab('Slope of peak weight loss') +
   theme_bw() +
   theme(legend.position = 'none') +
   labs(color = "Category")
 
(Fig_3A + Fig_3B) + Fig_3C + 
  plot_layout(ncol = 3, widths = c(1, 1, 2)) +
  plot_annotation(tag_levels = "A") -> Figure_3
 
## export TIFF
ragg::agg_tiff(filename = "Figure_3.tiff",
              width = 6, height = 2.5, units = 'in', res = 600,
              scaling = 0.8, compression = "lzw")
plot(Figure_3)
invisible(dev.off())
 
########################################################################################
 
###graphs for proposed figure 4 for svd equivalents
###first lets do weight
 
weights_long <- fullMorbidity %>%
  filter(expt == "DC" | expt == "RD" | expt == "path") %>%
  pivot_longer(cols = c(matches("^wt_\\d+_n")),
               names_to = "Day",
               values_to = "Weight") %>%
  group_by(Ferret) %>%
  # mutate(mean_wt = mean(Weight, na.rm = TRUE),
  #        SD = sd(Weight, na.rm = TRUE),
  #        n = n(),
  #        SE = SD / sqrt(n),
  #        ci_lower = mean_wt - 1.96 * SE,
  #        ci_upper = mean_wt + 1.96 * SE) %>%
  ungroup()
 
###order the day observations correctly
weights_long$Day <- factor(weights_long$Day, levels = c("wt_0_n", "wt_1_n", "wt_2_n", "wt_3_n",
                                                        "wt_4_n", "wt_5_n", "wt_6_n", "wt_7_n", "wt_8_n",
                                                        "wt_9_n", "wt_10_n", "wt_11_n", "wt_12_n", 
                                                        "wt_13_n", "wt_14_n"))
 
Fig_4A <- weights_long %>%
  filter(!is.na(lethal)) %>%
  subset(Day %in% c("wt_0_n", "wt_1_n", "wt_2_n", "wt_3_n",
                    "wt_4_n", "wt_5_n", "wt_6_n", "wt_7_n", "wt_8_n",
                    "wt_9_n", "wt_10_n", "wt_11_n", "wt_12_n")) %>%
  ggplot(aes(x = Day, y = Weight, color = lethal, group = Ferret)) +
  geom_line(linewidth = 0.8, alpha = 0.04) +
  stat_summary(fun = mean, geom = "line", size = 1.5, aes(group = lethal), linetype = "solid") +
  labs(x = "Day post-inoculation", y = "% weight change") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = 'bottom') +
  scale_y_continuous(limits = c(-30, 5, breaks = c(-30, -25, -20, -15, -10, -5, 0, 5))) +
  scale_color_manual(values = c(
    colorspace::lighten("#979136", 0.35),
    colorspace::darken("#5E2421",  0.9))) +
  scale_fill_manual(values = c("#979136", "#5E2421")) +
  labs(x = "Day post-inoculation") +
  scale_x_discrete(labels = c("wt_0_n" = "0", "wt_1_n" = "1", "wt_2_n" = "2", "wt_3_n" = "3", 
                              "wt_4_n" = "4", "wt_5_n" = "5", "wt_6_n" = "6", "wt_7_n" = "7", 
                              "wt_8_n" = "8", "wt_9_n" = "9", "wt_10_n" = "10", "wt_11_n" = "11",
                              "wt_12_n" = "12", "wt_13_n" = "13", "wt_14_n" = "14"))
 
### root mean squares stats

Fig_4B <- fullMorbidity %>%
  filter(!is.na(lethal)) %>%
  filter(expt == "DC" | expt == "RD" | expt == "path") %>%
  ggplot(aes(x = lethal, shape = lethal, y = wt_rmssd_rate, 
             color = lethal, fill = lethal)) +
  geom_violin(alpha = 0.5) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), 
              aes(group = lethal),
              size = 1.5, alpha = 0.5) + 
  ggpubr::geom_pwc(method = 'wilcox_test', label = 'p.adj.signif', 
                   step.increase = 0.07, hide.ns = TRUE, vjust = 0.5) +
  labs(x = "Treatment Group", y = "RMSSD (weight)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  scale_color_manual(values = c(
    colorspace::lighten("#979136", 0.35),
    colorspace::darken("#5E2421",  0.9))) +
  scale_fill_manual(values = c("#979136", "#5E2421")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_discrete(labels = c("no" = "not lethal", "yes" = "lethal"))

## ACF1

Fig_4C <- fullMorbidity %>%
  filter(!is.na(lethal)) %>%
  filter(expt == "DC" | expt == "RD" | expt == "path") %>%
  ggplot(aes(x = lethal, shape = lethal, y = wt_acf1, 
             color = lethal, fill = lethal)) +
  geom_violin(alpha = 0.5) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), 
              aes(group = lethal),
              size = 1.5, alpha = 0.5) + 
  ggpubr::geom_pwc(method = 'wilcox_test', label = 'p.adj.signif', 
                   step.increase = 0.07, hide.ns = TRUE, vjust = 0.5) +
  labs(x = "Treatment Group", y = "ACF (weight)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  scale_color_manual(values = c(
    colorspace::lighten("#979136", 0.35),
    colorspace::darken("#5E2421",  0.9))) +
  scale_fill_manual(values = c("#979136", "#5E2421")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_discrete(labels = c("no" = "not lethal", "yes" = "lethal"))

### temp
temps_long <- fullMorbidity %>%
  filter(expt == "DC" | expt == "RD" | expt == "path") %>%
  pivot_longer(cols = c(matches("^tp_\\d+_n")),
               names_to = "Day",
               values_to = "Temp") %>%
  group_by(Ferret) %>%
  ungroup()

###order the day observations correctly
temps_long$Day <- factor(temps_long$Day, levels = c("tp_0_n", "tp_1_n", "tp_2_n", "tp_3_n",
                                                    "tp_4_n", "tp_5_n", "tp_6_n", "tp_7_n", "tp_8_n",
                                                    "tp_9_n", "tp_10_n", "tp_11_n", "tp_12_n",
                                                        "tp_13_n", "tp_14_n"))

Fig_4D <- temps_long %>%
  filter(!is.na(lethal)) %>%
  subset(Day %in% c("tp_0_n", "tp_1_n", "tp_2_n", "tp_3_n",
                    "tp_4_n", "tp_5_n", "tp_6_n", "tp_7_n", "tp_8_n",
                    "tp_9_n", "tp_10_n", "tp_11_n", "tp_12_n")) %>%
  ggplot(aes(x = Day, y = Temp, color = lethal, group = Ferret)) +
  geom_line(linewidth = 0.8, alpha = 0.04) +
  stat_summary(fun = mean, geom = "line", size = 1.5, aes(group = lethal), linetype = "solid") +
  labs(x = "Day post-inoculation", y = expression("" * degree * "C change")) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = 'none') +
  scale_y_continuous(limits = c(-5, 5, breaks = c(-5, 0, 5))) +
  scale_color_manual(values = c(
    colorspace::lighten("#979136", 0.35),
    colorspace::darken("#5E2421",  0.9))) +
  scale_fill_manual(values = c("#979136", "#5E2421")) +
  labs(x = "Day post-inoculation") +
  scale_x_discrete(labels = c("tp_0_n" = "0", "tp_1_n" = "1", "tp_2_n" = "2", "tp_3_n" = "3", 
                              "tp_4_n" = "4", "tp_5_n" = "5", "tp_6_n" = "6", "tp_7_n" = "7", 
                              "tp_8_n" = "8", "tp_9_n" = "9", "tp_10_n" = "10", "tp_11_n" = "11",
                              "tp_12_n" = "12", "tp_13_n" = "13", "tp_14_n" = "14"))
 

### root mean squares stats

Fig_4E <- fullMorbidity %>%
  filter(!is.na(lethal)) %>%
  filter(expt == "DC" | expt == "RD" | expt == "path") %>%
  ggplot(aes(x = lethal, shape = lethal, y = temp_rmssd_rate, 
             color = lethal, fill = lethal)) +
  geom_violin(alpha = 0.5) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), 
              aes(group = lethal),
              size = 1.5, alpha = 0.5) + 
  ggpubr::geom_pwc(method = 'wilcox_test', label = 'p.adj.signif', 
                   step.increase = 0.07, hide.ns = TRUE, vjust = 0.5) +
  labs(x = "Treatment Group", y = "RMSSD (temperature)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  scale_color_manual(values = c(
    colorspace::lighten("#979136", 0.35),
    colorspace::darken("#5E2421",  0.9))) +
  scale_fill_manual(values = c("#979136", "#5E2421")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_discrete(labels = c("no" = "not lethal", "yes" = "lethal"))

## and now ACF1

Fig_4F <- fullMorbidity %>%
  filter(!is.na(lethal)) %>%
  filter(expt == "DC" | expt == "RD" | expt == "path") %>%
  ggplot(aes(x = lethal, shape = lethal, y = temp_acf1, 
             color = lethal, fill = lethal)) +
  geom_violin(alpha = 0.5) + 
  geom_jitter(position = position_jitter(width = 0.1, height = 0), 
              aes(group = lethal),
              size = 1.5, alpha = 0.5) + 
  ggpubr::geom_pwc(method = 'wilcox_test', label = 'p.adj.signif', 
                   step.increase = 0.07, hide.ns = TRUE, vjust = 0.5) +
  labs(x = "Treatment Group", y = "ACF (temperature)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  scale_color_manual(values = c(
    colorspace::lighten("#979136", 0.35),
    colorspace::darken("#5E2421",  0.9))) +
  scale_fill_manual(values = c("#979136", "#5E2421")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_discrete(labels = c("no" = "not lethal", "yes" = "lethal"))


Fig_4A + Fig_4B + Fig_4C + Fig_4D + Fig_4E + Fig_4F +
  plot_layout(ncol = 3, widths = c(2,1,1)) +
  plot_annotation(tag_levels = "A") -> Figure_4

## export TIFF
ragg::agg_tiff(filename = "Figure_4.tiff",
               width = 6, height = 4, units = 'in', res = 600,
               scaling = 0.8, compression = "lzw")
plot(Figure_4)
invisible(dev.off())

########################################################################################

### Figure 5

## summarize one CV metric for multiple variables
summarize_rcvm <- function(data, vars) {
  data %>%
    filter(NW_typical == "yes") %>%
    pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "val") %>%
    group_by(wt_loss_high, variable) %>%
    reframe(
      med = median(val, na.rm = TRUE),
      mad_val = mad(val, constant = 1, na.rm = TRUE),
      rcvm = ifelse(med == 0 | is.na(med), NA_real_, abs(1.4826 * (mad_val / med) * 100))) %>% drop_na()}

#summarize_rcvm(fullMorbidity, vars = c("wt_loss_slope_max", "wt_rmssd_rate", "wt_acf1"))

## use function to plot, swapping out multiple vars

## Fig5A
summarize_rcvm(data = fullMorbidity, 
               vars = c("max_wt_loss_pct", "wt_loss_slope_max", 
                        "wt_rmssd_rate", "wt_acf1")) %>% 
  mutate(variable = factor(variable, levels = c("max_wt_loss_pct", "wt_loss_slope_max", 
                                                "wt_rmssd_rate", "wt_acf1"))) %>%
  ggplot(aes(x = wt_loss_high, y = rcvm, fill = variable)) +
  ## Amacero colors from the poisonfrogs R package - Peru
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 20, ymax = 30,
           fill = "#C6E74B", alpha = 0.3) +
  geom_hline(yintercept = 25, linetype = 'dashed',linewidth = 1.2, color = '#C6E74B') +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(labels = c(max_wt_loss_pct = "Max Loss %", 
                               wt_loss_slope_max = "Max Slope", 
                               wt_rmssd_rate = "RMSSD", 
                               wt_acf1 = "ACF1"),
                    values = c("max_wt_loss_pct" = "#14374A", 
                               "wt_loss_slope_max" = "#3981BE", 
                               "wt_rmssd_rate" = "#64B936", 
                               "wt_acf1" = "#7D2A1D")) +
  labs(x = "Weight Loss High", y = "Coefficient of Variation - Weight", fill = element_blank()) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'bottom',
        panel.spacing = unit(0, "pt"),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = alpha("#C6E74B", 0.3), color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1)) -> cv_weight

## Fig5B
summarize_rcvm(data = fullMorbidity, 
               vars = c("temp_5", "temp_rmssd_rate", "temp_acf1")) %>%
  mutate(variable = factor(variable, levels = c("temp_5", "temp_rmssd_rate", "temp_acf1"))) %>%
  ggplot(aes(x = wt_loss_high, y = rcvm, fill = variable)) +
  ## Amacero colors from the poisonfrogs R package - Peru
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 20, ymax = 30,
           fill = "#C6E74B", alpha = 0.3) +
  geom_hline(yintercept = 25, linetype = 'dashed',linewidth = 1.2, color = '#C6E74B') +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(labels = c(temp_5 = "Temp Day 1-5", 
                               temp_rmssd_rate = "RMSSD", 
                               temp_acf1 = "ACF1"),
                    values = c("temp_5" = "#14374A", 
                               "temp_rmssd_rate" = "#64B936", 
                               "temp_acf1" = "#7D2A1D")) +
  labs(x = "Weight Loss High", y = "Coefficient of Variation - Temp", fill = element_blank()) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'bottom',
        panel.spacing = unit(0, "pt"),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = alpha("#C6E74B", 0.3), color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1)) -> cv_temp


Figure_5 <- cv_weight / cv_temp + 
  plot_annotation(tag_levels = "A")

## export TIFF
ragg::agg_tiff(filename = "Figure_5.tiff",
               width = 4, height = 6, units = 'in', res = 600,
               scaling = 0.85, compression = "lzw")
plot(Figure_5)
invisible(dev.off())

########################################################################################
########################################################################################
### End of Code
########################################################################################
########################################################################################
