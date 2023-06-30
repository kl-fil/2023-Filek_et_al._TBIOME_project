# Dataviz script for TurtleBiome Holobiont 2023 manuscript
## Processing and visualizing data related to sequencing results for
## endozoic and tank water samples (ITS2 region of V4 18S rRNA gene)
## Borna V. Vuković performed initial analyses
##--------------------load packages--------------------

library(tidyverse) #v2.0.0
library(ggplot2) #v3.4.2
library(qiime2R) #qiime to R import qza v0.99.6
library(patchwork) #collect multiple ggplots together v1.1.2
library(ggrepel) #repel text/labels v0.9.3
library(ggforce) #v0.4.1
library(ggtext) #v0.1.2
library(vegan) #v2.6-4
library(ggdendro) #v0.1.23
library(pairwiseAdonis) # ???
library(ggpubr) #0.6.0
library(pheatmap) #1.0.12
library(ggdist) #v3.3.0

#set global theme for ggplot
theme_set(theme_light(base_size = 11, base_family = "Arial")) #<< font!

#vectors for colors per sampling site
c("#93624d", "#ff9999", "#89d6dc")
c("#462f25","#750000", "#16474b")

# "#<<" comment indicates the line of code that needs to be changed when using 
# a different dataset for the same code

##---- make directory for visualizations ----
dir.create("r_output/ITS2_endo/", recursive = TRUE)

##---- Load metadata and data-----
# load metadata in tsv format
metadata_ITS <- readr::read_tsv("master_metadata/16SendoV34.tsv")
metadata_ITS <- metadata_ITS %>% mutate(SampleID = gsub("16S", "ITS", SampleID)) %>%
  filter(SampleSite %in% c("CLOACA", "TANK WATER", "NEGCTRL"))

##---- Alpha diversity ----

shannonITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results/shannon_vector.qza")
f_pdITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results/faith_pd_vector.qza")
evennessITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results/evenness_vector.qza")
observed_asvITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results/observed_features_vector.qza")

alpha_values_meta_ITS <- 
  metadata_ITS %>%
  left_join(rownames_to_column(shannonITS$data, "SampleID")) %>%
  left_join(rownames_to_column(f_pdITS$data, "SampleID")) %>%
  left_join(rownames_to_column(evennessITS$data, "SampleID")) %>%
  left_join(rownames_to_column(observed_asvITS$data, "SampleID")) %>%
  drop_na(observed_features)

alpha_values_meta_ITS$SampleSite <- ordered(alpha_values_meta_ITS$SampleSite,
                                            levels = c("CLOACA", "TANK WATER"))

##----- Set theme for alpha div plots -----
theme_alpha <- theme(
  axis.text.x = element_text(angle = 0, vjust = 0.5, 
                             hjust = 0.5, size = 11, 
                             color = "grey30"), 
  strip.text.x = element_text(size = 11, color = "gray30", face = "bold"),
  strip.background = element_rect(color = "gray20"),
  panel.grid = element_line(colour = "gray80"),
  panel.grid.major.x = element_blank(),
  panel.border = element_rect(color = "gray30"),
  axis.ticks.x = element_line(colour = "grey30"), # << visibility
  axis.ticks.y = element_line(colour = "grey30"),
  axis.title.x = element_text(size = 12, color = "black"),
  #axis.title.x = element_blank(), #<< visibility
  legend.text = element_markdown(size = 10, color = "grey30"),
  axis.title.y = element_text(size = 12, color = "black"),
  axis.text.y = element_text(size = 10, color = "grey30"),
  legend.title = element_text(size = 10, color = "black"),
  legend.spacing.y = unit(0.05, 'cm'),
  legend.key.size = unit(0.5, 'cm'),
  plot.tag = element_text(size = 12)
)

theme_alpha_pub <- theme(
  axis.text.x = element_text(angle = 0, vjust = 0.5, 
                             hjust = 0.5, size = 8, 
                             color = "black"), 
  strip.text.x = element_text(size = 8, color = "black", face = "bold"),
  strip.background = element_rect(color = "black"),
  panel.grid = element_line(colour = "gray90"),
  panel.grid.major.x = element_blank(),
  panel.border = element_rect(color = "black"),
  axis.ticks.x = element_line(colour = "black"), # << visibility
  axis.ticks.y = element_line(colour = "black"),
  axis.title.x = element_text(size = 10, color = "black"),
  #axis.title.x = element_blank(), #<< visibility
  legend.text = element_markdown(size = 8, color = "black"),
  axis.title.y = element_text(size = 10, color = "black"),
  axis.text.y = element_text(size = 8, color = "black"),
  legend.title = element_text(size = 8, color = "black"),
  legend.spacing.y = unit(0.05, 'cm'),
  legend.key.size = unit(0.5, 'cm'),
  plot.tag = element_text(size = 10)
)

theme_pcoa1_pub <- theme(
  panel.border = element_rect(fill = "transparent", color = "black"),
  plot.subtitle = element_text(size = 11),
  axis.title.x = element_text(size = 11, color = "black"),
  axis.title.y = element_text(size = 11, color = "black"),
  axis.text.x = element_text(size = 9, color = "black"),
  axis.text.y = element_text(size = 9, color = "black"),
  legend.title = element_text(size = 11, color = "black"),
  legend.text = element_markdown(size = 10, color = "black"),
  legend.spacing.y = unit(0.05, 'cm'),
  legend.key.size = unit(0.5, 'cm'),
  panel.grid.minor.x = element_line(linewidth = 0, color = "transparent"),
  panel.grid.minor.y = element_line(linewidth = 0, color = "transparent"),
  panel.grid.major.x = element_line(linewidth = 0, color = "transparent"),
  panel.grid.major.y = element_line(linewidth = 0, color = "transparent"),
  axis.line.x.bottom = element_line(linewidth = 0.3, color = "black"),
  axis.line.y.left = element_line(linewidth = 0.3, color = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.tag = element_text(size = 11)
)

##----- Alpha div plots -----

observed_alpha_rain_ITS <-
  alpha_values_meta_ITS %>%
  drop_na(SampleSite) %>%
  ggplot(aes(x = SampleSite, y = observed_features)) +
  ggdist::stat_halfeye(aes(fill = SampleSite),
                       adjust = .5, 
                       width = .6, 
                       .width = 0, 
                       justification = -.3, 
                       point_colour = NA
  ) + 
  geom_boxplot(aes(fill = SampleSite),
               alpha = 0.5,
               width = .25, 
               outlier.shape = NA,
               color = c("#462f25", "#16474b")
  ) +
  geom_point(aes(fill = SampleSite),
             shape = 21,
             color = "black",
             size = 1.2,
             alpha = 1,
             position = position_jitter(
               seed = 1, width = .1
             )
  ) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              "TANK WATER" = "Tank \nwater\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Observed features (ASVs)", #<<
    tag = "a)" #<<
  ) + 
  theme_alpha_pub +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

faith_alpha_rain_ITS <-
  alpha_values_meta_ITS %>%
  drop_na(SampleSite) %>%
  ggplot(aes(x = SampleSite, y = faith_pd)) +
  ggdist::stat_halfeye(aes(fill = SampleSite),
                       adjust = .5, 
                       width = .6, 
                       .width = 0, 
                       justification = -.3, 
                       point_colour = NA
  ) + 
  geom_boxplot(aes(fill = SampleSite),
               alpha = 0.5,
               width = .25, 
               outlier.shape = NA,
               color = c("#462f25", "#16474b")
  ) +
  geom_point(aes(fill = SampleSite),
             shape = 21,
             color = "black",
             size = 1.2,
             alpha = 1,
             position = position_jitter(
               seed = 1, width = .1
             )
  ) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              "TANK WATER" = "Tank \nwater\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Faith's PD", #<<
    tag = "b)" #<<
  ) + 
  theme_alpha_pub +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

shannon_alpha_rain_ITS <-
  alpha_values_meta_ITS %>%
  drop_na(SampleSite) %>%
  ggplot(aes(x = SampleSite, y = shannon_entropy)) +
  ggdist::stat_halfeye(aes(fill = SampleSite),
                       adjust = .5, 
                       width = .6, 
                       .width = 0, 
                       justification = -.3, 
                       point_colour = NA
  ) + 
  geom_boxplot(aes(fill = SampleSite),
               alpha = 0.5,
               width = .25, 
               outlier.shape = NA,
               color = c("#462f25", "#16474b")
  ) +
  geom_point(aes(fill = SampleSite),
             shape = 21,
             color = "black",
             size = 1.2,
             alpha = 1,
             position = position_jitter(
               seed = 1, width = .1
             )
  ) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              "TANK WATER" = "Tank \nwater\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Shannon's entropy", #<<
    tag = "c)" #<<
  ) + 
  theme_alpha_pub +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

evenness_alpha_rain_ITS <-
  alpha_values_meta_ITS %>%
  drop_na(SampleSite) %>%
  ggplot(aes(x = SampleSite, y = pielou_evenness)) +
  ggdist::stat_halfeye(aes(fill = SampleSite),
                       adjust = .5, 
                       width = .6, 
                       .width = 0, 
                       justification = -.3, 
                       point_colour = NA
  ) + 
  geom_boxplot(aes(fill = SampleSite),
               alpha = 0.5,
               width = .25, 
               outlier.shape = NA,
               color = c("#462f25", "#16474b")
  ) +
  geom_point(aes(fill = SampleSite),
             shape = 21,
             color = "black",
             size = 1.2,
             alpha = 1,
             position = position_jitter(
               seed = 1, width = .1
             )
  ) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              "TANK WATER" = "Tank \nwater\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Pielou's evenness", #<<
    tag = "d)" #<<
  ) + 
  theme_alpha_pub +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

rain_plots_ITS <- (observed_alpha_rain_ITS+faith_alpha_rain_ITS)/(shannon_alpha_rain_ITS+evenness_alpha_rain_ITS)

ggsave(
    "r_output/ITS2_endo/alpha-div-collected.pdf",
  plot = rain_plots_ITS,
  device = cairo_pdf,
  height = 150,
  width = 160,
  units = "mm"
)

##----- Kruskal-Wallis test -----

observed_summarystats_ITS <- group_by(alpha_values_meta_ITS, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(observed_features, na.rm = TRUE),
    sd = sd(observed_features, na.rm = TRUE),
    median = median(observed_features, na.rm = TRUE),
    IQR = IQR(observed_features, na.rm = TRUE)
  )

faith_summarystats_ITS <- group_by(alpha_values_meta_ITS, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(faith_pd, na.rm = TRUE),
    sd = sd(faith_pd, na.rm = TRUE),
    median = median(faith_pd, na.rm = TRUE),
    IQR = IQR(faith_pd, na.rm = TRUE)
  )

shannon_summarystats_ITS <- group_by(alpha_values_meta_ITS, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(shannon_entropy, na.rm = TRUE),
    sd = sd(shannon_entropy, na.rm = TRUE),
    median = median(shannon_entropy, na.rm = TRUE),
    IQR = IQR(shannon_entropy, na.rm = TRUE)
  )

evenness_summarystats_ITS <- group_by(alpha_values_meta_ITS, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(pielou_evenness, na.rm = TRUE),
    sd = sd(pielou_evenness, na.rm = TRUE),
    median = median(pielou_evenness, na.rm = TRUE),
    IQR = IQR(pielou_evenness, na.rm = TRUE)
  )

summary_stats_alpha_ITS <- list()
summary_stats_alpha_ITS[["Observed features"]] <- observed_summarystats_ITS
summary_stats_alpha_ITS[["Faith PD"]] <- faith_summarystats_ITS
summary_stats_alpha_ITS[["Shannon entropy"]] <- shannon_summarystats_ITS
summary_stats_alpha_ITS[["Pielou evenness"]] <- evenness_summarystats_ITS

write_tsv(do.call(cbind, summary_stats_alpha_ITS), "r_output/ITS2_endo/alpha_diversity_summary_stats.tsv")
capture.output(summary_stats_alpha_ITS, file = "r_output/ITS2_endo/alpha_diversity_summary_stats.txt")

##---- Kruskall-Wallis test & pairwise comparisons (Wilcox) ----
alpha_values_meta_ITS_noNA <- alpha_values_meta_ITS %>% drop_na(SampleSite)

krustats_obs_ITS <- list()
krustats_obs_ITS[["Kruskall-Wallis test"]] <- kruskal.test(observed_features ~ SampleSite, data = alpha_values_meta_ITS_noNA)
krustats_obs_ITS[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_meta_ITS_noNA$observed_features, alpha_values_meta_ITS_noNA$SampleSite,
                                                                           p.adjust.method = "BH")
capture.output(krustats_obs_ITS, file = "r_output/ITS2_endo/alpha_diversity-Kruskall-observed.txt")

krustats_faith_ITS <- list()
krustats_faith_ITS[["Kruskall-Wallis test"]] <- kruskal.test(faith_pd ~ SampleSite, data = alpha_values_meta_ITS_noNA)
krustats_faith_ITS[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_meta_ITS_noNA$faith_pd, alpha_values_meta_ITS_noNA$SampleSite,
                                                                             p.adjust.method = "BH")
capture.output(krustats_faith_ITS, file = "r_output/ITS2_endo/alpha_diversity-Kruskall-faith_pd.txt")

krustats_shan_ITS <- list()
krustats_shan_ITS[["Kruskall-Wallis test"]] <- kruskal.test(shannon_entropy ~ SampleSite, data = alpha_values_meta_ITS_noNA)
krustats_shan_ITS[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_meta_ITS_noNA$shannon_entropy, alpha_values_meta_ITS_noNA$SampleSite,
                                                                            p.adjust.method = "BH")
capture.output(krustats_shan_ITS, file = "r_output/ITS2_endo/alpha_diversity-Kruskall-shannon.txt")

krustats_even_ITS <- list()
krustats_even_ITS[["Kruskall-Wallis test"]] <- kruskal.test(pielou_evenness ~ SampleSite, data = alpha_values_meta_ITS_noNA)
krustats_even_ITS[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_meta_ITS_noNA$pielou_evenness, alpha_values_meta_ITS_noNA$SampleSite,
                                                                            p.adjust.method = "BH")
capture.output(krustats_even_ITS, file = "r_output/ITS2_endo/alpha_diversity-Kruskall-evenness.txt")

kruskall_cloaca_ITS <- list()
kruskall_cloaca_ITS[["HospStatus_faith"]] <- kruskal.test(faith_pd ~ HospStatus, 
                                                          data = alpha_values_meta_ITS_noNA %>%
                                                            filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["TurtleSex_faith"]] <- kruskal.test(faith_pd ~ TurtleSex, 
                                                         data = alpha_values_meta_ITS_noNA %>% 
                                                           filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["DetSex_faith"]] <- kruskal.test(faith_pd ~ DetSex, 
                                                      data = alpha_values_meta_ITS_noNA %>%
                                                        filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["AgeRange_faith"]] <- kruskal.test(faith_pd ~ AgeRange, 
                                                        data = alpha_values_meta_ITS_noNA %>%
                                                          filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["Age2_faith"]] <- kruskal.test(faith_pd ~ Age2, 
                                                    data = alpha_values_meta_ITS_noNA %>%
                                                      filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["HospStatus_observed"]] <- kruskal.test(observed_features ~ HospStatus, 
                                                             data = alpha_values_meta_ITS_noNA %>%
                                                               filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["TurtleSex_observed"]] <- kruskal.test(observed_features ~ TurtleSex, 
                                                            data = alpha_values_meta_ITS_noNA %>% 
                                                              filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["DetSex_observed"]] <- kruskal.test(observed_features ~ DetSex, 
                                                         data = alpha_values_meta_ITS_noNA %>%
                                                           filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["AgeRange_observed"]] <- kruskal.test(observed_features ~ AgeRange, 
                                                           data = alpha_values_meta_ITS_noNA %>%
                                                             filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["Age2_observed"]] <- kruskal.test(observed_features ~ Age2, 
                                                       data = alpha_values_meta_ITS_noNA %>%
                                                         filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["HospStatus_shannon"]] <- kruskal.test(shannon_entropy ~ HospStatus, 
                                                            data = alpha_values_meta_ITS_noNA %>%
                                                              filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["TurtleSex_shannon"]] <- kruskal.test(shannon_entropy ~ TurtleSex, 
                                                           data = alpha_values_meta_ITS_noNA %>% 
                                                             filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["DetSex_shannon"]] <- kruskal.test(shannon_entropy ~ DetSex, 
                                                        data = alpha_values_meta_ITS_noNA %>%
                                                          filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["AgeRange_shannon"]] <- kruskal.test(shannon_entropy ~ AgeRange, 
                                                          data = alpha_values_meta_ITS_noNA %>%
                                                            filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["Age2_shannon"]] <- kruskal.test(shannon_entropy ~ Age2, 
                                                      data = alpha_values_meta_ITS_noNA %>%
                                                        filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["HospStatus_evenness"]] <- kruskal.test(pielou_evenness ~ HospStatus, 
                                                             data = alpha_values_meta_ITS_noNA %>%
                                                               filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["TurtleSex_evenness"]] <- kruskal.test(pielou_evenness ~ TurtleSex, 
                                                            data = alpha_values_meta_ITS_noNA %>% 
                                                              filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["DetSex_evenness"]] <- kruskal.test(pielou_evenness ~ DetSex, 
                                                         data = alpha_values_meta_ITS_noNA %>%
                                                           filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["AgeRange_evenness"]] <- kruskal.test(pielou_evenness ~ AgeRange, 
                                                           data = alpha_values_meta_ITS_noNA %>%
                                                             filter(SampleSite == "CLOACA"))
kruskall_cloaca_ITS[["Age2_evenness"]] <- kruskal.test(pielou_evenness ~ Age2, 
                                                       data = alpha_values_meta_ITS_noNA %>%
                                                         filter(SampleSite == "CLOACA"))
capture.output(kruskall_cloaca_ITS, file = "r_output/ITS2_endo/alpha_diversity-Kruskall-cloaca.txt")

kruskall_wtr_ITS <- list()
kruskall_wtr_ITS[["HospStatus_faith"]] <- kruskal.test(faith_pd ~ HospStatus, 
                                                          data = alpha_values_meta_ITS_noNA %>%
                                                            filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["TurtleSex_faith"]] <- kruskal.test(faith_pd ~ TurtleSex, 
                                                         data = alpha_values_meta_ITS_noNA %>% 
                                                           filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["DetSex_faith"]] <- kruskal.test(faith_pd ~ DetSex, 
                                                      data = alpha_values_meta_ITS_noNA %>%
                                                        filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["AgeRange_faith"]] <- kruskal.test(faith_pd ~ AgeRange, 
                                                        data = alpha_values_meta_ITS_noNA %>%
                                                          filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["Age2_faith"]] <- kruskal.test(faith_pd ~ Age2, 
                                                    data = alpha_values_meta_ITS_noNA %>%
                                                      filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["HospStatus_observed"]] <- kruskal.test(observed_features ~ HospStatus, 
                                                             data = alpha_values_meta_ITS_noNA %>%
                                                               filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["TurtleSex_observed"]] <- kruskal.test(observed_features ~ TurtleSex, 
                                                            data = alpha_values_meta_ITS_noNA %>% 
                                                              filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["DetSex_observed"]] <- kruskal.test(observed_features ~ DetSex, 
                                                         data = alpha_values_meta_ITS_noNA %>%
                                                           filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["AgeRange_observed"]] <- kruskal.test(observed_features ~ AgeRange, 
                                                           data = alpha_values_meta_ITS_noNA %>%
                                                             filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["Age2_observed"]] <- kruskal.test(observed_features ~ Age2, 
                                                       data = alpha_values_meta_ITS_noNA %>%
                                                         filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["HospStatus_shannon"]] <- kruskal.test(shannon_entropy ~ HospStatus, 
                                                          data = alpha_values_meta_ITS_noNA %>%
                                                            filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["TurtleSex_shannon"]] <- kruskal.test(shannon_entropy ~ TurtleSex, 
                                                           data = alpha_values_meta_ITS_noNA %>% 
                                                             filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["DetSex_shannon"]] <- kruskal.test(shannon_entropy ~ DetSex, 
                                                        data = alpha_values_meta_ITS_noNA %>%
                                                          filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["AgeRange_shannon"]] <- kruskal.test(shannon_entropy ~ AgeRange, 
                                                          data = alpha_values_meta_ITS_noNA %>%
                                                            filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["Age2_shannon"]] <- kruskal.test(shannon_entropy ~ Age2, 
                                                      data = alpha_values_meta_ITS_noNA %>%
                                                        filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["HospStatus_evenness"]] <- kruskal.test(pielou_evenness ~ HospStatus, 
                                                             data = alpha_values_meta_ITS_noNA %>%
                                                               filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["TurtleSex_evenness"]] <- kruskal.test(pielou_evenness ~ TurtleSex, 
                                                            data = alpha_values_meta_ITS_noNA %>% 
                                                              filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["DetSex_evenness"]] <- kruskal.test(pielou_evenness ~ DetSex, 
                                                         data = alpha_values_meta_ITS_noNA %>%
                                                           filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["AgeRange_evenness"]] <- kruskal.test(pielou_evenness ~ AgeRange, 
                                                           data = alpha_values_meta_ITS_noNA %>%
                                                             filter(SampleSite == "TANK WATER"))
kruskall_wtr_ITS[["Age2_evenness"]] <- kruskal.test(pielou_evenness ~ Age2, 
                                                       data = alpha_values_meta_ITS_noNA %>%
                                                         filter(SampleSite == "TANK WATER"))
capture.output(kruskall_wtr_ITS, file = "r_output/ITS2_endo/alpha_diversity-Kruskall-tank water.txt")

#####PCOAS CLEAN UP
##----- Import data for PCoA plots -----

##----- Import data for PCoA plots -----

# load taxonomies
taxonomy_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/taxonomy.qza")

##----- Load data from qiime2 qza to R for PCoA/PCA plots
# load qiime pca and pcoa results.qza
# ITS data loading
bray_pcoa_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/bray_curtis_pcoa_results.qza")
unifrac_pcoa_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/unweighted_unifrac_pcoa_results.qza")
w_unifrac_pcoa_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/weighted_unifrac_pcoa_results.qza")
jaccard_pcoa_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/jaccard_pcoa_results.qza")
rAitch_pca_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/no-neg-ctrl-ordination.qza")

##----- PCoA plots

# define theme for pcoa/pca plots
theme_pcoa1 <- theme(
  panel.border = element_rect(fill = "transparent", color = "grey30"),
  plot.subtitle = element_text(size = 11),
  axis.title.x = element_text(size = 11, color = "black"),
  axis.title.y = element_text(size = 11, color = "black"),
  axis.text.x = element_text(size = 9, color = "grey30"),
  axis.text.y = element_text(size = 9, color = "grey30"),
  legend.title = element_text(size = 11, color = "black"),
  legend.text = element_markdown(size = 10, color = "grey30"),
  legend.spacing.y = unit(0.05, 'cm'),
  legend.key.size = unit(0.5, 'cm'),
  #legend.position = "bottom",
  #legend.box = "horizontal",
  panel.grid.minor.x = element_line(linewidth = 0, color = "transparent"),
  panel.grid.minor.y = element_line(linewidth = 0, color = "transparent"),
  panel.grid.major.x = element_line(linewidth = 0, color = "transparent"),
  panel.grid.major.y = element_line(linewidth = 0, color = "transparent"),
  axis.line.x.bottom = element_line(linewidth = 0.3, color = "grey30"),
  axis.line.y.left = element_line(linewidth = 0.3, color = "grey30"),
  axis.ticks = element_line(colour = "grey30")
)

bray1_plot_ITS <-
  bray_pcoa_ITS$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ITS) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  geom_line(aes(group = SamplingEvent), lty = 2, colour = "gray40") +
  geom_point(
  aes(fill = SampleSite, shape = SampleSite),
  alpha=0.8,
  size = 3
  ) +
  geom_point(
    aes(shape = SampleSite),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  # geom_text_repel(aes(label = SampleID), size = 1.75) +
  scale_shape_manual(values = c(23, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                NEGCTRL = "Negative control",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               "TANK WATER" = "Tank water")) +
  scale_color_manual(name = "Sample site:", 
                     values = c("#93624d", "#89d6dc"),
                     labels = c(CLOACA = "Cloacal",
                                "TANK WATER" = "Tank water")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*bray_pcoa_ITS$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*bray_pcoa_ITS$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Bray-Curtis",
       tag = "a)")

jacc1_plot_ITS <-
  jaccard_pcoa_ITS$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ITS) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  geom_line(aes(group = SamplingEvent), lty = 2, colour = "gray40") +
  geom_point(
    aes(fill = SampleSite, shape = SampleSite),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes(shape = SampleSite),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  # geom_text_repel(aes(label = SampleID), size = 1.75) +
  scale_shape_manual(values = c(23, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                NEGCTRL = "Negative control",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               "TANK WATER" = "Tank water")) +
  scale_color_manual(name = "Sample site:", 
                     values = c("#93624d", "#89d6dc"),
                     labels = c(CLOACA = "Cloacal",
                                "TANK WATER" = "Tank water")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*jaccard_pcoa_ITS$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*jaccard_pcoa_ITS$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Jaccard",
       tag = "b)")

unwunifrac1_plot_ITS <-
  unifrac_pcoa_ITS$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ITS) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  geom_line(aes(group = SamplingEvent), lty = 2, colour = "gray40") +
  geom_point(
    aes(fill = SampleSite, shape = SampleSite),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes(shape = SampleSite),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  # geom_text_repel(aes(label = SampleID), size = 1.75) +
  scale_shape_manual(values = c(23, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                NEGCTRL = "Negative control",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               "TANK WATER" = "Tank water")) +
  scale_color_manual(name = "Sample site:", 
                     values = c("#93624d", "#89d6dc"),
                     labels = c(CLOACA = "Cloacal",
                                "TANK WATER" = "Tank water")) +  guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*unifrac_pcoa_ITS$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*unifrac_pcoa_ITS$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Unweighted UniFrac",
       tag = "c)")

wunifrac1_plot_ITS <-
  w_unifrac_pcoa_ITS$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_ITS) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  geom_line(aes(group = SamplingEvent), lty = 2, colour = "gray40") +
  geom_point(
    aes(fill = SampleSite, shape = SampleSite),
    alpha=0.8,
    size = 3
  ) +
  geom_point(
    aes(shape = SampleSite),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  # geom_text_repel(aes(label = SampleID), size = 1.5) +
  scale_shape_manual(values = c(23, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                NEGCTRL = "Negative control",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               "TANK WATER" = "Tank water")) +
  scale_color_manual(name = "Sample site:", 
                     values = c("#93624d", "#89d6dc"),
                     labels = c(CLOACA = "Cloacal",
                                "TANK WATER" = "Tank water")) +  
    guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*w_unifrac_pcoa_ITS$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*w_unifrac_pcoa_ITS$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Weighted UniFrac",
       tag = "d)")

##---- PCA plots
# rename column in taxonomy so we can join it to our data
taxonomy_ITS <- rename(taxonomy_ITS$data, FeatureID = Feature.ID)

rPCA1_plot_ITS <- 
  ggplot() +
  geom_segment(data = rAitch_pca_ITS$data$Species %>% #<<
                 mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                 slice_max(n = 5, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                 mutate(PC1=PC1*3, PC2=PC2*3) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                 left_join(taxonomy_ITS), #<<
               aes(x=0, xend=PC1, y=0, yend=PC2),
               linewidth = 0.3,
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30") +
  geom_line(
    data = rAitch_pca_ITS$data$Vectors %>% #<<
      left_join(metadata_ITS), #<<
    aes(x = PC1, 
        y = PC2, 
        group = SamplingEvent), lty = 2, colour = "gray40") +
  geom_point(
    data = rAitch_pca_ITS$data$Vectors %>% #<<
      left_join(metadata_ITS), #<<
    aes(
      x = PC1, 
      y = PC2, 
      fill = SampleSite, 
      shape = SampleSite
    ),
    alpha = 0.8,
    size = 3
  ) +
  geom_point(
    data = rAitch_pca_ITS$data$Vectors %>% #<<
      left_join(metadata_ITS), #<<
    aes(
      x = PC1, 
      y = PC2, 
      shape = SampleSite
    ),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  # geom_text_repel(data = rAitch_pca_ITS$data$Vectors %>%
  #                   left_join(metadata_ITS),
  #                 aes(x = PC1, y = PC2, label = TurtleName), size = 1.75) +
    scale_shape_manual(values = c(23, 25), 
                       name = "Sample site:",
                       labels = c(CLOACA = "Cloacal",
                                  "TANK WATER" = "Tank water")) +
    scale_fill_manual(name = "Sample site:", 
                      values = c("#93624d", "#89d6dc"),
                      labels = c(CLOACA = "Cloacal",
                                 "TANK WATER" = "Tank water")) +
    scale_color_manual(name = "Sample site:", 
                       values = c("#93624d", "#89d6dc"),
                       labels = c(CLOACA = "Cloacal",
                                  "TANK WATER" = "Tank water")) + 
    geom_label_repel(data = rAitch_pca_ITS$data$Species %>% #<<
                     mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                     slice_max(n = 5, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                     mutate(PC1=PC1*3, PC2=PC2*3) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                     left_join(taxonomy_ITS) %>% #<<
                     mutate(Taxa = c("plain('family')~italic('Nectriaceae')",
                                     "plain('order Xylariales')~italic('')",
                                     "plain(' ')~italic('Preussia flanaganii')",
                                     "plain('genus')~italic('Rhizoctonia')",
                                     "plain('genus')~italic('Tetracladium')")),
                   aes(x = PC1, y = PC2, 
                       label = Taxa),
                   point.padding = 2,
                   alpha = 0.8,
                   colour = "black",
                   size = 3,
                   parse = TRUE) +
  xlab(paste("PC1 (",round(100*rAitch_pca_ITS$data$ProportionExplained[1],2),"%)")) + #<<
  ylab(paste("PC2 (",round(100*rAitch_pca_ITS$data$ProportionExplained[2],2),"%)")) + #<<
  theme_pcoa1 +
  labs(subtitle = NULL,
       tag = "e)")

PCoA_plots_collected_ITS <- ((bray1_plot_ITS + jacc1_plot_ITS) / (unwunifrac1_plot_ITS + wunifrac1_plot_ITS)) + plot_layout(guides = 'collect')

ggsave(
  filename = "r_output/ITS2_endo/pcoa_plots.pdf",
  plot = PCoA_plots_collected_ITS,
  device = cairo_pdf,
  height = 140,
  width = 160,
  units = "mm"
)

ggsave(
  filename = "r_output/ITS2_endo/rpca_plot.pdf",
  plot = rPCA1_plot_ITS,
  device = cairo_pdf,
  height = 140,
  width = 160,
  units = "mm"
)

unwunifrac_ITS_pub <- unwunifrac1_plot_ITS + theme_pcoa1_pub + labs(subtitle = NULL, tag = "a)")
rPCA_ITS_pub <- rPCA1_plot_ITS + theme_pcoa1_pub + labs(subtitle = NULL, tag = "b)")

supplement_fig1 <- unwunifrac_ITS_pub/rPCA_ITS_pub+plot_layout(guides = "collect")

ggsave(
  filename = "r_output/ITS2_endo/supp_beta.pdf",
  plot = supplement_fig1,
  device = cairo_pdf,
  height = 200,
  width = 140,
  units = "mm"
)

##----- ADONIS PERMANOVA -----
bray_dist_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/bray_curtis_distance_matrix.qza")
jacc_dist_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/jaccard_distance_matrix.qza")
uunif_dist_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/unweighted_unifrac_distance_matrix.qza")
wunif_dist_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/weighted_unifrac_distance_matrix.qza")
raitch_dist_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/no-neg-ctrl-distance-matrix.qza")

bray_dist_data_ITS <- bray_dist_ITS$data
jacc_dist_data_ITS <- jacc_dist_ITS$data
uunif_dist_data_ITS <- uunif_dist_ITS$data
wunif_dist_data_ITS <- wunif_dist_ITS$data
raitch_dist_data_ITS <- raitch_dist_ITS$data

##
removed_samples_ITS <- c("ITSNEGCTRL",
                        "ITS0118C",
                        "ITS0117W",
                        "ITS0146C")
removed_samples_ITS_aitch <- c("ITSNEGCTRL",
                               "ITS0146C")

adonis_metadata_ITS <- filter(metadata_ITS, !SampleID %in% removed_samples_ITS)
adonis_metadata_ITS_aitch <- filter(metadata_ITS, !SampleID %in% removed_samples_ITS_aitch)

# distances and metadata need to be sorted in the same way (in this cae alphabetically)
adonis_metadata_ITS <- adonis_metadata_ITS[order(adonis_metadata_ITS$SampleID),]
adonis_metadata_ITS_aitch <- adonis_metadata_ITS_aitch[order(adonis_metadata_ITS_aitch$SampleID),]

set.seed(123)
adonis_results_ITS <- list()
adonis_results_ITS[["Bray-Curtis ADONIS"]] <- adonis2(formula = bray_dist_data_ITS~SampleSite, data = adonis_metadata_ITS, permutations = 999, by = "terms")
adonis_results_ITS[["Jaccard ADONIS"]] <- adonis2(formula = jacc_dist_data_ITS~SampleSite, data = adonis_metadata_ITS, permutations = 999, by = "terms")
adonis_results_ITS[["unw. UniFrac ADONIS"]] <- adonis2(formula = uunif_dist_data_ITS~SampleSite, data = adonis_metadata_ITS, permutations = 999, by = "terms")
adonis_results_ITS[["w. UniFrac ADONIS"]] <- adonis2(formula = wunif_dist_data_ITS~SampleSite, data = adonis_metadata_ITS, permutations = 999, by = "terms")
adonis_results_ITS[["r. Aitchison ADONIS"]] <- adonis2(formula = raitch_dist_data_ITS~SampleSite, data = adonis_metadata_ITS_aitch, permutations = 999, by = "terms")
capture.output(adonis_results_ITS, file = "r_output/ITS2_endo/adonis_permanova-all_sample_sites.txt")

