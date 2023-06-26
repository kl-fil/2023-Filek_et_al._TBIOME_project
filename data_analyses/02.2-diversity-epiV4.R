## Dataviz script for TurtleBiome Holobiont 2023 manuscript
## Processing and visualizing data related to sequencing results for
## endozoicepizoic and tank water samples (16S rRNA gene, trimmed to V4 region)
## Author: Klara Filek

##--------------------load packages--------------------
library(tidyverse) #v2.0.0
library(ggplot2) #v3.4.2
library(qiime2R) #qiime to R import qza v0.99.6
library(patchwork) #collect multiple ggplots together v1.1.2
library(ggrepel) #repel text/labels v0.9.3
library(ggforce) #v0.4.1
library(ggtext) #v0.1.2
library(vegan) #v2.6-4
library(pairwiseAdonis) #v0.4.1
library(ggpubr) #v0.6.0s

#set global theme for ggplot
theme_set(theme_light(base_size = 11, base_family = "Arial")) #<< font!

#vectors for colors per sampling site
c("#93624d", "#ff9999", "#3BB848", "#89d6dc")
c("#462f25","#750000", "#143e18", "#16474b")

# "#<<" comment indicates the line of code that needs to be changed when using 
# a different dataset for the same code

##---- make directory for visualizations ----
dir.create("r_output/16S_epi_v4_merged/", recursive = TRUE)

##---- Load metadata and data-----
# load metadata in tsv format
metadata_epiv4 <- readr::read_tsv("master_metadata/16S_all_metadata.tsv")

##---- Alpha diversity ----

shannonv4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/shannon_vector.qza")
f_pdv4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/faith_pd_vector.qza")
evennessv4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/evenness_vector.qza")
observed_asvv4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/observed_features_vector.qza")

alpha_values_metav4 <- 
  metadata_epiv4 %>%
  left_join(rownames_to_column(shannonv4$data, "SampleID")) %>%
  left_join(rename(f_pdv4$data, "SampleID" = "V1", "faith_PD" = "V2")) %>%
  left_join(rownames_to_column(evennessv4$data, "SampleID")) %>%
  left_join(rownames_to_column(observed_asvv4$data, "SampleID")) %>%
  drop_na(observed_features)

alpha_values_metav4$SampleSite <- ordered(alpha_values_metav4$SampleSite,
                                          levels = c("CLOACA", "ORAL", "CARAPACE", "TANK WATER"))

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

##----- Alpha div plots -----

observed_alpha_rainv4 <-
  alpha_values_metav4 %>%
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
               color = c("#462f25","#750000", "#143e18", "#16474b")
  ) +
  geom_point(aes(fill = SampleSite),
             shape = 21,
             color = "black",
             size = 1.5,
             alpha = .6,
             position = position_jitter(
               seed = 1, width = .1
             )
  ) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
                              CARAPACE = "Carapace \n(n=15)",
                              "TANK WATER" = "Tank \nwater\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Observed features (ASVs)", #<<
    tag = "a)" #<<
  ) + 
  theme_alpha +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

faith_alpha_rainv4 <-
  alpha_values_metav4 %>%
  ggplot(aes(x = SampleSite, y = faith_PD)) +
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
               color = c("#462f25","#750000", "#143e18", "#16474b")
  ) +
  geom_point(aes(fill = SampleSite),
             shape = 21,
             color = "black",
             size = 1.5,
             alpha = .6,
             position = position_jitter(
               seed = 1, width = .1
             )
  ) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
                              CARAPACE = "Carapace \n(n=15)",
                              "TANK WATER" = "Tank \nwater\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Faith's phylogenetic diversity", #<<
    tag = "b)" #<<
  ) + 
  theme_alpha +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

shannon_alpha_rainv4 <-
  alpha_values_metav4 %>%
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
               color = c("#462f25","#750000", "#143e18", "#16474b")
  ) +
  geom_point(aes(fill = SampleSite),
             shape = 21,
             color = "black",
             size = 1.5,
             alpha = .6,
             position = position_jitter(
               seed = 1, width = .1
             )
  ) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
                              CARAPACE = "Carapace \n(n=15)",
                              "TANK WATER" = "Tank \nwater\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Shannon's entropy", #<<
    tag = "c)" #<<
  ) + 
  theme_alpha +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

evenness_alpha_rainv4 <-
  alpha_values_metav4 %>%
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
               color = c("#462f25","#750000", "#143e18", "#16474b")
  ) +
  geom_point(aes(fill = SampleSite),
             shape = 21,
             color = "black",
             size = 1.5,
             alpha = .6,
             position = position_jitter(
               seed = 1, width = .1
             )
  ) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
                              CARAPACE = "Carapace \n(n=15)",
                              "TANK WATER" = "Tank \nwater\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Pielou's evenness", #<<
    tag = "d)" #<<
  ) + 
  theme_alpha +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

rain_plotsv4 <- (observed_alpha_rainv4+faith_alpha_rainv4)/(shannon_alpha_rainv4+evenness_alpha_rainv4)

ggsave(
  "r_output/16S_epi_v4_merged/alpha-div-collected.pdf",
  plot = rain_plotsv4,
  device = cairo_pdf,
  height = 150,
  width = 168,
  units = "mm"
)

#----- Alpha diversity statistics -----

alpha_values_metav4$SampleSite <- ordered(alpha_values_metav4$SampleSite,
                                        levels = c("CLOACA", "ORAL", "CARAPACE", "TANK WATER"))

observed_summarystats_v4 <- group_by(alpha_values_metav4, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(observed_features, na.rm = TRUE),
    sd = sd(observed_features, na.rm = TRUE),
    median = median(observed_features, na.rm = TRUE),
    IQR = IQR(observed_features, na.rm = TRUE)
  )

faith_summarystats_v4 <- group_by(alpha_values_metav4, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(faith_PD, na.rm = TRUE),
    sd = sd(faith_PD, na.rm = TRUE),
    median = median(faith_PD, na.rm = TRUE),
    IQR = IQR(faith_PD, na.rm = TRUE)
  )

shannon_summarystats_v4 <- group_by(alpha_values_metav4, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(shannon_entropy, na.rm = TRUE),
    sd = sd(shannon_entropy, na.rm = TRUE),
    median = median(shannon_entropy, na.rm = TRUE),
    IQR = IQR(shannon_entropy, na.rm = TRUE)
  )

evenness_summarystats_v4 <- group_by(alpha_values_metav4, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(pielou_evenness, na.rm = TRUE),
    sd = sd(pielou_evenness, na.rm = TRUE),
    median = median(pielou_evenness, na.rm = TRUE),
    IQR = IQR(pielou_evenness, na.rm = TRUE)
  )

summary_stats_alpha_v4 <- list()
summary_stats_alpha_v4[["Observed features"]] <- observed_summarystats_v4
summary_stats_alpha_v4[["Faith PD"]] <- faith_summarystats_v4
summary_stats_alpha_v4[["Shannon entropy"]] <- shannon_summarystats_v4
summary_stats_alpha_v4[["Pielou evenness"]] <- evenness_summarystats_v4

write_tsv(do.call(cbind, summary_stats_alpha_v4), "r_output/16S_epi_v4_merged/alpha_diversity_summary_stats.tsv")
capture.output(summary_stats_alpha_v4, file = "r_output/16S_epi_v4_merged/alpha_diversity_summary_stats.txt")

##---- Kruskall-Wallis test & pairwise comparisons (Wilcox)
krustats_obs_v4 <- list()
krustats_obs_v4[["Kruskall-Wallis test"]] <- kruskal.test(observed_features ~ SampleSite, data = alpha_values_metav4)
krustats_obs_v4[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_metav4$observed_features, alpha_values_metav4$SampleSite,
                                                                            p.adjust.method = "BH")
capture.output(krustats_obs_v4, file = "r_output/16S_epi_v4_merged/alpha_diversity-Kruskall-observed.txt")

krustats_faith_v4 <- list()
krustats_faith_v4[["Kruskall-Wallis test"]] <- kruskal.test(faith_PD ~ SampleSite, data = alpha_values_metav4)
krustats_faith_v4[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_metav4$faith_PD, alpha_values_metav4$SampleSite,
                                                                              p.adjust.method = "BH")
capture.output(krustats_faith_v4, file = "r_output/16S_epi_v4_merged/alpha_diversity-Kruskall-faith_pd.txt")

krustats_shan_v4 <- list()
krustats_shan_v4[["Kruskall-Wallis test"]] <- kruskal.test(shannon_entropy ~ SampleSite, data = alpha_values_metav4)
krustats_shan_v4[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_metav4$shannon_entropy, alpha_values_metav4$SampleSite,
                                                                             p.adjust.method = "BH")
capture.output(krustats_shan_v4, file = "r_output/16S_epi_v4_merged/alpha_diversity-Kruskall-shannon.txt")

krustats_even_v4 <- list()
krustats_even_v4[["Kruskall-Wallis test"]] <- kruskal.test(pielou_evenness ~ SampleSite, data = alpha_values_metav4)
krustats_even_v4[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_metav4$pielou_evenness, alpha_values_metav4$SampleSite,
                                                                             p.adjust.method = "BH")
capture.output(krustats_even_v4, file = "r_output/16S_epi_v4_merged/alpha_diversity-Kruskall-evenness.txt")

#----- Correlation between alpha diversity and CCL per sample site -----

pears_faith_clov4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("CLOACA")) %>%
  ggplot(aes(x = CCL, y = faith_PD)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "CLOACA"), method=lm , fill = "gray80", color="#462f25", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "CLOACA"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#93624d", shape = 21, size = 2) + 
  labs(title = "Cloaca",
       x = "CCL (cm)",
       y = "Faith's PD",
       tag = "a)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank())

pears_faith_orlv4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("ORAL")) %>%
  ggplot(aes(x = CCL, y = faith_PD)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "ORAL"), method=lm , fill = "gray80", color="#750000", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "ORAL"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#ff9999", shape = 21, size = 2) + 
  labs(title = "Oral",
       x = "CCL (cm)",
       y = "Faith's PD",
       tag = "b)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pears_faith_carv4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("CARAPACE")) %>%
  ggplot(aes(x = CCL, y = faith_PD)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "CARAPACE"), method=lm , fill = "gray80", color="#143e18", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "CARAPACE"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#3BB848", shape = 21, size = 2) + 
  labs(title = "Carapace",
       x = "CCL (cm)",
       y = "Faith's PD",
       tag = "c)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pears_faith_wtrv4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("TANK WATER")) %>%
  ggplot(aes(x = CCL, y = faith_PD)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "TANK WATER"), method=lm , fill = "gray80", color="#16474b", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "TANK WATER"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#89d6dc", shape = 21, size = 2) + 
  labs(title = "Tank Water",
       x = "CCL (cm)",
       y = "Faith's PD",
       tag = "d)") +
  theme_alpha+ 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pears_obs_clov4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("CLOACA")) %>%
  ggplot(aes(x = CCL, y = observed_features)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "CLOACA"), method=lm , fill = "gray80", color="#462f25", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "CLOACA"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#93624d", shape = 21, size = 2) + 
  labs(title = "Cloaca",
       x = "CCL (cm)",
       y = "ASVs",
       tag = "e)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank())

pears_obs_orlv4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("ORAL")) %>%
  ggplot(aes(x = CCL, y = observed_features)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "ORAL"), method=lm , fill = "gray80", color="#750000", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "ORAL"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#ff9999", shape = 21, size = 2) + 
  labs(title = "Oral",
       x = "CCL (cm)",
       y = "ASVs",
       tag = "f)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pears_obs_carv4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("CARAPACE")) %>%
  ggplot(aes(x = CCL, y = observed_features)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "CARAPACE"), method=lm , fill = "gray80", color="#143e18", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "CARAPACE"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#3BB848", shape = 21, size = 2) + 
  labs(title = "Carapace",
       x = "CCL (cm)",
       y = "ASVs",
       tag = "g)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pears_obs_wtrv4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("TANK WATER")) %>%
  ggplot(aes(x = CCL, y = observed_features)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "TANK WATER"), method=lm , fill = "gray80", color="#16474b", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "TANK WATER"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#89d6dc", shape = 21, size = 2) + 
  labs(title = "Tank Water",
       x = "CCL (cm)",
       y = "ASVs",
       tag = "h)") +
  theme_alpha+ 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pears_sha_clov4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("CLOACA")) %>%
  ggplot(aes(x = CCL, y = shannon_entropy)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "CLOACA"), method=lm , fill = "gray80", color="#462f25", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "CLOACA"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#93624d", shape = 21, size = 2) + 
  labs(title = "Cloaca",
       x = "CCL (cm)",
       y = "Shannon's",
       tag = "i)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line())

pears_sha_orlv4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("ORAL")) %>%
  ggplot(aes(x = CCL, y = shannon_entropy)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "ORAL"), method=lm , fill = "gray80", color="#750000", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "ORAL"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#ff9999", shape = 21, size = 2) + 
  labs(title = "Oral",
       x = "CCL (cm)",
       y = "Shannon's",
       tag = "j)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.y = element_blank())

pears_sha_carv4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("CARAPACE")) %>%
  ggplot(aes(x = CCL, y = shannon_entropy)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "CARAPACE"), method=lm , fill = "gray80", color="#143e18", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "CARAPACE"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#3BB848", shape = 21, size = 2) + 
  labs(title = "Carapace",
       x = "CCL (cm)",
       y = "Shannon's",
       tag = "k)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.y = element_blank())

pears_sha_wtrv4 <- alpha_values_metav4 %>%
  filter(SampleSite %in% c("TANK WATER")) %>%
  ggplot(aes(x = CCL, y = shannon_entropy)) +
  geom_smooth(data = alpha_values_metav4 %>%
                filter(SampleSite== "TANK WATER"), method=lm , fill = "gray80", color="#16474b", se=TRUE) +
  stat_cor(data = alpha_values_metav4 %>%
             filter(SampleSite== "TANK WATER"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#89d6dc", shape = 21, size = 2) + 
  labs(title = "Tank Water",
       x = "CCL (cm)",
       y = "Shannon's",
       tag = "l)") +
  theme_alpha+ 
  theme(panel.grid.major.x = element_line(),
        axis.title.y = element_blank())

pears_alphav4 <- pears_faith_clov4+pears_faith_orlv4+pears_faith_carv4+pears_faith_wtrv4+
  pears_obs_clov4+pears_obs_orlv4+pears_obs_carv4+pears_obs_wtrv4+ 
  pears_sha_clov4+pears_sha_orlv4+pears_sha_carv4+pears_sha_wtrv4+plot_layout(ncol = 4)

ggsave(
  "r_output/16S_epi_v4_merged/pears-alpha-collected.pdf",
  plot = pears_alphav4,
  device = cairo_pdf,
  height = 210,
  width = 250,
  units = "mm"
)

pearson_corr_cloaca_v4 <- list()
pearson_corr_cloaca_v4[["Observed features"]] <- with(alpha_values_metav4 %>%
                                                     filter(SampleSite == "CLOACA"),cor.test(CCL,observed_features, method = "pearson"))
pearson_corr_cloaca_v4[["Faith's PD"]] <- with(alpha_values_metav4 %>%
                                              filter(SampleSite == "CLOACA"),cor.test(CCL,faith_PD, method = "pearson"))
pearson_corr_cloaca_v4[["Shannon's entropy"]] <- with(alpha_values_metav4 %>%
                                                     filter(SampleSite == "CLOACA"),cor.test(CCL,shannon_entropy, method = "pearson"))
pearson_corr_cloaca_v4[["Pielou's evenness"]] <- with(alpha_values_metav4 %>%
                                                     filter(SampleSite == "CLOACA"),cor.test(CCL,pielou_evenness, method = "pearson"))
capture.output(pearson_corr_cloaca_v4, file = "r_output/16S_epi_v4_merged/alpha_diversity-pearson-cloaca.txt")

pearson_corr_oral_v4 <- list()
pearson_corr_oral_v4[["Observed features"]] <- with(alpha_values_metav4 %>%
                                                   filter(SampleSite == "ORAL"),cor.test(CCL,observed_features, method = "pearson"))
pearson_corr_oral_v4[["Faith's PD"]] <- with(alpha_values_metav4 %>%
                                            filter(SampleSite == "ORAL"),cor.test(CCL,faith_PD, method = "pearson"))
pearson_corr_oral_v4[["Shannon's entropy"]] <- with(alpha_values_metav4 %>%
                                                   filter(SampleSite == "ORAL"),cor.test(CCL,shannon_entropy, method = "pearson"))
pearson_corr_oral_v4[["Pielou's evenness"]] <- with(alpha_values_metav4 %>%
                                                   filter(SampleSite == "ORAL"),cor.test(CCL,pielou_evenness, method = "pearson"))
capture.output(pearson_corr_oral_v4, file = "r_output/16S_epi_v4_merged/alpha_diversity-pearson-oral.txt")

pearson_corr_tw_v4 <- list()
pearson_corr_tw_v4[["Observed features"]] <- with(alpha_values_metav4 %>%
                                                 filter(SampleSite == "TANK WATER"),cor.test(CCL,observed_features, method = "pearson"))
pearson_corr_tw_v4[["Faith's PD"]] <- with(alpha_values_metav4 %>%
                                          filter(SampleSite == "TANK WATER"),cor.test(CCL,faith_PD, method = "pearson"))
pearson_corr_tw_v4[["Shannon's entropy"]] <- with(alpha_values_metav4 %>%
                                                 filter(SampleSite == "TANK WATER"),cor.test(CCL,shannon_entropy, method = "pearson"))
pearson_corr_tw_v4[["Pielou's evenness"]] <- with(alpha_values_metav4 %>%
                                                 filter(SampleSite == "TANK WATER"),cor.test(CCL,pielou_evenness, method = "pearson"))
capture.output(pearson_corr_tw_v4, file = "r_output/16S_epi_v4_merged/alpha_diversity-pearson-tank_water.txt")

pearson_corr_car_v4 <- list()
pearson_corr_car_v4[["Observed features"]] <- with(alpha_values_metav4 %>%
                                                    filter(SampleSite == "CARAPACE"),cor.test(CCL,observed_features, method = "pearson"))
pearson_corr_car_v4[["Faith's PD"]] <- with(alpha_values_metav4 %>%
                                             filter(SampleSite == "CARAPACE"),cor.test(CCL,faith_PD, method = "pearson"))
pearson_corr_car_v4[["Shannon's entropy"]] <- with(alpha_values_metav4 %>%
                                                    filter(SampleSite == "CARAPACE"),cor.test(CCL,shannon_entropy, method = "pearson"))
pearson_corr_car_v4[["Pielou's evenness"]] <- with(alpha_values_metav4 %>%
                                                    filter(SampleSite == "CARAPACE"),cor.test(CCL,pielou_evenness, method = "pearson"))
capture.output(pearson_corr_car_v4, file = "r_output/16S_epi_v4_merged/alpha_diversity-pearson-carapace.txt")

##----- Import data for PCoA plots -----

# load taxonomies
taxonomy_epiv4 <- read_qza("qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza")
metadata_epiv4$SampleSite <- ordered(metadata_epiv4$SampleSite,
                                     levels = c("CLOACA", "ORAL", "CARAPACE", "TANK WATER"))
##---- Load data from qiime2 qza to R for PCoA/PCA plots
# load qiime pca and pcoa results.qza
# 16S data loading (from core-metrics-results-merged_16S-0-with-phyla-no-mitochondria-no-chloroplast-filtered-phylogeny )

bray_pcoa_v4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/bray_curtis_pcoa_results.qza")
unifrac_pcoa_v4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/unweighted_unifrac_pcoa_results.qza")
w_unifrac_pcoa_v4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/weighted_unifrac_pcoa_results.qza")
jaccard_pcoa_v4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/jaccard_pcoa_results.qza")
rAitch_pca_v4 <- read_qza("qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination-29000.qza")

##---- PCoA plots

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

bray1v4_plot <-
  bray_pcoa_v4$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_epiv4) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
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
  scale_shape_manual(values = c(23, 21, 22, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                CARAPACE = "Carapace",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*bray_pcoa_v4$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*bray_pcoa_v4$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Bray-Curtis",
       tag = "a)")

jacc1v4_plot <-
  jaccard_pcoa_v4$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_epiv4) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
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
  scale_shape_manual(values = c(23, 21, 22, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                CARAPACE = "Carapace",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*jaccard_pcoa_v4$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*jaccard_pcoa_v4$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Jaccard",
       tag = "b)")

unwunifrac1v4_plot <-
  unifrac_pcoa_v4$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_epiv4) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
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
  scale_shape_manual(values = c(23, 21, 22, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                CARAPACE = "Carapace",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  xlab(paste("PC1 (", round(100*unifrac_pcoa_v4$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*unifrac_pcoa_v4$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Unweighted UniFrac",
       tag = "c)")

wunifrac1v4_plot <-
  w_unifrac_pcoa_v4$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_epiv4) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
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
  scale_shape_manual(values = c(23, 21, 22, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                CARAPACE = "Carapace",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  #scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  xlab(paste("PC1 (", round(100*w_unifrac_pcoa_v4$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*w_unifrac_pcoa_v4$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Weighted UniFrac",
       tag = "d)")


##---- PCA plots
# rename column in taxonomy so we can join it to our data
taxonomy_epiv4 <- rename(taxonomy_epiv4$data, FeatureID = Feature.ID)

rpca_v41_plot <- 
ggplot() +
  geom_segment(data = rAitch_pca_v4$data$Species %>% #<<
                 mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                 slice_max(n = 8, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                 mutate(PC1=PC1*3.5, PC2=PC2*3.5) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                 left_join(taxonomy_epiv4), #<<
               aes(x=0, xend=PC1, y=0, yend=PC2),
               linewidth = 0.3,
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30") +
  geom_point(
    data = rAitch_pca_v4$data$Vectors %>% #<<
      left_join(metadata_epiv4), #<<
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
    data = rAitch_pca_v4$data$Vectors %>% #<<
      left_join(metadata_epiv4), #<<
    aes(
      x = PC1, 
      y = PC2, 
      shape = SampleSite
    ),
    size = 3,
    color = "black",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(23, 21, 22, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                CARAPACE = "Carapace",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")) +
  geom_label_repel(data = rAitch_pca_v4$data$Species %>% #<<
                     mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                     slice_max(n = 8, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                     mutate(PC1=PC1*3.5, PC2=PC2*3.5) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                     left_join(taxonomy_epiv4) %>% #<<
                     filter(FeatureID != "af46a1c0a6a2689dfb2631bb616076e1") %>%
                     mutate(Taxa = c("italic('Pseudoalteromonas')~plain('sp.')",
                                     "italic('Shewanella')~plain('sp.')",
                                     "italic('Vibrio')~plain('sp.')",
                                     "italic('Vibrio')~plain('sp.')",
                                     "italic('Vibrio')~plain('sp.')",
                                     #"plain('')~italic('')",
                                     "plain('unc.')~italic('Cardiobacteriaceae')",
                                     "plain('unc.')~italic('Marinifilum')")),
                   aes(x = PC1, y = PC2, 
                       label = Taxa),
                   point.padding = 1,
                   alpha = 0.8,
                   colour = "black",
                   size = 2,
                   parse = TRUE) +
  geom_label_repel(data = rAitch_pca_v4$data$Species %>% #<<
                     mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                     slice_max(n = 8, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                     mutate(PC1=PC1*3.5, PC2=PC2*3.5) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                     left_join(taxonomy_epiv4) %>% #<<
                     filter(FeatureID == "af46a1c0a6a2689dfb2631bb616076e1") %>%
                     mutate(Taxa = c("plain('unc.')~italic('Saccharospirillaceae')")),
                   aes(x = PC1, y = PC2, 
                       label = Taxa),
                   point.padding = 1,
                   nudge_x = -0.02,
                   nudge_y = -0.01,
                   alpha = 0.8,
                   colour = "black",
                   size = 2,
                   parse = TRUE) +
  xlab(paste("PC1 (",round(100*rAitch_pca_v4$data$ProportionExplained[1],2),"%)")) + #<<
  ylab(paste("PC2 (",round(100*rAitch_pca_v4$data$ProportionExplained[2],2),"%)")) + #<<
  theme_pcoa1 +
  labs(subtitle = "Robust Aitchison",
       tag = "e)")

PCoA_plots_collected_v4 <- ((bray1v4_plot + jacc1v4_plot) / (unwunifrac1v4_plot + wunifrac1v4_plot)) + plot_layout(guides = 'collect')

ggsave(filename = "r_output/16S_epi_v4_merged/pcoa_plots_v4.pdf",
       plot = PCoA_plots_collected_v4,
       device = cairo_pdf,
       height = 140,
       width = 160,
       units = "mm"
)

ggsave(filename = "r_output/16S_epi_v4_merged/rpca_plot_v4.pdf",
       plot = rpca_v41_plot,
       device = cairo_pdf,
       height = 130,
       width = 160,
       units = "mm"
)

##----- ADONIS PERMANOVA -----
bray_dist_v4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/bray_curtis_distance_matrix.qza")
jacc_dist_v4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/jaccard_distance_matrix.qza")
uunif_dist_v4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/unweighted_unifrac_distance_matrix.qza")
wunif_dist_v4 <- read_qza("qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/weighted_unifrac_distance_matrix.qza")
#raitch missing one sample in matrix
veg_table <- read_qza("qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza")
veg_table1 <- veg_table$data
tveg_table1 <- t(veg_table1)
rAitch_vegan_dist <- vegdist(x = tveg_table1,
                             method = "robust.aitchison")

bray_dist_data_v4 <- bray_dist_v4$data
jacc_dist_data_v4 <- jacc_dist_v4$data
uunif_dist_data_v4 <- uunif_dist_v4$data
wunif_dist_data_v4 <- wunif_dist_v4$data

##
removed_samples_v4 <- c("16SNEGCTRL",
                        "16S0094O",
                        "16S0113C",
                        "16S0118O",
                        "16S0118W",
                        "16S0119C",
                        "16S0092O",
                        "16S0064O")

adonis_metadata_v4 <- filter(metadata_epiv4, !SampleID %in% removed_samples_v4)

# distances and metadata need to be sorted in the same way (in this cae alphabetically)
adonis_metadata_v4 <- adonis_metadata_v4[order(adonis_metadata_v4$SampleID),]

set.seed(123)
adonis_results_v4 <- list()
adonis_results_v4[["Bray-Curtis ADONIS"]] <- adonis2(formula = bray_dist_data_v4~SampleSite, data = adonis_metadata_v4, permutations = 999, by = "terms")
adonis_results_v4[["Bray-Curtis Pairwise"]] <- pairwiseAdonis::pairwise.adonis2(bray_dist_data_v4~SampleSite, data = adonis_metadata_v4, permutations = 999)
adonis_results_v4[["Jaccard ADONIS"]] <- adonis2(formula = jacc_dist_data_v4~SampleSite, data = adonis_metadata_v4, permutations = 999, by = "terms")
adonis_results_v4[["Jaccard Pairwise"]] <- pairwiseAdonis::pairwise.adonis2(jacc_dist_data_v4~SampleSite, data = adonis_metadata_v4, permutations = 999)
adonis_results_v4[["unw. UniFrac ADONIS"]] <- adonis2(formula = uunif_dist_data_v4~SampleSite, data = adonis_metadata_v4, permutations = 999, by = "terms")
adonis_results_v4[["unw. UniFrac Pairwise"]] <- pairwiseAdonis::pairwise.adonis2(uunif_dist_data_v4~SampleSite, data = adonis_metadata_v4, permutations = 999)
adonis_results_v4[["w. UniFrac ADONIS"]] <- adonis2(formula = wunif_dist_data_v4~SampleSite, data = adonis_metadata_v4, permutations = 999, by = "terms")
adonis_results_v4[["w. UniFrac Pairwise"]] <- pairwiseAdonis::pairwise.adonis2(wunif_dist_data_v4~SampleSite, data = adonis_metadata_v4, permutations = 999)
adonis_results_v4[["r. Aitchison ADONIS"]] <- adonis2(formula = rAitch_vegan_dist~SampleSite, data = adonis_metadata_v4, permutations = 999, by = "terms")
adonis_results_v4[["r. Aitchison Pairwise"]] <- pairwiseAdonis::pairwise.adonis2(rAitch_vegan_dist~SampleSite, data = adonis_metadata_v4, permutations = 999)
capture.output(adonis_results_v4, file = "r_output/16S_epi_v4_merged/adonis_permanova-all_sample_sites.txt")

#----- Plots for publication
# Define themes
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
  plot.tag = element_text(size = 11)
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
  #legend.position = "bottom",
  #legend.box = "horizontal",
  panel.grid.minor.x = element_line(linewidth = 0, color = "transparent"),
  panel.grid.minor.y = element_line(linewidth = 0, color = "transparent"),
  panel.grid.major.x = element_line(linewidth = 0, color = "transparent"),
  panel.grid.major.y = element_line(linewidth = 0, color = "transparent"),
  axis.line.x.bottom = element_line(linewidth = 0.3, color = "black"),
  axis.line.y.left = element_line(linewidth = 0.3, color = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.tag = element_text(size = 11)
)

# Alpha diversity - Faith's PD and Shannon's entropy

faith_alpha_rain_pub_v4 <-
  alpha_values_metav4 %>%
  ggplot(aes(x = SampleSite, y = faith_PD)) + #<<
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
               color = c("#462f25","#750000", "#143e18", "#16474b")
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
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
                              CARAPACE = "Carapace \n(n=15)",
                              "TANK WATER" = "Tank \nwater\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Faith's PD", #<<
    tag = "a)" #<<
  ) + 
  theme_alpha_pub +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

shannon_alpha_rain_pub_v4 <-
  alpha_values_metav4 %>%
  ggplot(aes(x = SampleSite, y = shannon_entropy)) + #<<
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
               color = c("#462f25","#750000", "#143e18", "#16474b")
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
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
                              CARAPACE = "Carapace \n(n=15)",
                              "TANK WATER" = "Tank \nwater\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Shannon's entropy", #<<
    tag = "b)" #<<
  ) + 
  theme_alpha_pub +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

rPCA1_plot_v4_pub <-
ggplot() +
  geom_segment(data = rAitch_pca_v4$data$Species %>% #<<
                 mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                 slice_max(n = 8, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                 mutate(PC1=PC1*3.5, PC2=PC2*3.5) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                 left_join(taxonomy_epiv4), #<<
               aes(x=0, xend=PC1, y=0, yend=PC2),
               linewidth = 0.3,
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30") +
  geom_point(
    data = rAitch_pca_v4$data$Vectors %>% #<<
      left_join(metadata_epiv4), #<<
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
    data = rAitch_pca_v4$data$Vectors %>% #<<
      left_join(metadata_epiv4), #<<
    aes(
      x = PC1, 
      y = PC2, 
      shape = SampleSite
    ),
    size = 3,
    color = "black",
    fill = "transparent"
  ) +
  scale_shape_manual(values = c(23, 21, 22, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                CARAPACE = "Carapace",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#3BB848", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               CARAPACE = "Carapace",
                               "TANK WATER" = "Tank water")) +
  geom_label_repel(data = rAitch_pca_v4$data$Species %>% #<<
                     mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                     slice_max(n = 8, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                     mutate(PC1=PC1*3.5, PC2=PC2*3.5) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                     left_join(taxonomy_epiv4) %>% #<<
                     filter(FeatureID != "af46a1c0a6a2689dfb2631bb616076e1") %>%
                     mutate(Taxa = c("italic('Pseudoalteromonas')~plain('sp.')",
                                     "italic('Shewanella')~plain('sp.')",
                                     "italic('Vibrio')~plain('sp.')",
                                     "italic('Vibrio')~plain('sp.')",
                                     "italic('Vibrio')~plain('sp.')",
                                     "plain('unc.')~italic('Cardiobacteriaceae')",
                                     "plain('unc.')~italic('Marinifilum')")),
                   aes(x = PC1, y = PC2, 
                       label = Taxa),
                   point.padding = 2,
                   alpha = 0.8,
                   colour = "black",
                   size = 3,
                   parse = TRUE) +
  geom_label_repel(data = rAitch_pca_v4$data$Species %>% #<<
                     mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                     slice_max(n = 8, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                     mutate(PC1=PC1*3.5, PC2=PC2*3.5) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                     left_join(taxonomy_epiv4) %>% #<<
                     filter(FeatureID == "af46a1c0a6a2689dfb2631bb616076e1") %>%
                     mutate(Taxa = c("plain('unc.')~italic('Saccharospirillaceae')")),
                   aes(x = PC1, y = PC2, 
                       label = Taxa),
                   point.padding = 3,
                   nudge_x = -0.02,
                   nudge_y = -0.01,
                   alpha = 0.8,
                   colour = "black",
                   size = 3,
                   parse = TRUE) +
  xlab(paste("PC1 (",round(100*rAitch_pca_v4$data$ProportionExplained[1],2),"%)")) + #<<
  ylab(paste("PC2 (",round(100*rAitch_pca_v4$data$ProportionExplained[2],2),"%)")) + #<<
  theme_pcoa1_pub +
  labs(subtitle = NULL,
       tag = "c)")

rPCA1_plot_v4_pub2 <-
  rPCA1_plot_v4_pub + theme(legend.position = c(0.2, 0.8),
                                       legend.box.background = element_rect(linetype = 3, linewidth = 1))

pubfig_epi <- faith_alpha_rain_pub_v4+shannon_alpha_rain_pub_v4 +
  rPCA1_plot_v4_pub2 + 
  plot_layout(design = "
                       ACC
                       BCC
                       ")

ggsave("r_output/16S_epi_v4_merged/fig5_pub.pdf",
       plot = pubfig_epi,
       device=cairo_pdf,
       height = 120,
       width = 190,
       units = "mm")

#----- Session info -----
sessioninfo_v4 <- sessionInfo()
capture.output(sessioninfo, file = "r_output/16S_epi_v4_merged/sessioninfo.txt")
