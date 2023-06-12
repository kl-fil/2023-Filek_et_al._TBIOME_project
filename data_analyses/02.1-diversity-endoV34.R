## Dataviz script for TurtleBiome Holobiont 2023 manuscript
## Processing and visualizing data related to sequencing results for
## endozoic and tank water samples (16S rRNA gene, V34 region)
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
library(pairwiseAdonis) # ???
library(ggpubr) #0.6.0

#set global theme for ggplot
theme_set(theme_light(base_size = 11, base_family = "Arial")) #<< font!

#vectors for colors per sampling site
c("#93624d", "#ff9999", "#89d6dc")
c("#462f25","#750000", "#16474b")

# "#<<" comment indicates the line of code that needs to be changed when using 
# a different dataset for the same code

##---- make directory for visualizations ----
dir.create("r_output/16S_endo_v34/", recursive = TRUE)

##---- Load metadata and data-----
# load metadata in tsv format
metadata_endov34 <- readr::read_tsv("master_metadata/16SendoV34.tsv")

##---- Alpha diversity ----

shannonv34 <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/shannon_vector.qza")
f_pdv34 <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/faith_pd_vector.qza")
evennessv34 <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/evenness_vector.qza")
observed_asvv34 <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/observed_features_vector.qza")

alpha_values_meta <- 
  metadata_endov34 %>%
  left_join(rownames_to_column(shannonv34$data, "SampleID")) %>%
  left_join(rename(f_pdv34$data, "SampleID" = "V1", "faith_PD" = "V2")) %>%
  left_join(rownames_to_column(evennessv34$data, "SampleID")) %>%
  left_join(rownames_to_column(observed_asvv34$data, "SampleID")) %>%
  drop_na(observed_features)

##----- Set theme for alpha div plots
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

##----- Alpha div plots

observed_alpha_rain_v34 <-
alpha_values_meta %>%
  ggplot(aes(x = SampleSite, y = observed_features)) + #<< shannon_entropy|observed_features|faith_PD|pielou_evenness
  ggdist::stat_halfeye(aes(fill = SampleSite),
                       adjust = .5, 
                       width = .6, 
                       .width = 0, 
                       justification = -.3, 
                       point_colour = NA
  ) + 
  geom_boxplot(aes(fill = SampleSite, alpha = 0.5),
               width = .25, 
               outlier.shape = NA,
               color = c("#462f25","#750000", "#16474b")
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
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
                              "TANK WATER" = "Tank water\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Observed features (ASVs)", #<<
    tag = "a)" #<<
  ) + 
  theme_alpha +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

faith_alpha_rain_v34 <-
alpha_values_meta %>%
  ggplot(aes(x = SampleSite, y = faith_PD)) + #<< shannon_entropy|observed_features|faith_PD|pielou_evenness
  ggdist::stat_halfeye(aes(fill = SampleSite),
                       adjust = .5, 
                       width = .6, 
                       .width = 0, 
                       justification = -.3, 
                       point_colour = NA
  ) + 
  geom_boxplot(aes(fill = SampleSite, alpha = 0.5),
               width = .25, 
               outlier.shape = NA,
               color = c("#462f25","#750000", "#16474b")
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
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
                              "TANK WATER" = "Tank water\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Faith's phylogenetic diversity", #<<
    tag = "b)" #<<
  ) + 
  theme_alpha +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

shannon_alpha_rain_v34 <-
alpha_values_meta %>%
  ggplot(aes(x = SampleSite, y = shannon_entropy)) + #<< shannon_entropy|observed_features|faith_PD|pielou_evenness
  ggdist::stat_halfeye(aes(fill = SampleSite),
                       adjust = .5, 
                       width = .6, 
                       .width = 0, 
                       justification = -.3, 
                       point_colour = NA
  ) + 
  geom_boxplot(aes(fill = SampleSite, alpha = 0.5),
               width = .25, 
               outlier.shape = NA,
               color = c("#462f25","#750000", "#16474b")
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
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
                              "TANK WATER" = "Tank water\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Shannon's entropy", #<<
    tag = "c)" #<<
  ) + 
  theme_alpha +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

evenness_alpha_rain_v34 <-
  alpha_values_meta %>%
  ggplot(aes(x = SampleSite, y = pielou_evenness)) + #<< shannon_entropy|observed_features|faith_PD|pielou_evenness
  ggdist::stat_halfeye(aes(fill = SampleSite),
                       adjust = .5, 
                       width = .6, 
                       .width = 0, 
                       justification = -.3, 
                       point_colour = NA
  ) + 
  geom_boxplot(aes(fill = SampleSite, alpha = 0.5),
               width = .25, 
               outlier.shape = NA,
               color = c("#462f25","#750000", "#16474b")
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
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
                              "TANK WATER" = "Tank water\n(n=8)")
  ) +
  labs(
    x = "Sample Site",
    y = "Pielou's evenness", #<<
    tag = "d)" #<<
  ) + 
  theme_alpha +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

rain_plots_v34 <- observed_alpha_rain_v34 + faith_alpha_rain_v34 + shannon_alpha_rain_v34 + evenness_alpha_rain_v34 + plot_layout(ncol = 2)

ggsave(
  "r_output/16S_endo_v34/alpha-div-collected.pdf",
  plot = rain_plots_v34,
  device = cairo_pdf,
  height = 150,
  width = 168,
  units = "mm"
)

#----- Alpha diversity statistics -----

alpha_values_meta$SampleSite <- ordered(alpha_values_meta$SampleSite,
                                        levels = c("CLOACA", "ORAL", "TANK WATER"))

observed_summarystats_v34 <- group_by(alpha_values_meta, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(observed_features, na.rm = TRUE),
    sd = sd(observed_features, na.rm = TRUE),
    median = median(observed_features, na.rm = TRUE),
    IQR = IQR(observed_features, na.rm = TRUE)
  )

faith_summarystats_v34 <- group_by(alpha_values_meta, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(faith_PD, na.rm = TRUE),
    sd = sd(faith_PD, na.rm = TRUE),
    median = median(faith_PD, na.rm = TRUE),
    IQR = IQR(faith_PD, na.rm = TRUE)
  )

shannon_summarystats_v34 <- group_by(alpha_values_meta, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(shannon_entropy, na.rm = TRUE),
    sd = sd(shannon_entropy, na.rm = TRUE),
    median = median(shannon_entropy, na.rm = TRUE),
    IQR = IQR(shannon_entropy, na.rm = TRUE)
  )

evenness_summarystats_v34 <- group_by(alpha_values_meta, SampleSite) %>%
  summarise(
    count = n(),
    mean = mean(pielou_evenness, na.rm = TRUE),
    sd = sd(pielou_evenness, na.rm = TRUE),
    median = median(pielou_evenness, na.rm = TRUE),
    IQR = IQR(pielou_evenness, na.rm = TRUE)
  )

summary_stats_alpha_v34 <- list()
summary_stats_alpha_v34[["Observed features"]] <- observed_summarystats_v34
summary_stats_alpha_v34[["Faith PD"]] <- faith_summarystats_v34
summary_stats_alpha_v34[["Shannon entropy"]] <- shannon_summarystats_v34
summary_stats_alpha_v34[["Pielou evenness"]] <- evenness_summarystats_v34

write_tsv(do.call(cbind, summary_stats_alpha_v34), "r_output/16S_endo_v34/alpha_diversity_summary_stats.tsv")

##---- Kruskall-Wallis test & pairwise comparisons (Wilcox)
krustats_obs_v34 <- list()
krustats_obs_v34[["Kruskall-Wallis test"]] <- kruskal.test(observed_features ~ SampleSite, data = alpha_values_meta)
krustats_obs_v34[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_meta$observed_features, alpha_values_meta$SampleSite,
                                                                            p.adjust.method = "BH")
capture.output(krustats_obs_v34, file = "r_output/16S_endo_v34/alpha_diversity-Kruskall-observed.txt")

krustats_faith_v34 <- list()
krustats_faith_v34[["Kruskall-Wallis test"]] <- kruskal.test(faith_PD ~ SampleSite, data = alpha_values_meta)
krustats_faith_v34[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_meta$faith_PD, alpha_values_meta$SampleSite,
                                                                              p.adjust.method = "BH")
capture.output(krustats_faith_v34, file = "r_output/16S_endo_v34/alpha_diversity-Kruskall-faith_pd.txt")

krustats_shan_v34 <- list()
krustats_shan_v34[["Kruskall-Wallis test"]] <- kruskal.test(shannon_entropy ~ SampleSite, data = alpha_values_meta)
krustats_shan_v34[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_meta$shannon_entropy, alpha_values_meta$SampleSite,
                                                                             p.adjust.method = "BH")
capture.output(krustats_shan_v34, file = "r_output/16S_endo_v34/alpha_diversity-Kruskall-shannon.txt")

krustats_even_v34 <- list()
krustats_even_v34[["Kruskall-Wallis test"]] <- kruskal.test(pielou_evenness ~ SampleSite, data = alpha_values_meta)
krustats_even_v34[["Pairwise comparisons (Wilcox)"]] <- pairwise.wilcox.test(alpha_values_meta$pielou_evenness, alpha_values_meta$SampleSite,
                                                                             p.adjust.method = "BH")
capture.output(krustats_even_v34, file = "r_output/16S_endo_v34/alpha_diversity-Kruskall-evenness.txt")

##----- Kruskal-Wallis test per sample sites -----

# Cloacal samples
kruskall_cloaca <- list()
kruskall_cloaca[["HospStatus"]] <- kruskal.test(faith_PD ~ HospStatus, 
                                                data = alpha_values_meta %>%
                                                  filter(SampleSite == "CLOACA"))
kruskall_cloaca[["TurtleSex"]] <- kruskal.test(faith_PD ~ TurtleSex, 
                                               data = alpha_values_meta %>% 
                                                 filter(SampleSite == "CLOACA"))
kruskall_cloaca[["DetSex"]] <- kruskal.test(faith_PD ~ DetSex, 
                                            data = alpha_values_meta %>%
                                              filter(SampleSite == "CLOACA"))
kruskall_cloaca[["AgeRange"]] <- kruskal.test(faith_PD ~ AgeRange, 
                                              data = alpha_values_meta %>%
                                                filter(SampleSite == "CLOACA"))
kruskall_cloaca[["Age2"]] <- kruskal.test(faith_PD ~ Age2, 
                                          data = alpha_values_meta %>%
                                            filter(SampleSite == "CLOACA"))
capture.output(kruskall_cloaca, file = "r_output/16S_endo_v34/alpha_diversity-Kruskall-cloaca.txt")

# Oral samples
kruskall_oral <- list()
kruskall_oral[["HospStatus"]] <- kruskal.test(faith_PD ~ HospStatus, 
                                              data = alpha_values_meta %>%
                                                filter(SampleSite == "ORAL"))
kruskall_oral[["TurtleSex"]] <- kruskal.test(faith_PD ~ TurtleSex, 
                                             data = alpha_values_meta %>% 
                                               filter(SampleSite == "ORAL"))
kruskall_oral[["DetSex"]] <- kruskal.test(faith_PD ~ DetSex, 
                                          data = alpha_values_meta %>%
                                            filter(SampleSite == "ORAL"))
kruskall_oral[["AgeRange"]] <- kruskal.test(faith_PD ~ AgeRange, 
                                            data = alpha_values_meta %>%
                                              filter(SampleSite == "ORAL"))
kruskall_oral[["Age2"]] <- kruskal.test(faith_PD ~ Age2, 
                                        data = alpha_values_meta %>%
                                          filter(SampleSite == "ORAL"))
capture.output(kruskall_oral, file = "r_output/16S_endo_v34/alpha_diversity-Kruskall-oral.txt")

# Tank water samples
kruskall_wtr <- list()
kruskall_wtr[["HospStatus"]] <- kruskal.test(observed_features ~ HospStatus, 
                                             data = alpha_values_meta %>%
                                               filter(SampleSite == "TANK WATER"))
kruskall_wtr[["TurtleSex"]] <- kruskal.test(observed_features ~ TurtleSex, 
                                            data = alpha_values_meta %>% 
                                              filter(SampleSite == "TANK WATER"))
kruskall_wtr[["DetSex"]] <- kruskal.test(observed_features ~ DetSex, 
                                         data = alpha_values_meta %>%
                                           filter(SampleSite == "TANK WATER"))
kruskall_wtr[["AgeRange"]] <- kruskal.test(observed_features ~ AgeRange, 
                                           data = alpha_values_meta %>%
                                             filter(SampleSite == "TANK WATER"))
kruskall_wtr[["Age2"]] <- kruskal.test(observed_features ~ Age2, 
                                       data = alpha_values_meta %>%
                                         filter(SampleSite == "TANK WATER"))
capture.output(kruskall_wtr, file = "r_output/16S_endo_v34/alpha_diversity-Kruskall-tank_water.txt")

## Pairwise for interesting categories
# cloaca - Faith's PD, observed features, Shannon's entropy - sex and age

alpha_values_meta_cloaca <- alpha_values_meta %>%
  filter(SampleSite == "CLOACA")
pwilcox_cloaca <- list()
pwilcox_cloaca[["faithPD TurtleSex"]] <-
  pairwise.wilcox.test(alpha_values_meta_cloaca$faith_PD, alpha_values_meta_cloaca$TurtleSex, p.adjust.method = "BH")
pwilcox_cloaca[["faithPD DetSex"]] <-
  pairwise.wilcox.test(alpha_values_meta_cloaca$faith_PD, alpha_values_meta_cloaca$DetSex, p.adjust.method = "BH")
pwilcox_cloaca[["faithPD AgeRange"]] <- 
  pairwise.wilcox.test(alpha_values_meta_cloaca$faith_PD, alpha_values_meta_cloaca$AgeRange, p.adjust.method = "BH")
pwilcox_cloaca[["faithPD Age2"]] <- 
  pairwise.wilcox.test(alpha_values_meta_cloaca$faith_PD, alpha_values_meta_cloaca$Age2, p.adjust.method = "BH")
pwilcox_cloaca[["shannon TurtleSex"]] <-
  pairwise.wilcox.test(alpha_values_meta_cloaca$shannon_entropy, alpha_values_meta_cloaca$TurtleSex, p.adjust.method = "BH")
pwilcox_cloaca[["shannon DetSex"]] <-
  pairwise.wilcox.test(alpha_values_meta_cloaca$shannon_entropy, alpha_values_meta_cloaca$DetSex, p.adjust.method = "BH")
pwilcox_cloaca[["shannon AgeRange"]] <- 
  pairwise.wilcox.test(alpha_values_meta_cloaca$shannon_entropy, alpha_values_meta_cloaca$AgeRange, p.adjust.method = "BH")
pwilcox_cloaca[["shannon Age2"]] <- 
  pairwise.wilcox.test(alpha_values_meta_cloaca$shannon_entropy, alpha_values_meta_cloaca$Age2, p.adjust.method = "BH")
pwilcox_cloaca[["observed TurtleSex"]] <-
  pairwise.wilcox.test(alpha_values_meta_cloaca$observed_features, alpha_values_meta_cloaca$TurtleSex, p.adjust.method = "BH")
pwilcox_cloaca[["observed DetSex"]] <-
  pairwise.wilcox.test(alpha_values_meta_cloaca$observed_features, alpha_values_meta_cloaca$DetSex, p.adjust.method = "BH")
pwilcox_cloaca[["observed AgeRange"]] <- 
  pairwise.wilcox.test(alpha_values_meta_cloaca$observed_features, alpha_values_meta_cloaca$AgeRange, p.adjust.method = "BH")
pwilcox_cloaca[["observed Age2"]] <- 
  pairwise.wilcox.test(alpha_values_meta_cloaca$observed_features, alpha_values_meta_cloaca$Age2, p.adjust.method = "BH")
capture.output(pwilcox_cloaca, file = "r_output/16S_endo_v34/alpha_diversity-pairwise-cloaca.txt")

#----- Correlation between alpha diversity and CCL per sample site -----

pears_faith_clo_v34 <- alpha_values_meta %>%
  filter(SampleSite %in% c("CLOACA")) %>%
  ggplot(aes(x = CCL, y = faith_PD)) +
  geom_smooth(data = alpha_values_meta %>%
                filter(SampleSite== "CLOACA"), method=lm , fill = "gray80", color="#462f25", se=TRUE) +
  stat_cor(data = alpha_values_meta %>%
             filter(SampleSite== "CLOACA"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#93624d", shape = 21, size = 2) + 
  labs(title = "Cloaca",
       x = "CCL (cm)",
       y = "Faith's PD",
       tag = "a)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank())

pears_faith_orl_v34 <- alpha_values_meta %>%
  filter(SampleSite %in% c("ORAL")) %>%
  ggplot(aes(x = CCL, y = faith_PD)) +
  geom_smooth(data = alpha_values_meta %>%
                filter(SampleSite== "ORAL"), method=lm , fill = "gray80", color="#750000", se=TRUE) +
  stat_cor(data = alpha_values_meta %>%
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

pears_faith_wtr_v34 <- alpha_values_meta %>%
  filter(SampleSite %in% c("TANK WATER")) %>%
  ggplot(aes(x = CCL, y = faith_PD)) +
  geom_smooth(data = alpha_values_meta %>%
                filter(SampleSite== "TANK WATER"), method=lm , fill = "gray80", color="#16474b", se=TRUE) +
  stat_cor(data = alpha_values_meta %>%
             filter(SampleSite== "TANK WATER"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#89d6dc", shape = 21, size = 2) + 
  labs(title = "Tank Water",
       x = "CCL (cm)",
       y = "Faith's PD",
       tag = "c)") +
  theme_alpha+ 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pears_obs_clo_v34 <- alpha_values_meta %>%
  filter(SampleSite %in% c("CLOACA")) %>%
  ggplot(aes(x = CCL, y = observed_features)) +
  geom_smooth(data = alpha_values_meta %>%
                filter(SampleSite== "CLOACA"), method=lm , fill = "gray80", color="#462f25", se=TRUE) +
  stat_cor(data = alpha_values_meta %>%
             filter(SampleSite== "CLOACA"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#93624d", shape = 21, size = 2) + 
  labs(title = "Cloaca",
       x = "CCL (cm)",
       y = "ASVs",
       tag = "d)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank())

pears_obs_orl_v34 <- alpha_values_meta %>%
  filter(SampleSite %in% c("ORAL")) %>%
  ggplot(aes(x = CCL, y = observed_features)) +
  geom_smooth(data = alpha_values_meta %>%
                filter(SampleSite== "ORAL"), method=lm , fill = "gray80", color="#750000", se=TRUE) +
  stat_cor(data = alpha_values_meta %>%
             filter(SampleSite== "ORAL"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#ff9999", shape = 21, size = 2) + 
  labs(title = "Oral",
       x = "CCL (cm)",
       y = "ASVs",
       tag = "e)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pears_obs_wtr_v34 <- alpha_values_meta %>%
  filter(SampleSite %in% c("TANK WATER")) %>%
  ggplot(aes(x = CCL, y = observed_features)) +
  geom_smooth(data = alpha_values_meta %>%
                filter(SampleSite== "TANK WATER"), method=lm , fill = "gray80", color="#16474b", se=TRUE) +
  stat_cor(data = alpha_values_meta %>%
             filter(SampleSite== "TANK WATER"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#89d6dc", shape = 21, size = 2) + 
  labs(title = "Tank Water",
       x = "CCL (cm)",
       y = "ASVs",
       tag = "f)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pears_sha_clo_v34 <- alpha_values_meta %>%
  filter(SampleSite %in% c("CLOACA")) %>%
  ggplot(aes(x = CCL, y = shannon_entropy)) +
  geom_smooth(data = alpha_values_meta %>%
                filter(SampleSite== "CLOACA"), method=lm , fill = "gray80", color="#462f25", se=TRUE) +
  stat_cor(data = alpha_values_meta %>%
             filter(SampleSite== "CLOACA"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#93624d", shape = 21, size = 2) + 
  labs(title = "Cloaca",
       x = "CCL (cm)",
       y = "Shannon's",
       tag = "g)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line())

pears_sha_orl_v34 <- alpha_values_meta %>%
  filter(SampleSite %in% c("ORAL")) %>%
  ggplot(aes(x = CCL, y = shannon_entropy)) +
  geom_smooth(data = alpha_values_meta %>%
                filter(SampleSite== "ORAL"), method=lm , fill = "gray80", color="#750000", se=TRUE) +
  stat_cor(data = alpha_values_meta %>%
             filter(SampleSite== "ORAL"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#ff9999", shape = 21, size = 2) + 
  labs(title = "Oral",
       x = "CCL (cm)",
       y = "Shannon's",
       tag = "h)") +
  theme_alpha + 
  theme(panel.grid.major.x = element_line(),
        axis.title.y = element_blank())

pears_sha_wtr_v34 <- alpha_values_meta %>%
  filter(SampleSite %in% c("TANK WATER")) %>%
  ggplot(aes(x = CCL, y = shannon_entropy)) +
  geom_smooth(data = alpha_values_meta %>%
                filter(SampleSite== "TANK WATER"), method=lm , fill = "gray80", color="#16474b", se=TRUE) +
  stat_cor(data = alpha_values_meta %>%
             filter(SampleSite== "TANK WATER"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#89d6dc", shape = 21, size = 2) + 
  labs(title = "Tank Water",
       x = "CCL (cm)",
       y = "Shannon's",
       tag = "i)") +
  theme_alpha+ 
  theme(panel.grid.major.x = element_line(),
        axis.title.y = element_blank())

pears_alpha_v34<-pears_faith_clo_v34+pears_faith_orl_v34+pears_faith_wtr_v34+
  pears_obs_clo_v34+pears_obs_orl_v34+pears_obs_wtr_v34+
  pears_sha_clo_v34+pears_sha_orl_v34+pears_sha_wtr_v34+plot_layout(ncol = 3)

ggsave(
  "r_output/16S_endo_v34/pears-alpha-collected.pdf",
  plot = pears_alpha_v34,
  device = cairo_pdf,
  height = 170,
  width = 200,
  units = "mm"
)

pearson_corr_cloaca <- list()
pearson_corr_cloaca[["Observed features"]] <- with(alpha_values_meta %>%
                                                     filter(SampleSite == "CLOACA"),cor.test(CCL,observed_features, method = "pearson"))
pearson_corr_cloaca[["Faith's PD"]] <- with(alpha_values_meta %>%
                                              filter(SampleSite == "CLOACA"),cor.test(CCL,faith_PD, method = "pearson"))
pearson_corr_cloaca[["Shannon's entropy"]] <- with(alpha_values_meta %>%
                                                     filter(SampleSite == "CLOACA"),cor.test(CCL,shannon_entropy, method = "pearson"))
pearson_corr_cloaca[["Pielou's evenness"]] <- with(alpha_values_meta %>%
                                                     filter(SampleSite == "CLOACA"),cor.test(CCL,pielou_evenness, method = "pearson"))
capture.output(pearson_corr_cloaca, file = "r_output/16S_endo_v34/alpha_diversity-pearson-cloaca.txt")


pearson_corr_oral <- list()
pearson_corr_oral[["Observed features"]] <- with(alpha_values_meta %>%
                                                   filter(SampleSite == "ORAL"),cor.test(CCL,observed_features, method = "pearson"))
pearson_corr_oral[["Faith's PD"]] <- with(alpha_values_meta %>%
                                              filter(SampleSite == "ORAL"),cor.test(CCL,faith_PD, method = "pearson"))
pearson_corr_oral[["Shannon's entropy"]] <- with(alpha_values_meta %>%
                                                     filter(SampleSite == "ORAL"),cor.test(CCL,shannon_entropy, method = "pearson"))
pearson_corr_oral[["Pielou's evenness"]] <- with(alpha_values_meta %>%
                                                     filter(SampleSite == "ORAL"),cor.test(CCL,pielou_evenness, method = "pearson"))
capture.output(pearson_corr_oral, file = "r_output/16S_endo_v34/alpha_diversity-pearson-oral.txt")

pearson_corr_tw <- list()
pearson_corr_tw[["Observed features"]] <- with(alpha_values_meta %>%
                                                   filter(SampleSite == "TANK WATER"),cor.test(CCL,observed_features, method = "pearson"))
pearson_corr_tw[["Faith's PD"]] <- with(alpha_values_meta %>%
                                            filter(SampleSite == "TANK WATER"),cor.test(CCL,faith_PD, method = "pearson"))
pearson_corr_tw[["Shannon's entropy"]] <- with(alpha_values_meta %>%
                                                   filter(SampleSite == "TANK WATER"),cor.test(CCL,shannon_entropy, method = "pearson"))
pearson_corr_tw[["Pielou's evenness"]] <- with(alpha_values_meta %>%
                                                   filter(SampleSite == "TANK WATER"),cor.test(CCL,pielou_evenness, method = "pearson"))
capture.output(pearson_corr_tw, file = "r_output/16S_endo_v34/alpha_diversity-pearson-tank_water.txt")

##----- Import data for PCoA plots -----

# load taxonomies
taxonomy_endov34 <- read_qza("qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza")

##----- Load data from qiime2 qza to R for PCoA/PCA plots
# load qiime pca and pcoa results.qza (in_data_16 folder)
# 16S data loading (from core-metrics-results-merged_16S-0-with-phyla-no-mitochondria-no-chloroplast-filtered-phylogeny )
# at SeqTry3 folder of 16S analyses
bray_pcoa <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/bray_curtis_pcoa_results.qza")
unifrac_pcoa <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/unweighted_unifrac_pcoa_results.qza")
w_unifrac_pcoa <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/weighted_unifrac_pcoa_results.qza")
jaccard_pcoa <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/jaccard_pcoa_results.qza")
rAitch_pca <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-ordination-28000.qza")

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

bray1_plot_v34 <-
bray_pcoa$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_endov34) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(x = PC1, y = PC2, fill = SampleSite, color = SampleSite), 
               type="norm",
               level = 0.75,
               alpha = 0.1, 
               geom = "polygon",
               show.legend = FALSE) +
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
  # geom_mark_ellipse(expand = 0.02, aes(color = SampleSite)) +
  # geom_text_repel(aes(label = SampleID), size = 1.75) +
  scale_shape_manual(values = c(23, 21, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")) +
  scale_color_manual(name = "Sample site:", 
                     values = c("#93624d", "#ff9999", "#89d6dc"),
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  #scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  xlab(paste("PC1 (", round(100*bray_pcoa$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*bray_pcoa$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Bray-Curtis",
       tag = "a)")

#jacc1_plot_v34 <-
jaccard_pcoa$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_endov34) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(x = PC1, y = PC2, fill = SampleSite, color = SampleSite), 
               type="norm",
               level = 0.75,
               alpha = 0.1, 
               geom = "polygon",
               show.legend = FALSE) +
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
  # geom_mark_ellipse(expand = 0.02, aes(color = Genus)) +
  # geom_text_repel(aes(label = SampleID), size = 1.75) +
  scale_shape_manual(values = c(23, 21, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")) +
  scale_color_manual(name = "Sample site:", 
                     values = c("#93624d", "#ff9999", "#89d6dc"),
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  #scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  xlab(paste("PC1 (", round(100*jaccard_pcoa$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*jaccard_pcoa$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Jaccard",
       tag = "b)")

#unwunifrac1_plot_v34 <-
unifrac_pcoa$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_endov34) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(x = PC1, y = PC2, fill = SampleSite, color = SampleSite), 
               type="norm",
               level = 0.75,
               alpha = 0.1, 
               geom = "polygon",
               show.legend = FALSE) +
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
  # geom_mark_ellipse(expand = 0.02, aes(color = Genus)) +
  # geom_text_repel(aes(label = SampleID), size = 1.75) +
  scale_shape_manual(values = c(23, 21, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")) +
  scale_color_manual(name = "Sample site:", 
                     values = c("#93624d", "#ff9999", "#89d6dc"),
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  #scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  xlab(paste("PC1 (", round(100*unifrac_pcoa$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*unifrac_pcoa$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Unweighted UniFrac",
       tag = "c)")

#wunifrac1_plot_v34 <-
w_unifrac_pcoa$data$Vectors %>% #<<
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata_endov34) %>%  #<<
  ggplot(aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(x = PC1, y = PC2, fill = SampleSite, color = SampleSite), 
               type="norm",
               level = 0.75,
               alpha = 0.1, 
               geom = "polygon",
               show.legend = FALSE) +
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
  # geom_mark_ellipse(expand = 0.02, aes(color = Genus)) +
  # geom_text_repel(aes(label = SampleID), size = 1.5) +
  scale_shape_manual(values = c(23, 21, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")) +
  scale_color_manual(name = "Sample site:", 
                     values = c("#93624d", "#ff9999", "#89d6dc"),
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  guides(fill = guide_legend(byrow = TRUE, direction = "vertical"),
         shape = guide_legend(byrow = TRUE, direction = "vertical")) +
  #scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  xlab(paste("PC1 (", round(100*w_unifrac_pcoa$data$ProportionExplained[1],2),"%)")) +  #<<
  ylab(paste("PC2 (", round(100*w_unifrac_pcoa$data$ProportionExplained[2],2),"%)")) +  #<<
  theme_pcoa1 + #<< set theme
  labs(subtitle = "Weighted UniFrac",
       tag = "d)")


##---- PCA plots
# rename column in taxonomy so we can join it to our data
taxonomy_endov34 <- rename(taxonomy_endov34$data, FeatureID = Feature.ID)

rPCA1_plot_v34 <- 
ggplot() +
  # stat_ellipse(data = rAitch_pca$data$Vectors %>% #<<
  #                left_join(metadata_endov34),
  #              aes(x = PC1, y = PC2, fill = SampleSite, color = SampleSite),
  #              type="norm",
  #              level = 0.75,
  #              alpha = 0.1,
  #              geom = "polygon",
  #              show.legend = FALSE) +
  geom_segment(data = rAitch_pca$data$Species %>% #<<
                 mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                 slice_max(n = 5, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                 mutate(PC1=PC1*3, PC2=PC2*3) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                 left_join(taxonomy_endov34), #<<
               aes(x=0, xend=PC1, y=0, yend=PC2),
               linewidth = 0.3,
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30") +
  geom_point(
    data = rAitch_pca$data$Vectors %>% #<<
      left_join(metadata_endov34), #<<
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
    data = rAitch_pca$data$Vectors %>% #<<
      left_join(metadata_endov34), #<<
    aes(
      x = PC1, 
      y = PC2, 
      shape = SampleSite
    ),
    size = 3,
    color = "gray30",
    fill = "transparent"
  ) +
  # geom_text_repel(data = rAitch_pca$data$Vectors %>%
  #                   left_join(metadata_endov34),
  #                 aes(x = PC1, y = PC2, label = TurtleName), size = 1.75) +
  scale_shape_manual(values = c(23, 21, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")) +
  scale_color_manual(name = "Sample site:", 
                     values = c("#93624d", "#ff9999", "#89d6dc"),
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  geom_label_repel(data = rAitch_pca$data$Species %>% #<<
                     mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                     slice_max(n = 5, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                     mutate(PC1=PC1*3, PC2=PC2*3) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                     left_join(taxonomy_endov34) %>% #<<
                     mutate(Taxa = c("plain('order Oceanospirillales')~italic('')",
                                     "plain('')~italic('Shewanella algae')",
                                     "plain('NS3a marine group')~italic('')",
                                     "italic('Vibrio')~plain('sp.')",
                                     "plain('class Gammaproteobacteria')~italic('')")),
                   aes(x = PC1, y = PC2, 
                       label = Taxa),
                   point.padding = 2,
                   alpha = 0.8,
                   colour = "black",
                   size = 3,
                   parse = TRUE) +
  # guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
  #        fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  xlab(paste("PC1 (",round(100*rAitch_pca$data$ProportionExplained[1],2),"%)")) + #<<
  ylab(paste("PC2 (",round(100*rAitch_pca$data$ProportionExplained[2],2),"%)")) + #<<
  theme_pcoa1 +
  labs(subtitle = NULL,
       tag = "e)")

PCoA_plots_collected_v34 <- ((bray1_plot_v34 + jacc1_plot_v34) / (unwunifrac1_plot_v34 + wunifrac1_plot_v34)) + plot_layout(guides = 'collect')

rPCA1_plot_v34

ggsave(
  filename = "r_output/16S_endo_v34/rpca_plot.pdf",
  #plot = PCoA_plots_collected_v34,
  device = cairo_pdf,
  height = 140,
  width = 160,
  units = "mm"
)

ggsave(filename = "r_output/16S_endo_v34/rpca_plot_turtlenames.pdf",
       #plot = rPCA1_plot_v34,
       device = cairo_pdf,
       height = 140,
       width = 160,
       units = "mm"
)

# Pearson and beta vs CCL
pclo5 <- ggplot(raitch_pcoa_clo$data$Vectors %>%
                  left_join(alpha_values_meta),
                aes(
                  x = CCL, 
                  y = PC1
                )) +
  geom_smooth(method=lm , fill = "gray80", color="black", se=TRUE) +
  stat_cor(method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(aes(fill=AgeRange, size=CCL),
             #size = 4,
             alpha = 0.8,
             color = "black",
             shape = 21
  ) +
  guides(fill = guide_legend(override.aes=list(shape=21, size = 4)),
         shape = guide_legend(override.aes=list(size = 4))) +
  labs(title = "Cloaca - r. Aitchison~CCL")

(pclo1 + pclo2+pclo3+pclo4+pclo5) + plot_layout(guides = "collect", ncol = 2)


ggsave(filename = "r_output/16S_endo_v34/pearson_pc1s_clo.pdf",
       #plot = rPCA1_plot,
       device = cairo_pdf,
       height = 200,
       width = 250,
       units = "mm"
)

##----- ADONIS -----
dir.create(path = "r_output/16S_endo_v34/adonis_stats", recursive = TRUE)

colnames(metadata_endov34)

bray_dist_v34 <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/bray_curtis_distance_matrix.qza")
jacc_dist_v34 <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/jaccard_distance_matrix.qza")
uunif_dist_v34 <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/unweighted_unifrac_distance_matrix.qza")
wunif_dist_v34 <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/weighted_unifrac_distance_matrix.qza")
raitc_dist_v34 <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix-28000.qza")

bray_dist_data_v34 <- bray_dist_v34$data
jacc_dist_data_v34 <- jacc_dist_v34$data
uunif_dist_data_v34 <- uunif_dist_v34$data
wunif_dist_data_v34 <- wunif_dist_v34$data
raitc_dist_data_v34 <- raitc_dist_v34$data
##

adonis_metadata_v34 <- filter(metadata_endov34, !SampleID %in% removed_samples_v34)
removed_samples_v34 <- c("16SNEGCTRL",
                     "16S0094O",
                     "16S0113C",
                     "16S0118O",
                     "16S0118W",
                     "16S0119C",
                     "16S0092O",
                     "16S0064O")
# distances and metadata need to be sorted in the same way (in this cae alphabetically)
adonis_metadata_v34 <- adonis_metadata_v34[order(adonis_metadata_v34$SampleID),]

set.seed(123)
adonis2(formula = bray_dist_data_v34~SampleSite, data = adonis_metadata_v34, permutations = 999, by = "terms")
adonis2(formula = jacc_dist_data_v34~SampleSite, data = adonis_metadata_v34, permutations = 999, by = "terms")
adonis2(formula = uunif_dist_data_v34~SampleSite, data = adonis_metadata_v34, permutations = 999, by = "terms")
adonis2(formula = wunif_dist_data_v34~SampleSite, data = adonis_metadata_v34, permutations = 999, by = "terms")
adonis2(formula = raitc_dist_data_v34~SampleSite, data = adonis_metadata_v34, permutations = 999, by = "terms")
pairwiseAdonis::pairwise.adonis2(bray_dist_data_v34~SampleSite, data = adonis_metadata_v34, permutations = 999)
pairwiseAdonis::pairwise.adonis2(jacc_dist_data_v34~SampleSite, data = adonis_metadata_v34, permutations = 999)
pairwiseAdonis::pairwise.adonis2(uunif_dist_data_v34~SampleSite, data = adonis_metadata_v34, permutations = 999)
pairwiseAdonis::pairwise.adonis2(wunif_dist_data_v34~SampleSite, data = adonis_metadata_v34, permutations = 999)
pairwiseAdonis::pairwise.adonis2(raitc_dist_data_v34~SampleSite, data = adonis_metadata_v34, permutations = 999)

#Adonis2 cloaca
bray_dist_clo <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/bray_curtis_distance_matrix.qza")
jacc_dist_clo <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/jaccard_distance_matrix.qza")
uunif_dist_clo <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/unweighted_unifrac_distance_matrix.qza")
wunif_dist_clo <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/weighted_unifrac_distance_matrix.qza")
raitc_dist_clo <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix-28000.qza")

bray_dist_data_clo <- bray_dist_clo$data
jacc_dist_data_clo <- jacc_dist_clo$data
uunif_dist_data_clo <- uunif_dist_clo$data
wunif_dist_data_clo <- wunif_dist_clo$data
raitc_dist_data_clo <- raitc_dist_clo$data

removed_samples <- c("16SNEGCTRL",
                     "16S0094O",
                     "16S0113C",
                     "16S0118O",
                     "16S0118W",
                     "16S0119C",
                     "16S0092O",
                     "16S0064O")
adonis_metadata_clo <- metadata_endov34 %>%
  filter(!SampleID %in% removed_samples) %>%
  filter(SampleSite == "CLOACA")

# distances and metadata need to be sorted in the same way (in this cae alphabetically)
adonis_metadata_clo <- adonis_metadata_clo[order(adonis_metadata_clo$SampleID),]

adonis2(formula = bray_dist_data_clo~TurtleSex, data = adonis_metadata_clo, permutations = 999, by = "terms")
adonis2(formula = jacc_dist_data_clo~TurtleSex, data = adonis_metadata_clo, permutations = 999, by = "terms")
adonis2(formula = uunif_dist_data_clo~TurtleSex, data = adonis_metadata_clo, permutations = 999, by = "terms")
adonis2(formula = wunif_dist_data_clo~TurtleSex, data = adonis_metadata_clo, permutations = 999, by = "terms")
#adonis2(formula = raitc_dist_data_clo~HospStatus, data = adonis_metadata_clo, permutations = 999, by = "terms") #error, lost sample

#aitch from raw counts D28000
clo_table <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal.qza")
clo_counts <- t(clo_table$data)
clo_raitch_dist <- vegdist(clo_counts, method = "robust.aitchison")

adonis2(formula = clo_raitch_dist~TurtleSex*AgeRange, data = adonis_metadata_clo, permutations = 999, by = "terms")
pairwiseAdonis::pairwise.adonis2(bray_dist_data_clo~TurtleSex*AgeRange, data = adonis_metadata_clo, permutations = 999)

pairwiseAdonis::pairwise.adonis2(bray_dist_data_clo~SampleSite, data = adonis_metadata_clo, permutations = 999)
pairwiseAdonis::pairwise.adonis2(jacc_dist_data_clo~SampleSite, data = adonis_metadata_clo, permutations = 999)
pairwiseAdonis::pairwise.adonis2(uunif_dist_data_clo~SampleSite, data = adonis_metadata_clo, permutations = 999)
pairwiseAdonis::pairwise.adonis2(wunif_dist_data_clo~SampleSite, data = adonis_metadata_clo, permutations = 999)
pairwiseAdonis::pairwise.adonis2(raitc_dist_data_clo~SampleSite, data = adonis_metadata_clo, permutations = 999)

#Adonis2 oral
bray_dist_orl <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/bray_curtis_distance_matrix.qza")
jacc_dist_orl <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/jaccard_distance_matrix.qza")
uunif_dist_orl <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/unweighted_unifrac_distance_matrix.qza")
wunif_dist_orl <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/weighted_unifrac_distance_matrix.qza")
raitc_dist_orl <- read_qza("qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix-28000.qza")

bray_dist_data_orl <- bray_dist_orl$data
jacc_dist_data_orl <- jacc_dist_orl$data
uunif_dist_data_orl <- uunif_dist_orl$data
wunif_dist_data_orl <- wunif_dist_orl$data
raitc_dist_data_orl <- raitc_dist_orl$data

removed_samples <- c("16SNEGCTRL",
                     "16S0094O",
                     "16S0113C",
                     "16S0118O",
                     "16S0118W",
                     "16S0119C",
                     "16S0092O",
                     "16S0064O")
adonis_metadata_orl <- metadata_endov34 %>%
  filter(!SampleID %in% removed_samples) %>%
  filter(SampleSite == "ORAL")

# distances and metadata need to be sorted in the same way (in this cae alphabetically)
adonis_metadata_orl <- adonis_metadata_orl[order(adonis_metadata_orl$SampleID),]

adonis2(formula = bray_dist_data_orl~HospStatus, data = adonis_metadata_orl, permutations = 999, by = "terms")
adonis2(formula = jacc_dist_data_orl~HospStatus, data = adonis_metadata_orl, permutations = 999, by = "terms")
adonis2(formula = uunif_dist_data_orl~HospStatus, data = adonis_metadata_orl, permutations = 999, by = "terms")
adonis2(formula = wunif_dist_data_orl~HospStatus, data = adonis_metadata_orl, permutations = 999, by = "terms")
adonis2(formula = raitc_dist_data_orl~HospStatus, data = adonis_metadata_orl, permutations = 999, by = "terms") 

pairwiseAdonis::pairwise.adonis2(bray_dist_data_orl~HospStatus, data = adonis_metadata_orl, permutations = 999)
pairwiseAdonis::pairwise.adonis2(jacc_dist_data_orl~HospStatus, data = adonis_metadata_orl, permutations = 999)
pairwiseAdonis::pairwise.adonis2(uunif_dist_data_orl~HospStatus, data = adonis_metadata_orl, permutations = 999)
pairwiseAdonis::pairwise.adonis2(wunif_dist_data_orl~HospStatus, data = adonis_metadata_orl, permutations = 999)
pairwiseAdonis::pairwise.adonis2(raitc_dist_data_orl~HospStatus, data = adonis_metadata_orl, permutations = 999)


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
faith_alpha_rain_pub <-
alpha_values_meta %>%
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
               color = c("#462f25","#750000", "#16474b")
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
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
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

shannon_alpha_rain_pub <-
  alpha_values_meta %>%
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
               color = c("#462f25","#750000", "#16474b")
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
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")
  ) +
  scale_x_discrete(labels = c(CLOACA = "Cloacal\n(n=19)",
                              ORAL = "Oral\n(n=16)",
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

# Pearson correlation for cloacal samples and alpha diversity
pears_faith_clo_pub <- alpha_values_meta %>%
  filter(SampleSite %in% c("CLOACA")) %>%
  ggplot(aes(x = CCL, y = faith_PD)) +
  geom_smooth(data = alpha_values_meta %>%
                filter(SampleSite== "CLOACA"), method=lm , fill = "gray80", color="#462f25", se=TRUE) +
  stat_cor(data = alpha_values_meta %>%
             filter(SampleSite== "CLOACA"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_point(fill = "#93624d", shape = 21, size = 2) + 
  labs(x = "CCL (cm)",
       y = "Faith's PD",
       tag = "c)") +
  theme_alpha_pub + 
  theme(panel.grid.major.x = element_line())

pears_sha_clo_pub2 <- alpha_values_meta %>%
  filter(SampleSite %in% c("CLOACA")) %>%
  ggplot(aes(x = CCL, y = shannon_entropy)) +
  geom_smooth(data = alpha_values_meta %>%
                filter(SampleSite== "CLOACA"), method=lm , fill = "gray80", color="#462f25", se=TRUE) +
  stat_cor(data = alpha_values_meta %>%
             filter(SampleSite== "CLOACA"), method="pearson", p.accuracy = 0.001, r.accuracy = 0.01,
           label.y = 3.5) +
  geom_point(fill = "#93624d", shape = 21, size = 2) + 
  labs(x = "CCL (cm)",
       y = "Shannon's entropy",
       tag = "c)") +
  theme_alpha_pub + 
  theme(panel.grid.major.x = element_line())

# robust Aitchison (DEICODE) PCA plot
rPCA1_plot_pub2 <- 
  ggplot() +
  geom_segment(data = rAitch_pca$data$Species %>% #<<
                 mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                 slice_max(n = 5, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                 mutate(PC1=PC1*3, PC2=PC2*3) %>% #scale arrow linearly (look at emperor qiime vizualisation)
                 left_join(taxonomy_endov34), #<<
               aes(x=0, xend=PC1, y=0, yend=PC2),
               linewidth = 0.3,
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30") +
  geom_point(
    data = rAitch_pca$data$Vectors %>% #<<
      left_join(metadata_endov34), #<<
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
    data = rAitch_pca$data$Vectors %>% #<<
      left_join(metadata_endov34), #<<
    aes(
      x = PC1, 
      y = PC2, 
      shape = SampleSite
    ),
    size = 3,
    color = "black",
    fill = "transparent"
  ) +
  # geom_text_repel(data = rAitch_pca$data$Vectors %>%
  #                   left_join(metadata_endov34),
  #                 aes(x = PC1, y = PC2, label = TurtleName), size = 1.75) +
  scale_shape_manual(values = c(23, 21, 25), 
                     name = "Sample site:",
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  scale_fill_manual(name = "Sample site:", 
                    values = c("#93624d", "#ff9999", "#89d6dc"),
                    labels = c(CLOACA = "Cloacal",
                               ORAL = "Oral",
                               "TANK WATER" = "Tank water")) +
  scale_color_manual(name = "Sample site:", 
                     values = c("#93624d", "#ff9999", "#89d6dc"),
                     labels = c(CLOACA = "Cloacal",
                                ORAL = "Oral",
                                "TANK WATER" = "Tank water")) +
  geom_label_repel(data = rAitch_pca$data$Species %>% #<<
                     mutate(a = sqrt(PC1^2+PC2^2+PC3^2)) %>% #calculate the distance from the origin
                     slice_max(n = 5, order_by = a) %>% #keep X furthest away points - instead of top_n -- slice_max
                     mutate(PC1=PC1*3, PC2=PC2*3) %>% #scale arrow linearly, i guess it's ok (look at emperor)
                     left_join(taxonomy_endov34) %>% #<<
                     mutate(Taxa = c("plain('order Oceanospirillales')~italic('')",
                                     "plain('')~italic('Shewanella algae')",
                                     "plain('NS3a marine group')~italic('')",
                                     "italic('Vibrio')~plain('sp.')",
                                     "plain('class Gammaproteobacteria')~italic('')")),
                   aes(x = PC1, y = PC2, 
                       label = Taxa),
                   point.padding = 2,
                   alpha = 0.8,
                   colour = "black",
                   size = 3,
                   parse = TRUE) +
  # guides(shape = guide_legend(byrow = TRUE, direction = "vertical", order = 1), 
  #        fill = guide_legend(byrow = TRUE, override.aes=list(shape=21), direction = "vertical", order = 2)) +
  xlab(paste("PC1 (",round(100*rAitch_pca$data$ProportionExplained[1],2),"%)")) + #<<
  ylab(paste("PC2 (",round(100*rAitch_pca$data$ProportionExplained[2],2),"%)")) + #<<
  theme_pcoa1_pub +
  labs(subtitle = NULL,
       tag = "d)")

rPCA1_plot22 <-rPCA1_plot_pub2 + theme(legend.position = c(0.8, 0.8),
                        legend.box.background = element_rect(linetype = 3, linewidth = 1))

faith_alpha_rain_pub+shannon_alpha_rain_pub +
  pears_sha_clo_pub2 + rPCA1_plot22 + 
  plot_layout(design = "
                       ABDDD
                       CCDDD
                       ")

ggsave("r_output/16S_endo_v34/fig2_pub.pdf",
       device=cairo_pdf,
       height = 120,
       width = 220,
       units = "mm")

#----- Session info -----
sessionInfo()