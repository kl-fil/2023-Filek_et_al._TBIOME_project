## Dataviz script for TurtleBiome Holobiont 2023 manuscript
## Processing and visualizing data related to sequencing results for
## endozoic and tank water samples (16S rRNA gene, V34 region)
## Author: Klara Filek

##--------------------load packages--------------------
library(tidyverse)
library(ggplot2)
library(qiime2R)
library(ANCOMBC)
library(tidyverse)
library(phyloseq)
library(ggtext)
library(qiime2R)

#set global theme for ggplot
theme_set(theme_light(base_size = 11, base_family = "Arial")) #<< font!

# "#<<" comment indicates the line of code that needs to be changed when using 
# a different dataset for the same code

##---- make directory for visualizations ----
dir.create("r_output/16S_endo_v34/", recursive = TRUE)

##---- Load metadata and data-----
# load metadata in tsv format
metadata_endov34 <- readr::read_tsv("master_metadata/16SendoV34.tsv")

##---- Import taxonomy and counts (extracted from qiime2 barplots) ----

# generate csv files with abundances from mono-source/taxonomy/taxa-bar-plots.qzv 
# save each level as csv in the mono-source/taxonomy folder
# import csv files taxa as columns, samples as rows
taxadata_2_endov34 <- read.csv("qiime2_output/16S_endo_v34_2021/taxonomy/level-2.csv", row.names = 1) #row names as sample id
taxadata_3_endov34 <- read.csv("qiime2_output/16S_endo_v34_2021/taxonomy/level-3.csv", row.names = 1) #row names as sample id
taxadata_4_endov34 <- read.csv("qiime2_output/16S_endo_v34_2021/taxonomy/level-4.csv", row.names = 1) #row names as sample id
taxadata_5_endov34 <- read.csv("qiime2_output/16S_endo_v34_2021/taxonomy/level-5.csv", row.names = 1) #row names as sample id
taxadata_6_endov34 <- read.csv("qiime2_output/16S_endo_v34_2021/taxonomy/level-6.csv", row.names = 1) #row names as sample id
taxadata_7_endov34 <- read.csv("qiime2_output/16S_endo_v34_2021/taxonomy/level-7.csv", row.names = 1) #row names as sample id

##---- Barplots - monocultures and source samples ----

# define general theme to use for plots
theme_barcharts <- theme(
  axis.text.x = element_text(angle = 90, vjust = 0.5, 
                             hjust = 0.9, size = 11, 
                             color = "grey30"), 
  strip.text.x = element_text(size = 11, color = "gray30", face = "bold"),
  strip.background = element_rect(color = "gray20"),
  panel.grid = element_line(colour = "transparent"),
  panel.grid.major.x = element_line(color = "gray30", linetype = 3),
  panel.border = element_rect(color = "gray30"),
  axis.ticks.x = element_line(colour = "grey30"), # << visibility
  axis.ticks.y = element_line(colour = "grey30"),
  axis.title.x = element_text(size = 11, color = "black"),
  #axis.title.x = element_blank(), #<< visibility
  axis.title.y = element_text(size = 11, color = "black"),
  axis.text.y = element_text(size = 9, color = "grey30"),
  legend.title = element_text(size = 9, color = "black"),
  legend.text = element_text(size = 8, color = "grey30"),
  legend.spacing.y = unit(0.05, 'cm'),
  legend.key.size = unit(0.3, 'cm'),
  plot.tag = element_text(size = 12)
)

# remove samples with <4500 reads
removed_samples_v34 <- c("16SNEGCTRL",
                     "16S0094O",
                     "16S0113C",
                     "16S0118O",
                     "16S0118W",
                     "16S0119C",
                     "16S0092O",
                     "16S0064O")

taxaplot_2_v34 <- 
taxadata_2_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  select(-which(names(.) %in% names(which(apply(., 2, max) <= 0.05)))) %>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.)) %>% #add row with other as 1- sum of abundances
  rename_with(~ gsub("d__Bacteria.p__", " ", .x)) %>%
  t() %>% #transpose
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Phyla (> 5%)", palette = "Paired") +
  scale_y_continuous(breaks = c(0, 0.25, .50, .75, 1), 
                     labels = function(x) x*100, expand = c(0, 0), 
                     limits = c(0, 1)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts

taxaplot_plancto <-
taxadata_5_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.001)))) %>% #select taxa above a percentage / remove taxa below percentage
  #mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Planctomycetota")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Planctomycetota.c__", "", .x)) %>%
  rename_with(~ gsub("028H05.P.BN.P5.o__028H05.P.BN.P5.f__", "", .x)) %>%
  rename_with(~ gsub("OM190.o__OM190.f__", "", .x)) %>%
  rename_with(~ gsub("Phycisphaerae.o__Phycisphaerales.f__Phycisphaeraceae",
                     "order Phycisphaerales;  
                     family *Phycisphaeraceae*", .x)) %>%
  rename_with(~ gsub("Planctomycetes.o__",
                     "order ", .x)) %>%
  rename_with(~ gsub("Gemmatales.f__Gemmataceae",
                     "Gemmatales;  
                     family *Gemmataceae*", .x)) %>%
  rename_with(~ gsub("Pirellulales.f__Pirellulaceae",
                     "Pirellulales;  
                     family *Pirellulaceae*", .x)) %>%
  rename_with(~ gsub("Planctomycetales.f__Rubinisphaeraceae",
                     "Planctomycetales;  
                     family *Rubinisphaeraceae*", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  #filter(SampleSite == "ORAL" | SampleSite == "CLOACA") %>%
  #mutate(SamplePlotName = fct_relevel(SamplePlotName, guni_d_order_labels_ms_dm)) %>%
  #mutate(ExpGen = paste(TypeExp, Genus, sep = "_")) %>%
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Planctomycetota <0.1%", palette = "Paired") +
  # scale_fill_manual(name = "Alphaproteobacteria (> 5%)", values = c("#A6CEE3", "#1F78B4",                                                                                #                                                                    "#B2DF8A", "#33A02C",
  #                                                                   "#FB9A99", "#E31A1C", 
  #                                                                   "#FDBF6F", "#FF7F00",
  #                                                                   "#CAB2D6", "#6A3D9A",
  #                                                                   "#FFFF99", "#B15928",
  #                                                                   "#A6CEE3", "#1F78B4",
  #                                                                   "#B2DF8A", "#E31A1C", 
  #                                                                   "#FDBF6F", "#111111")) +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.04)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  theme(legend.text = element_markdown(),
        legend.position = "bottom")

taxaplot_alpha <-
  taxadata_4_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  #mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Alphaproteobacteria")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria", "", .x)) %>%
  rename_with(~ gsub(".o__", "", .x)) %>%
  # rename_with(~ gsub("Gammaproteobacteria_Incertae_Sedis", "Incertae Sedis", .x)) %>%
  # rename_with(~ gsub(".__", " unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  #scale_fill_brewer(name = "Gammaproteobacteria", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.4)) + #set scale size #<<
  scale_fill_manual(name = "Alphaproteobacteria <1%",
                    values = c("#A6CEE3", "#1F78B4",
                               "#B2DF8A", "#33A02C",
                               "#FB9A99", "#E31A1C",
                               "#FDBF6F", "#FF7F00",
                               "#CAB2D6", "#6A3D9A",
                               "#FFFF99", "#B15928",
                               "#A6CEE3", "#1F78B4",
                               "#B2DF8A", "#E31A1C",
                               "#FDBF6F", "#111111",
                               "grey"
                    )) +
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  theme(legend.text = element_markdown(),
        legend.position = "bottom")

taxaplot_gamma <-
taxadata_4_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.05)))) %>% #select taxa above a percentage / remove taxa below percentage
  #mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Gammaproteobacteria")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria", "", .x)) %>%
  rename_with(~ gsub(".o__", "", .x)) %>%
  rename_with(~ gsub("Gammaproteobacteria_Incertae_Sedis", "Incertae Sedis", .x)) %>%
  rename_with(~ gsub(".__", " unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  #scale_fill_brewer(name = "Gammaproteobacteria", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.8)) + #set scale size #<<
  scale_fill_manual(name = "Gammaproteobacteria <5%",
                    values = c("#A6CEE3", "#1F78B4",
                               "#B2DF8A", "#33A02C",
                               "#FB9A99", "#E31A1C",
                               "#FDBF6F", "#FF7F00",
                               "#CAB2D6", "#6A3D9A",
                               "#FFFF99", "#B15928",
                               "#A6CEE3", "#1F78B4",
                               "#B2DF8A", "#E31A1C",
                               "#FDBF6F", "#111111",
                               "grey"
                    )) +
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  theme(legend.text = element_markdown(),
        legend.position = "bottom")

taxaplot_firmicutes <-
taxadata_4_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  #mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Firmicutes")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Firmicutes.c__", "", .x)) %>%
  rename_with(~ gsub(".o__", "; ", .x)) %>%
  #rename_with(~ gsub(".__", "unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Firmicutes <1%", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.4)) + #set scale size #<<
  # scale_fill_manual(name = "Gammaproteobacteria",
  #                   values = c("#A6CEE3", "#1F78B4",
  #                              "#B2DF8A", "#33A02C",
  #                              "#FB9A99", "#E31A1C",
  #                              "#FDBF6F", "#FF7F00",
  #                              "#CAB2D6", "#6A3D9A",
  #                              "#FFFF99", "#B15928",
  #                              "#A6CEE3", "#1F78B4",
  #                              "#B2DF8A", "#E31A1C",
  #                              "#FDBF6F", "#111111",
  #                              #"#B2DF8A", "#33A02C",
#                              "grey"
#                   )) +
labs(
  x = "Sample ID",
  y = "Relative abundance (%)"
) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  theme(legend.text = element_markdown(),
        legend.position = "bottom")

taxaplot_patescibacteria <-
taxadata_4_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Patescibacteria")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Patescibacteria.c__", "", .x)) %>%
  rename_with(~ gsub(".o__", "; ", .x)) %>%
  #rename_with(~ gsub(".__", "unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Patescibacteria <1%", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.1)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
    ) +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme_barcharts +
    theme(legend.text = element_markdown(),
          legend.position = "bottom")

taxaplot_spirochaetota <-
taxadata_4_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  #select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Spirochaetota")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Spirochaetota.c__", "", .x)) %>%
  rename_with(~ gsub(".o__", "; ", .x)) %>%
  #rename_with(~ gsub(".__", "unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Spirochaetota", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.1)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme_barcharts +
    theme(legend.text = element_markdown(),
          legend.position = "bottom")

taxaplot_verrucomicrobiota <-
taxadata_3_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  #select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Verrucomicrobiota")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Verrucomicrobiota.c__", "", .x)) %>%
  rename_with(~ gsub("d__Bacteria.p__Verrucomicrobiota.__", "Verrucomicrobiota unassigned", .x)) %>%
  # rename_with(~ gsub(".o__", "", .x)) %>%
  #rename_with(~ gsub(".__", "unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Verrucomicrobiota", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.1)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme_barcharts +
    theme(legend.text = element_markdown(),
          legend.position = "bottom")

taxaplot_actinobacteria <-
  taxadata_4_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Actinobacteriota")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Actinobacteriota.c__", "", .x)) %>%
  rename_with(~ gsub(".o__", "; ", .x)) %>%
  rename_with(~ gsub(".__", " unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Actinobacteriota <1%", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.1)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  theme(legend.text = element_markdown(),
        legend.position = "bottom")

taxaplot_bacteroidota <-
  taxadata_4_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Bacteroidota")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Bacteroidota.c__", "", .x)) %>%
  rename_with(~ gsub(".o__", "; ", .x)) %>%
  # rename_with(~ gsub(".__", " unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Bacteroidota", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 1)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  theme(legend.text = element_markdown(),
        legend.position = "bottom")

taxaplot_bdellovibrionota <-
  taxadata_4_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  #select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Bdellovibrionota")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Bdellovibrionota.c__", "", .x)) %>%
  rename_with(~ gsub(".o__", "; ", .x)) %>%
  # rename_with(~ gsub(".__", " unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Bdellovibrionota", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.15)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  theme(legend.text = element_markdown(),
        legend.position = "bottom")

taxaplot_campilobacterota <-
  taxadata_5_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  #select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Campilobacterota")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Campilobacterota.c__Campylobacteria", "", .x)) %>%
  rename_with(~ gsub(".o__", "", .x)) %>%
  rename_with(~ gsub(".f__", "; ", .x)) %>%
  rename_with(~ gsub(".__", " unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Campylobacteria", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.25)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  theme(legend.text = element_markdown(),
        legend.position = "bottom")

taxaplot_cyanobacteria <-
  taxadata_4_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Cyanobacteria")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Cyanobacteria.c__", "", .x)) %>%
  rename_with(~ gsub(".o__", "; ", .x)) %>%
  # rename_with(~ gsub(".f__", "; ", .x)) %>%
  # rename_with(~ gsub(".__", " unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Cyanobacteria <1%", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.1)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  theme(legend.text = element_markdown(),
        legend.position = "bottom")

taxaplot_fusobacteriota <-
  taxadata_5_endov34 %>% #<<
  .[!(row.names(.) %in% removed_samples_v34),1:(ncol(.)-32)] %>% #remove metadata and samples with low reads <4500
  "/"(rowSums(.)) %>% #calculate proportions (rel. abund.) 
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01)))) %>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.)) %>% #select taxa above a percentage / remove taxa below percentage
  select(contains("Fusobacteriota")) %>% #select specific taxa #<<
  rename_with(~ gsub("d__Bacteria.p__Fusobacteriota.c__", "", .x)) %>%
  rename_with(~ gsub(".o__", "; ", .x)) %>%
  rename_with(~ gsub(".f__", "; ", .x)) %>%
  # rename_with(~ gsub(".__", " unassigned", .x)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("16S")), names_to = "SampleID", values_to = "abundance") %>% #<<
  left_join(metadata_endov34) %>% #<<
  ggplot(aes(x = SampleID, y = abundance)) +
  facet_grid(cols = vars(SampleSite), scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 0.8, color = "grey20", linewidth = 0.35) +
  scale_fill_brewer(name = "Fusobacteriota <1%", palette = "Paired") +
  scale_y_continuous(labels = function(x) x*100, 
                     expand = c(0, 0), 
                     limits = c(0, 0.15)) + #set scale size #<<
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)"
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme_barcharts +
  theme(legend.text = element_markdown(),
        legend.position = "bottom")

ggsave("r_output/16S_endo_v34/barchart_phyla.pdf",
       plot = taxaplot_2_v34,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")

ggsave("r_output/16S_endo_v34/barchart_planctomycetota.pdf",
       plot = taxaplot_plancto,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")

ggsave("r_output/16S_endo_v34/barchart_gammaproteobacteria.pdf",
       plot = taxaplot_gamma,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")
ggsave("r_output/16S_endo_v34/barchart_alphaproteobacteria.pdf",
       plot = taxaplot_alpha,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")
ggsave("r_output/16S_endo_v34/barchart_firmicutes.pdf",
       plot = taxaplot_firmicutes,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")

ggsave("r_output/16S_endo_v34/barchart_patescibacteria.pdf",
       plot = taxaplot_patescibacteria,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")

ggsave("r_output/16S_endo_v34/barchart_spirochaetota.pdf",
       plot = taxaplot_spirochaetota,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")

ggsave("r_output/16S_endo_v34/barchart_verrucomicrobiota.pdf",
       plot = taxaplot_verrucomicrobiota,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")

ggsave("r_output/16S_endo_v34/barchart_actinobacteriota.pdf",
       plot = taxaplot_actinobacteria,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")

ggsave("r_output/16S_endo_v34/barchart_bacteroidota.pdf",
       plot = taxaplot_bacteroidota,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")

ggsave("r_output/16S_endo_v34/barchart_bdellovibrionota.pdf",
       plot = taxaplot_bdellovibrionota,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")

ggsave("r_output/16S_endo_v34/barchart_campilobacterota.pdf",
       plot = taxaplot_campilobacterota,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")

ggsave("r_output/16S_endo_v34/barchart_cyano.pdf",
       plot = taxaplot_cyanobacteria,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")

ggsave("r_output/16S_endo_v34/barchart_fusobacteriota.pdf",
       plot = taxaplot_fusobacteriota,
       device = cairo_pdf,
       height = 100,
       width = 220,
       units = "mm")



##----- Differential abundance tests -----
# ANCOM-BC2
# read qiime to phyloseq object
endo_phy_data <- qza_to_phyloseq(features = "qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qza",
                                 tree = "qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-rooted-tree.qza",
                                 taxonomy = "qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza",
                                 metadata = "master_metadata/16SendoV34.tsv")
endo_phy_data_genus <- qza_to_phyloseq(features = "qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-6.qza",
                                       metadata = "master_metadata/16SendoV34.tsv")
endo_phy_data_species <- qza_to_phyloseq(features = "qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-7.qza",
                                         metadata = "master_metadata/16SendoV34.tsv")

# import asv numbers for rep seqs
asvs <- readxl::read_xlsx("qiime2_output/16S_endo_v34_2021/merged-no_filter/metadata-asvs.xlsx")
asvs <- asvs %>% select("taxon", "ASVno", "Taxon")

# set seed
set.seed(123)

# ancombc2 on all samples grouped by samplesite (ASV and Species and Genus)
output_all_asv <- ancombc2(data = endo_phy_data,
                           fix_formula = "SampleSite",
                           group = "SampleSite",
                           global = TRUE,
                           pairwise = TRUE,
                           n_cl = 2)
# on genus
output_genusnew <- ancombc2(data = endo_phy_data_genus,
                            fix_formula = "SampleSite",
                            group = "SampleSite",
                            global = TRUE,
                            pairwise = TRUE,
                            n_cl = 2,
                            verbose = TRUE)
# on species
output_speciesnew <- ancombc2(data = endo_phy_data_species,
                              fix_formula = "SampleSite",
                              group = "SampleSite",
                              global = TRUE,
                              pairwise = TRUE,
                              n_cl = 2)

write_tsv(output_all_asv$res %>%
            left_join(output_all_asv$pseudo_sens_tab) %>%
            left_join(asvs), 
          file = "r_output/16S_endo_v34/ancombc2_results_asvs.tsv")
write_tsv(output_all_asv$res_pair %>%
            left_join(output_all_asv$pseudo_sens_tab) %>%
            left_join(asvs), 
          file = "r_output/16S_endo_v34/ancombc2_results_pairwise_asvs.tsv")
write_tsv(output_all_asv$res_global %>%
            left_join(output_all_asv$pseudo_sens_tab) %>%
            left_join(asvs), 
          file = "r_output/16S_endo_v34/ancombc2_results_global_asvs.tsv")

write_tsv(output_speciesnew$res %>%
            left_join(output_speciesnew$pseudo_sens_tab), 
          file = "r_output/16S_endo_v34/ancombc2_results_spec.tsv")
write_tsv(output_speciesnew$res_pair %>%
            left_join(output_speciesnew$pseudo_sens_tab), 
          file = "r_output/16S_endo_v34/ancombc2_results_pairwise_spec.tsv")
write_tsv(output_speciesnew$res_global %>%
            left_join(output_speciesnew$pseudo_sens_tab), 
          file = "r_output/16S_endo_v34/ancombc2_results_global_spec.tsv")

write_tsv(output_genusnew$res %>%
            left_join(output_genusnew$pseudo_sens_tab), 
          file = "r_output/16S_endo_v34/ancombc2_results_gen.tsv")
write_tsv(output_genusnew$res_pair %>%
            left_join(output_genusnew$pseudo_sens_tab), 
          file = "r_output/16S_endo_v34/ancombc2_results_pairwise_gen.tsv")
write_tsv(output_genusnew$res_global %>%
            left_join(output_genusnew$pseudo_sens_tab), 
          file = "r_output/16S_endo_v34/ancombc2_results_global_gen.tsv")

#----- Session info -----
sessioninfo_da <- sessionInfo()
capture.output(sessioninfo, file = "r_output/16S_endo_v34/sessioninfo_da.txt")

