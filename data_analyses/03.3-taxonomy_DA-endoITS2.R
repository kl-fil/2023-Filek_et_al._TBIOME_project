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
library(patchwork) #collect multiple ggplots together v1.1.2

#colors
c("#93624d", "#89d6dc")
c("#462f25", "#16474b")

#set global theme for ggplot
theme_set(theme_light(base_size = 11, base_family = "Arial")) #<< font!

# "#<<" comment indicates the line of code that needs to be changed when using 
# a different dataset for the same code

##---- make directory for visualizations ----
dir.create("r_output/ITS2_endo/", recursive = TRUE)

##---- Load metadata and data-----
# load metadata in tsv format
metadata_ITS <- readr::read_tsv("master_metadata/16SendoV34.tsv") %>% 
  mutate(SampleID = gsub("16S", "ITS", SampleID)) %>%
  filter(!SampleSite %in% "ORAL") %>%
  filter(!SampleID %in% "ITS0146C")

##---- Taxa bar plots ----

# define general theme to use for plots
theme_barcharts <- theme(
  strip.text.x = element_text(size = 8, color = "black", face = "bold"),
  strip.background = element_rect(color = "black", fill = "gray90"),
  panel.grid = element_line(colour = "transparent"),
  panel.grid.major.x = element_line(color = "gray30", linetype = 3),
  panel.border = element_rect(color = "black"),
  axis.ticks.x = element_line(colour = "black"), # << visibility
  axis.ticks.y = element_line(colour = "black"),
  axis.title.x = element_text(size = 10, color = "black"),
  axis.text.x = element_text(angle = 60, vjust = 1, 
                             hjust = 1, size = 8, 
                             color = "black"), 
  axis.title.y = element_text(size = 10, color = "black"),
  axis.text.y = element_text(size = 8, color = "black"),
  legend.title = element_text(size = 8, color = "black"),
  legend.text = element_text(size = 8, color = "black"),
  legend.spacing = unit(0.1, 'cm'),
  legend.key.size = unit(0.3, 'cm'),
  plot.tag = element_text(size = 10)
)

taxadata_2_ITS <- read.csv("qiime2_output/ITS2_analysis-outputs/taxonomy_level_2.csv", row.names = 1)
taxadata_2_ITS <- taxadata_2_ITS[1:29,1:16] %>% #izvlacenje prvih 16 stupaca koji sadrze samo podatke o abundanciji taksona i prvih 29 redaka da izbacimo negativnu kontrolu
  mutate(" Unidentified Fungi" = k__Fungi.__ + k__Fungi.p__unidentified) %>%
  select(- c(k__Fungi.__,k__Fungi.p__unidentified))           #spajanje neidentificiranih fungi u jedan stupac i brisanje prethodnih stupaca u kojima su bili odvojeni

phyla_barplot_ITS <-
  taxadata_2_ITS %>%
  "/"(rowSums(.)) %>%
  select(-which(names(.) %in% names(which(apply(., 2, max) < 0.01))))%>% #select taxa above a percentage / remove taxa below percentage
  mutate(Other = 1 - rowSums(.))%>%   #other column - taxa below 1% abundance
  rename_with(~ gsub("k__Fungi.p__", " ", .x))%>%
  t() %>%# #transpose
  as.data.frame() %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = contains(c("ITS")), names_to = "SampleID", values_to = "abundance") %>% 
  left_join(metadata_ITS) %>%
  ggplot(aes(x = SampleID, y = abundance)) + 
  geom_bar(aes(fill = Taxa), position = "stack", stat = "identity", width = 1, color = "black", size = 0.2) +
  scale_fill_brewer(name = "Phyla (> 1%)", palette = "Set3") +
  scale_y_continuous(breaks = c(0, 0.25, .50, .75, 1), 
                     labels = function(x) x*100, expand = c(0, 0), 
                     limits = c(0, 1)) +
  theme_barcharts +
  facet_grid(cols = vars(SampleSite),scales = "free", space = "free") + 
  labs(
    x = "Sample ID",
    y = "Relative abundance (%)", #<<
    tag = "d)" #<<
  )
  
faith_alpha_rain_ITS + labs(tag = "a)")

fig4_pub <- (faith_alpha_rain_ITS + labs(tag = "a)")) + (shannon_alpha_rain_ITS + labs(tag = "b)")) + phyla_barplot_ITS + 
  plot_layout(design = "AB###
                        CCCCC")

ggsave("r_output/ITS2_endo/fig4_pub.pdf",
       plot = fig4_pub,
       device=cairo_pdf,
       height = 120,
       width = 175,
       units = "mm")

##----- Differential abundance tests -----
# ANCOM-BC2
phytable_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/ANCOMBC2/feature-frequency-filtered-table.qza")
phytable_spec_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/ANCOMBC2/table-lvl-7.qza")
phytable_gen_ITS <- read_qza("qiime2_output/ITS2_analysis-outputs/ANCOMBC2/table-lvl-6.qza")

colnames(phytable$data)
asvtable <- phytable$data
class(asvtable)

OTU_ITS = otu_table(phytable_ITS$data, taxa_are_rows = TRUE)
asv_phy_ITS <- phyloseq(OTU_ITS)
OTU_spec_ITS <- otu_table(phytable_spec_ITS$data, taxa_are_rows = TRUE)
spec_phy_ITS <- phyloseq(OTU_spec_ITS)
OTU_gen_ITS <- otu_table(phytable_gen_ITS$data, taxa_are_rows = TRUE)
gen_phy_ITS <- phyloseq(OTU_gen_ITS)


df_sampleinfo_ITS <- as.data.frame(metadata_ITS)
rownames(df_sampleinfo_ITS) <- df_sampleinfo_ITS$SampleID
SAMPDATA_ITS <- sample_data(df_sampleinfo_ITS)
endo_phy_data_ITS <- merge_phyloseq(asv_phy_ITS, SAMPDATA_ITS)
endo_phy_data_spec_ITS <- merge_phyloseq(spec_phy_ITS, SAMPDATA_ITS)
endo_phy_data_gen_ITS <- merge_phyloseq(gen_phy_ITS, SAMPDATA_ITS)

#import asv numbers for representative sequences
asvs_ITS <- readxl::read_xlsx("qiime2_output/ITS2_analysis-outputs/merged-table-seq/metadata_asvno.xlsx")
asvs_ITS <- asvs_ITS %>% select("taxon", "ASVno", "Taxon")

output_all_asv_ITS <- ancombc2(data = endo_phy_data_ITS,
                               fix_formula = "SampleSite",
                               n_cl = 2)

output_all_spec_ITS <- ancombc2(data = endo_phy_data_spec_ITS,
                                fix_formula = "SampleSite",
                                n_cl = 2)

output_all_gen_ITS <- ancombc2(data = endo_phy_data_gen_ITS,
                               fix_formula = "SampleSite",
                               n_cl = 2)

write_tsv(output_all_asv_ITS$res %>%
            left_join(output_all_asv_ITS$pseudo_sens_tab) %>%
            left_join(asvs_ITS), 
          file = "r_output/ITS2_endo/ancombc2_results_asvs.tsv")

write_tsv(output_all_spec_ITS$res %>%
            left_join(output_all_spec_ITS$pseudo_sens_tab), 
          file = "r_output/ITS2_endo/ancombc2_results_spec.tsv")

write_tsv(output_all_gen_ITS$res %>%
            left_join(output_all_gen_ITS$pseudo_sens_tab), 
          file = "r_output/ITS2_endo/ancombc2_results_gen.tsv")