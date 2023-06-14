## Dataviz script for TurtleBiome Holobiont 2023 manuscript
## Processing and visualizing data related to sequencing results for
## endozoic and epizoic samples (16S rRNA gene, trimmed to V4 region)
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
library(pheatmap)
library(compositions)
library(RColorBrewer)

#set global theme for ggplot
theme_set(theme_light(base_size = 11, base_family = "Arial")) #<< font!

# "#<<" comment indicates the line of code that needs to be changed when using 
# a different dataset for the same code

##---- make directory for visualizations ----
dir.create("r_output/16S_epi_v4_merged/", recursive = TRUE)

##---- Load metadata and data-----
# load metadata in tsv format
metadata_epiv4 <- readr::read_tsv("master_metadata/16S_all_metadata.tsv")

#colors
c("#93624d", "#ff9999", "#3BB848", "#89d6dc")
c("#462f25","#750000", "#143e18", "#16474b")

##----- Differential abundance tests -----
# ANCOM-BC2
#read qiime to phyloseq object
phytable <- read_qza("qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza")
phytable_spec <- read_qza("qiime2_output/16S_v4_merged/taxonomy_collapsed/table-7.qza")
phytable_gen <- read_qza("qiime2_output/16S_v4_merged/taxonomy_collapsed/table-6.qza")
colnames(phytable$data)
asvtable <- phytable$data
class(asvtable)

OTU = otu_table(asvtable, taxa_are_rows = TRUE)
testphy <- phyloseq(OTU)
OTU_spec <- otu_table(phytable_spec$data, taxa_are_rows = TRUE)
spec_phy <- phyloseq(OTU_spec)
OTU_gen <- otu_table(phytable_gen$data, taxa_are_rows = TRUE)
gen_phy <- phyloseq(OTU_gen)


df_sampleinfo <- as.data.frame(metadata_epiv4)
rownames(df_sampleinfo) <- df_sampleinfo$SampleID
SAMPDATA <- sample_data(df_sampleinfo)
endo_phy_data_v4_2 <- merge_phyloseq(testphy, SAMPDATA)
epi_phy_data_spec_v4 <- merge_phyloseq(spec_phy, SAMPDATA)
epi_phy_data_gen_v4 <- merge_phyloseq(gen_phy, SAMPDATA)

#import asv numbers for rep seqs
asvs_v4 <- readxl::read_xlsx("qiime2_output/16S_v4_merged/merged-no_filter/metadata_asvs.xlsx")
asvs_v4 <- asvs_v4 %>% select("taxon", "ASVno", "Taxon")

#set seed
set.seed(123)

output_all_asv_v4_zero <- ancombc2(data = endo_phy_data_v4_2,
                                   fix_formula = "SampleSite",
                                   group = "SampleSite",
                                   global = TRUE,
                                   pairwise = TRUE,
                                   n_cl = 2,
                                   struc_zero = TRUE)

output_all_spec_v4_zero <- ancombc2(data = epi_phy_data_spec_v4,
                                    fix_formula = "SampleSite",
                                    group = "SampleSite",
                                    global = TRUE,
                                    pairwise = TRUE,
                                    n_cl = 2,
                                    struc_zero = TRUE)

output_all_gen_v4_zero <- ancombc2(data = epi_phy_data_gen_v4,
                                   fix_formula = "SampleSite",
                                   group = "SampleSite",
                                   global = TRUE,
                                   pairwise = TRUE,
                                   n_cl = 2,
                                   struc_zero = TRUE)

asvs_v4 <- readxl::read_xlsx("qiime2_output/16S_v4_merged/merged-no_filter/metadata_asvs.xlsx")
asvs_v4 <- asvs_v4 %>% select("taxon", "ASVno", "Taxon")

output_all_asv_v4
output_all_asv_v4_zero
output_all_spec_v4
output_all_spec_v4_zero
output_all_gen_v4
output_all_gen_v4_zero

#save results with pseudocount sensitivity and struct. zeros
write_tsv(output_all_asv_v4_zero$res %>%
            left_join(output_all_asv_v4_zero$pseudo_sens_tab) %>%
            left_join(output_all_asv_v4_zero$zero_ind) %>%
            left_join(asvs_v4), 
          file = "r_output/16S_epi_v4_merged/ancombc2-zero_results_asvs.tsv")
write_tsv(output_all_asv_v4_zero$res_pair %>%
            left_join(output_all_asv_v4_zero$pseudo_sens_tab) %>%
            left_join(output_all_asv_v4_zero$zero_ind) %>%
            left_join(asvs_v4), 
          file = "r_output/16S_epi_v4_merged/ancombc2-zero_results_pairwise_asvs.tsv")
write_tsv(output_all_asv_v4_zero$res_global %>%
            left_join(output_all_asv_v4_zero$pseudo_sens_tab) %>%
            left_join(output_all_asv_v4_zero$zero_ind) %>%
            left_join(asvs_v4), 
          file = "r_output/16S_epi_v4_merged/ancombc2-zero_results_global_asvs.tsv")

write_tsv(output_all_spec_v4_zero$res %>%
            left_join(output_all_spec_v4_zero$pseudo_sens_tab) %>%
            left_join(output_all_spec_v4_zero$zero_ind), 
          file = "r_output/16S_epi_v4_merged/ancombc2-zero_results_spec.tsv")
write_tsv(output_all_spec_v4_zero$res_pair %>%
            left_join(output_all_spec_v4_zero$pseudo_sens_tab) %>%
            left_join(output_all_spec_v4_zero$zero_ind),  
          file = "r_output/16S_epi_v4_merged/ancombc2-zero_results_pairwise_spec.tsv")
write_tsv(output_all_spec_v4_zero$res_global %>%
            left_join(output_all_spec_v4_zero$pseudo_sens_tab) %>%
            left_join(output_all_spec_v4_zero$zero_ind),  
          file = "r_output/16S_epi_v4_merged/ancombc2-zero_results_global_spec.tsv")

write_tsv(output_all_gen_v4_zero$res %>%
            left_join(output_all_gen_v4_zero$pseudo_sens_tab) %>%
            left_join(output_all_gen_v4_zero$zero_ind), 
          file = "r_output/16S_epi_v4_merged/ancombc2-zero_results_gen.tsv")
write_tsv(output_all_gen_v4_zero$res_pair %>%
            left_join(output_all_gen_v4_zero$pseudo_sens_tab) %>%
            left_join(output_all_gen_v4_zero$zero_ind), 
          file = "r_output/16S_epi_v4_merged/ancombc2-zero_results_pairwise_gen.tsv")
write_tsv(output_all_gen_v4_zero$res_global %>%
            left_join(output_all_gen_v4_zero$pseudo_sens_tab) %>%
            left_join(output_all_gen_v4_zero$zero_ind),
          file = "r_output/16S_epi_v4_merged/ancombc2-zero_results_global_gen.tsv")

##----- Heatmaps -----
# generate csv files with abundances from mono-source/taxonomy/taxa-bar-plots.qzv 
# save each level as csv in the mono-source/taxonomy folder
# import csv files taxa as columns, samples as rows
taxadata_7_epiv4 <- read.csv("qiime2_output/16S_v4_merged/taxa_bar_plot_export/level-7.csv", row.names = 1) #row names as sample id

# remove samples with low reads
removed_samples_v4 <- c("16SNEGCTRL",
                        "16S0094O",
                        "16S0113C",
                        "16S0118O",
                        "16S0118W",
                        "16S0119C",
                        "16S0092O",
                        "16S0064O")

# generate relative abundance heatmap of differentially abundant taxa collapsed at species level
# relative abundance based on clr transformed counts

heatmap_df  <- as.data.frame(metadata_epiv4)
row.names(heatmap_df) <- heatmap_df$SampleID

# set information for columns (samples)
col_heat <- heatmap_df %>%
  select(c("SampleSite", ))

# remove low read samples
clr_taxa7_no_t <- clr(taxadata_7_epiv4 %>% #<<
                        .[!(row.names(.) %in% removed_samples_v4),1:(ncol(.)-32)])

clr_taxa7_no_t <- as.data.frame(t(clr_taxa7_no_t))

# filter clr transformed taxa counts to features detected as differentially abundant across sample sites by ANCOM-BC2
filtered_clr_taxa7_no_t <- clr_taxa7_no_t %>%
  filter(row.names(clr_taxa7_no_t) %in% c("d__Bacteria.__.__.__.__.__.__",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.__.__.__",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides.__",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Marinifilaceae.g__Marinifilum.s__uncultured_bacterium",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Rikenellaceae.g__Rikenellaceae_RC9_gut_group.s__uncultured_bacterium",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Cyclobacteriaceae.g__Ekhidna.s__Ekhidna_sp.",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Cyclobacteriaceae.g__uncultured.__",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Cyclobacteriaceae.g__uncultured.s__uncultured_bacterium",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Cryomorphaceae.g__uncultured.__",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.g__Algitalea.s__uncultured_bacterium",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.g__Maritimimonas.s__uncultured_bacterium",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.g__NS3a_marine_group.__",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.g__Pseudofulvibacter.__",
                                          "d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.g__Tenacibaculum.s__uncultured_bacterium",
                                          "d__Bacteria.p__Bacteroidota.c__Rhodothermia.o__Balneolales.f__Balneolaceae.g__Balneola.__",
                                          "d__Bacteria.p__Campilobacterota.c__Campylobacteria.o__Campylobacterales.f__Campylobacteraceae.g__Campylobacter.s__uncultured_bacterium",
                                          "d__Bacteria.p__Cyanobacteria.c__Cyanobacteriia.o__Phormidesmiales.f__Phormidesmiaceae.g__Acrophormium_PCC.7375.s__Leptolyngbya_sp.",
                                          "d__Bacteria.p__Deinococcota.c__Deinococci.o__Deinococcales.f__Trueperaceae.g__Truepera.s__uncultured_organism",
                                          "d__Bacteria.p__Dependentiae.c__Babeliae.o__Babeliales.__.__.__",
                                          "d__Bacteria.p__Dependentiae.c__Babeliae.o__Babeliales.f__Babeliales.g__Babeliales.s__uncultured_bacterium",
                                          "d__Bacteria.p__Patescibacteria.c__Gracilibacteria.o__Absconditabacteriales_.SR1..f__Absconditabacteriales_.SR1..g__Absconditabacteriales_.SR1..s__uncultured_bacterium",
                                          "d__Bacteria.p__Patescibacteria.c__Saccharimonadia.o__Saccharimonadales.f__Saccharimonadales.g__Saccharimonadales.s__uncultured_bacterium",
                                          "d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Pirellulales.f__Pirellulaceae.g__Blastopirellula.__",
                                          "d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Pirellulales.f__Pirellulaceae.g__Pir4_lineage.s__uncultured_bacterium",
                                          "d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Pirellulales.f__Pirellulaceae.g__Rhodopirellula.__",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.__.__.__.__",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Hyphomonadaceae.__.__",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Hyphomonadaceae.g__Litorimonas.__",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Kordiimonadales.f__Kordiimonadaceae.g__Kordiimonas.__",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Micavibrionales.f__Micavibrionaceae.g__uncultured.s__uncultured_organism",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.__.__",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.g__Ahrensia.s__uncultured_bacterium",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Stappiaceae.g__Labrenzia.__",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhodospirillales.f__Thalassospiraceae.g__Thalassospira.__",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.__.__",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.g__Altererythrobacter.__",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.g__Erythrobacter.__",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.g__Erythrobacter.s__uncultured_bacterium",
                                          "d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.g__Sphingomonas.__",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae.__.__",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae.g__Alteromonas.__",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae.g__Glaciecola.s__uncultured_bacterium",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Colwelliaceae.g__Thalassotalea.__",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Pseudoalteromonadaceae.g__Pseudoalteromonas.__",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Shewanellaceae.g__Shewanella.__",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Cardiobacteriales.f__Cardiobacteriaceae.g__Cardiobacterium.__",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Cardiobacteriales.f__Cardiobacteriaceae.g__uncultured.s__uncultured_bacterium",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Cellvibrionales.f__Cellvibrionaceae.g__Aestuariicella.s__uncultured_bacterium",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Cellvibrionales.f__Halieaceae.__.__",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Cellvibrionales.f__Halieaceae.g__Halioglobus.s__uncultured_bacterium",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Cellvibrionales.f__Spongiibacteraceae.__.__",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Coxiellales.f__Coxiellaceae.g__Coxiella.s__uncultured_bacterium",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Oceanospirillales.f__Marinomonadaceae.g__Marinomonas.s__uncultured_marine",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Oceanospirillales.f__Nitrincolaceae.g__uncultured.s__uncultured_bacterium",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Oceanospirillales.f__Oleiphilaceae.g__Oleiphilus.__",
                                          "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Moraxellaceae.g__uncultured.s__uncultured_bacterium",
                                          "d__Bacteria.p__Verrucomicrobiota.c__Kiritimatiellae.o__WCHB1.41.f__WCHB1.41.g__WCHB1.41.s__uncultured_bacterium",
                                          "d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Opitutales.f__Puniceicoccaceae.g__Lentimonas.s__uncultured_bacterium"))

# find row names to rename and make more readable
row_names_clr_taxa7_no_t <- row.names(filtered_clr_taxa7_no_t)
row.names(filtered_clr_taxa7_no_t) <- c("Patescibacteria; Saccharimonadales; uncultured Saccharimonadales",
                                        "Proteobacteria; Alphaproteobacteria",
                                        "Verrucomicrobiota; Puniceicoccaceae; uncultured Lentimonas",
                                        "Bacteroidota; Flavobacteriaceae; uncultured Tenacibaculum",
                                        "Planctomycetota; uncultured Pirellulaceae (Pir4 lineage)",
                                        "Proteobacteria; Hyphomonadaceae",
                                        "Proteobacteria; Oleiphilaceae; Oleiphilus",
                                        "Planctomycetota; Pirellulaceae; Blastopirellula",
                                        "Proteobacteria; Coxiellaceae; uncultured Coxiella",
                                        "Proteobacteria; Halieaceae",
                                        "Proteobacteria; Stappiaceae; Labrenzia",
                                        "Proteobacteria; uncultured Cardiobacteriaceae",
                                        "unassigned Bacteria",
                                        "Verrucomicrobiota; uncultured WCHB1.41",
                                        "Proteobacteria; Sphingomonadaceae",
                                        "Proteobacteria; Alteromonadaceae",
                                        "Bacteroidota; Balneolaceae; Balneola",
                                        "Bacteroidota; uncultured Cyclobacteriaceae 1",
                                        "Proteobacteria; uncultured Micavibrionaceae",
                                        "Proteobacteria; Alteromonadaceae; Alteromonas",
                                        "Bacteroidota; Rikenellaceae; uncultured Rikenellaceae RC9 gut group",
                                        "Cyanobacteria; Phormidesmiaceae; Leptolyngbya sp.",
                                        "Proteobacteria; Cellvibrionaceae; uncultured Aestuariicella",
                                        "Proteobacteria; Spongiibacteraceae",
                                        "Dependentiae; Babeliales; uncultured Babeliales",
                                        "Proteobacteria; Rhizobiaceae",
                                        "Bacteroidota; uncultured Cyclobacteriaceae 2",
                                        "Bacteroidota; uncultured Cryomorphaceae",
                                        "Proteobacteria; Sphingomonadaceae; Sphingomonas",
                                        "Proteobacteria; Pseudoalteromonadaceae; Pseudoalteromonas",
                                        "Bacteroidota; Flavobacteriaceae; NS3a marine group",
                                        "Bacteroidota; Flavobacteriaceae; uncultured Algitalea",
                                        "Patescibacteria; uncultured Absconditabacteriales (SR1)",
                                        "Proteobacteria; Sphingomonadaceae; Erythrobacter",
                                        "Bacteroidota; Bacteroidaceae; Bacteroides",
                                        "Proteobacteria; uncultured Nitrincolaceae",
                                        "Proteobacteria; Alteromonadaceae; uncultured Glaciecola",
                                        "Proteobacteria; Colwelliaceae; Thalassotalea",
                                        "Proteobacteria; Shewanellaceae; Shewanella",
                                        "Campilobacterota; Campylobacteraceae; uncultured Campylobacter",
                                        "Bacteroidota; Flavobacteriaceae; Pseudofulvibacter",
                                        "Dependentiae; Babeliales",
                                        "Deinococcota; Trueperaceae; uncultured Truepera",
                                        "Proteobacteria; Sphingomonadaceae; Altererythrobacter",
                                        "Bacteroidota; Flavobacteriaceae; uncultured Maritimimonas",
                                        "Planctomycetota; Pirellulaceae; Rhodopirellula",
                                        "Bacteroidota; Marinifilaceae; uncultured Marinifilum",
                                        "Proteobacteria; Kordiimonadaceae; Kordiimonas",
                                        "Bacteroidota; Cyclobacteriaceae; Ekhidna sp.",
                                        "Bacteroidota; Bacteroidales",
                                        "Proteobacteria; Rhizobiaceae; uncultured Ahrensia",
                                        "Proteobacteria; Thalassospiraceae; Thalassospira",
                                        "Proteobacteria; uncultured Moraxellaceae",
                                        "Proteobacteria; Halieaceae; uncultured Halioglobus",
                                        "Proteobacteria; Hyphomonadaceae; Litorimonas",
                                        "Proteobacteria; Sphingomonadaceae; uncultured Erythrobacter",
                                        "Proteobacteria; Marinomonadaceae; uncultured Marinomonas",
                                        "Proteobacteria; Cardiobacteriaceae; Cardiobacterium")

# determine maximum clr value for heatmap legend
rg <- max(abs(filtered_clr_taxa7_no_t))

# produce heatmap
clr7_ancombc_heatmap <- 
  pheatmap(mat = filtered_clr_taxa7_no_t,
           show_colnames = TRUE,
           show_rownames = TRUE,
           cluster_rows = TRUE,
           cutree_cols = 7,
           cutree_rows = 12,
           #cutree_rows = 10,
           clustering_method = "complete",
           annotation_col = heatmap_df %>%
             select(c("SampleSite")),
           color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "PRGn")))(100),
           breaks = seq(-rg, rg, length.out = 100),
           fontsize = 8,
           cellwidth = 8,
           cellheight = 8,
           annotation_colors = list(
             SampleSite = c(CARAPACE = "#3BB848",
                            CLOACA = "#93624d",
                            ORAL = "#ff9999",
                            "TANK WATER" = "#89d6dc")))

ggsave("r_output/16S_epi_v4_merged/clr_ancom-lvl7-heatmap.pdf",
       plot = clr7_ancombc_heatmap,
       device = cairo_pdf,
       height = 220,
       width = 320,
       units = "mm")
 