# -----------------------------------------------------------------------------# 
# QIIME2 16S NGS data analyses for TurtleBIOME project (Bosak lab, University
# of Zagreb)
# Paper: TBD holobiont
# By: Klara Filek
# Samples: cloacal, oral, and tank water samples of loggerhead sea turtles
# Objective: analyze endozoic V34 dataset independently
# Region: V3-V4 region of 16S rRNA gene
# Primers: 314F and 805R (Klindworth et al. 2013)
# Platform: Illumina MiSeq v3 (300x2 bp paired-end)
# Sequences available at European Nucleotide Archive under accession:
# - PRJEB62752
# ZENODO data DOI: 10.5281/zenodo.8054926
# Analyses performed in Qiime2 v. 2023.2
# ---System versions---
# Python version: 3.8.16
# QIIME 2 release: 2023.2
# QIIME 2 version: 2023.2.0
# q2cli version: 2023.2.0
# ---Installed plugins---
# alignment: 2023.2.0
# composition: 2023.2.0
# cutadapt: 2023.2.0
# dada2: 2023.2.0
# deblur: 2023.2.0
# deicode: 0.2.4
# demux: 2023.2.0
# diversity: 2023.2.0
# diversity-lib: 2023.2.0
# emperor: 2023.2.0
# empress: 1.2.0
# feature-classifier: 2023.2.0
# feature-table: 2023.2.0
# fragment-insertion: 2023.2.0
# gemelli: 0.0.8
# gneiss: 2023.2.0
# longitudinal: 2023.2.0
# metadata: 2023.2.0
# phylogeny: 2023.2.0
# quality-control: 2023.2.0
# quality-filter: 2023.2.0
# qurro: 0.8.0
# rescript: 2023.2.0
# sample-classifier: 2023.2.0
# taxa: 2023.2.0
# types: 2023.2.0
# vsearch: 2023.2.0
#-----------------------------------------------------------------------------#

# Note: prior to assigning taxonomy a classifier needs to be downloaded and put in classifier folder
# v34 SILVA classifier: made by Lucija Kanjer following Qiime 2 instructions - see Mendeley Data

#----- QIIME 2 environment activation in Conda
conda activate qiime2-2023.2

#----- Change to working directory and set up folders

cd data_analyses
mkdir -p input_data/16S_endo_v34_2021
# copy fastq.gz files from source (ENA) to input_data/16S_endo_v34_2021
# the sequences were processed by the sequencing company and non-biological sequences were removed 

# create directory for qiime2 output
mkdir -p qiime2_output/16S_endo_v34_2021/{demux-dada2,filtered,taxonomy,phylogeny}

#----- Import 16S V34 cloacal, oral, tank water data in Cassava 1.8 format

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input_data/16S_endo_v34_2021 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-demux-paired-end.qza
qiime demux summarize \
  --i-data qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-demux-paired-end.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-demux-paired-end.qzv

#----- Remove primers in sequences using cutadapt prior to DADA2, avoiding readthrough due to sequence length variability
#for V34
#  --p-adapter-f GGATTAGATACCCBDGTAGTC\ reverse primer reverse complement
#  --p-front-f CCTACGGGNGGCWGCAG\ forward primer 314F
#  --p-adapter-r CTGCWGCCNCCCGTAGG\ forward primer reverse complement
#  --p-front-r GACTACHVGGGTATCTAATCC\ reverse primer 805R

qiime cutadapt trim-paired \
  --p-cores 4 \
  --i-demultiplexed-sequences qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-demux-paired-end.qza \
  --p-adapter-f GGATTAGATACCCBDGTAGTC\
  --p-front-f CCTACGGGNGGCWGCAG\
  --p-adapter-r CTGCWGCCNCCCGTAGG\
  --p-front-r GACTACHVGGGTATCTAATCC\
  --o-trimmed-sequences qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-demux-trimmed.qza \
  --verbose

qiime demux summarize \
  --i-data qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-demux-trimmed.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-demux-trimmed.qzv

#----- Denoising sequences with DADA2

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-demux-trimmed.qza \
  --p-trunc-len-f 277 \
  --p-trunc-len-r 235 \
  --o-representative-sequences qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-rep-seqs-dada2-trimmed.qza \
  --o-table qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-table-dada2-trimmed.qza \
  --o-denoising-stats qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-stats-dada2-trimmed.qza \
  --p-n-threads 0

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-stats-dada2-trimmed.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-stats-dada2-trimmed.qzv

qiime feature-table tabulate-seqs \
  --i-data qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-rep-seqs-dada2-trimmed.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-rep-seqs-dada2-trimmed.qzv

qiime feature-table summarize \
  --i-table qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-table-dada2-trimmed.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-table-dada2-trimmed.qzv \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv

#----- Assigning taxonomy via SILVA v.138 
# classifier was made by Lucija Kanjer following the instructions for qiime2 and using primers 314F and 805R (q2-feature-classifier):
# qiime feature-classifier extract-reads \
# --i-sequences silva-138-99-seqs.qza \
# --p-f-primer CCTACGGGNGGCWGCAG \
# --p-r-primer GACTACHVGGGTATCTAATCC \
# --p-min-length 250 \
# --p-max-length 500 \
# --o-reads ref-seqs-V34.qza
# qiime feature-classifier fit-classifier-naive-bayes \
# --i-reference-reads ref-seqs-V34.qza \
# --i-reference-taxonomy silva-138-99-tax.qza \
# --o-classifier classifier-V34.qza

qiime feature-classifier classify-sklearn \
  --i-classifier classifiers/classifier-V34-LK.qza \
  --p-reads-per-batch 10000 \
  --p-n-jobs -2 \
  --i-reads qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-rep-seqs-dada2-trimmed.qza \
  --o-classification qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza &&

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qzv &&

qiime taxa barplot \
  --i-table qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-table-dada2-trimmed.qza \
  --i-taxonomy qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --o-visualization qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxa-bar-plots.qzv

#----- Align sequences using mafft and fasttree

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-rep-seqs-dada2-trimmed.qza \
  --p-n-threads auto \
  --o-alignment qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-aligned-rep-seqs.qza \
  --o-masked-alignment qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-masked-aligned-rep-seqs.qza \
  --o-tree qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-unrooted-tree.qza \
  --o-rooted-tree qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-rooted-tree.qza

qiime empress tree-plot \
    --i-tree qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-rooted-tree.qza \
    --m-feature-metadata-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
    --o-visualization qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-rooted-tree-viz-empress.qzv

#----- Filter based on taxonomy (filter mitochondria and chloroplast)

qiime taxa filter-table \
  --i-table qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-table-dada2-trimmed.qza \
  --i-taxonomy qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza

qiime feature-table summarize \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qzv \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv

# New taxa bar plot without mitochondria and chloroplast sequences

qiime taxa barplot \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --i-taxonomy qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-taxa-bar-plots-filtered.qzv

#----- Alpha rarefaction to determine approapriate depth for downstream analyses

mkdir qiime2_output/16S_endo_v34_2021/filtered/diversity

qiime diversity alpha-rarefaction \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --i-phylogeny qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-rooted-tree.qza \
  --p-max-depth 100000 \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/alpha-rarefaction.qzv

#----- Beta rarefaction to sequencing depth 28000 for jaccard|braycurtis|unweighted_unifrac|weighted_unifrac

mkdir qiime2_output/16S_endo_v34_2021/filtered/diversity/beta_rarefaction

qiime diversity beta-rarefaction \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --p-metric jaccard \
  --p-clustering-method upgma \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-sampling-depth 28000 \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/beta_rarefaction/jaccard-rarefied.qzv

qiime diversity beta-rarefaction \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --p-metric unweighted_unifrac \
  --p-clustering-method upgma \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-sampling-depth 28000 \
  --i-phylogeny qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-rooted-tree.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/beta_rarefaction/unweightedunifrac-rarefied.qzv

qiime diversity beta-rarefaction \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --p-metric braycurtis \
  --p-clustering-method upgma \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-sampling-depth 28000 \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/beta_rarefaction/braycurtis-rarefied.qzv

qiime diversity beta-rarefaction \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --p-metric weighted_unifrac \
  --p-clustering-method upgma \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-sampling-depth 28000 \
  --i-phylogeny qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-rooted-tree.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/beta_rarefaction/weightedunifrac-rarefied.qzv

#----- Calculate core metrics at depth 28000 -> D28000

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-rooted-tree.qza \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --p-sampling-depth 28000 \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --output-dir qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000

#----- Calculate alpha diversity group significance

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/shannon_vector.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/shannon_vector.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/faith_pd_vector.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/faith_pd_vector.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/evenness_vector.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/evenness_vector.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/observed_features_vector.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/observed_features_vector.qzv

#----- Compositional data analysis with DEICODE: robust Aitchison PCA (rPCA)

mkdir qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode

qiime deicode rpca \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 500 \
  --o-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-ordination.qza \
  --o-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix.qza

qiime emperor biplot \
    --i-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-ordination.qza \
    --m-sample-metadata-file master_metadata/16SendoV34.tsv \
    --m-feature-metadata-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
    --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-biplot.qzv \
    --p-number-of-features 8  

# Repeat rPCA with minimum sample count 1000 or more and then repeat PERMANOVA to test the effects of low read samples

qiime deicode rpca \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 1000 \
  --o-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-ordination-1000.qza \
  --o-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix-1000.qza

qiime emperor biplot \
    --i-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-ordination-1000.qza \
    --m-sample-metadata-file master_metadata/16SendoV34.tsv \
    --m-feature-metadata-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
    --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-biplot-1000.qzv \
    --p-number-of-features 8

qiime deicode rpca \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 4000 \
  --o-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-ordination-4000.qza \
  --o-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix-4000.qza

qiime emperor biplot \
    --i-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-ordination-4000.qza \
    --m-sample-metadata-file master_metadata/16SendoV34.tsv \
    --m-feature-metadata-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
    --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-biplot-4000.qzv \
    --p-number-of-features 8

qiime deicode rpca \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 28000 \
  --o-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-ordination-28000.qza \
  --o-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix-28000.qza

qiime emperor biplot \
    --i-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-ordination-28000.qza \
    --m-sample-metadata-file master_metadata/16SendoV34.tsv \
    --m-feature-metadata-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
    --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-biplot-28000.qzv \
    --p-number-of-features 8

#----- Filter samples with low reads <28000

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --p-min-frequency 28000 \
  --o-filtered-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qza

qiime feature-table summarize \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qzv \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv

#----- Filter tables based on sample sites

# Cloacal samples without low read samples
qiime feature-table filter-samples \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-where "[SampleSite]='CLOACA'" \
  --o-filtered-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal.qza

qiime feature-table summarize \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal.qzv \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv

# Oral samples without low read samples
qiime feature-table filter-samples \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-where "[SampleSite]='ORAL'" \
  --o-filtered-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-oral.qza

qiime feature-table summarize \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-oral.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-oral.qzv \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv

# Cloacal from all samples chl ands mitoh filtered (low read samples included)
qiime feature-table filter-samples \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-where "[SampleSite]='CLOACA'" \
  --o-filtered-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal-fltrd.qza

qiime feature-table summarize \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal-fltrd.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal-fltrd.qzv \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv

# Oral from all samples chl ands mitoh filtered (low read samples included)
qiime feature-table filter-samples \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-where "[SampleSite]='ORAL'" \
  --o-filtered-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-oral-fltrd.qza

qiime feature-table summarize \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-oral-fltrd.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-oral-fltrd.qzv \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv

# Water from all samples chl ands mitoh filtered
qiime feature-table filter-samples \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-where "[SampleSite]='TANK WATER'" \
  --o-filtered-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-water-fltrd.qza

qiime feature-table summarize \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-water-fltrd.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-water-fltrd.qzv \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv

# Core metrics D28000 cloacal
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-rooted-tree.qza \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal.qza \
  --p-sampling-depth 28000 \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --output-dir qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal

# Core metrics D28000 oral
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny qiime2_output/16S_endo_v34_2021/phylogeny/16Sendov34-rooted-tree.qza \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-oral.qza \
  --p-sampling-depth 28000 \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --output-dir qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral

# Robust Aitchison PCA on cloacal and oral samples (low read samples included)

mkdir qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal
mkdir qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral

qiime deicode rpca \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal-fltrd.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 500 \
  --o-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-ordination.qza \
  --o-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix.qza

qiime emperor biplot \
  --i-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-ordination.qza \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv \
  --m-feature-metadata-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-biplot.qzv \
  --p-number-of-features 8

qiime deicode rpca \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-oral-fltrd.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 500 \
  --o-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-ordination.qza \
  --o-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix.qza

qiime emperor biplot \
  --i-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-ordination.qza \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv \
  --m-feature-metadata-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-biplot.qzv \
  --p-number-of-features 8

# Robust Aitchison PCA on cloacal and oral samples (low read samples excluded)

qiime deicode rpca \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 28000 \
  --o-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-ordination-28000.qza \
  --o-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix-28000.qza

qiime emperor biplot \
  --i-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-ordination.qza \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv \
  --m-feature-metadata-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-biplot-28000.qzv \
  --p-number-of-features 8

qiime deicode rpca \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-oral.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 28000 \
  --o-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-ordination-28000.qza \
  --o-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix-28000.qza

qiime emperor biplot \
  --i-biplot qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-ordination.qza \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv \
  --m-feature-metadata-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-biplot-28000.qzv \
  --p-number-of-features 8

#----- PERMANOVA testing

# SampleSite clo vs- oral vs- water on complete samples

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column SampleSite \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/permanova-bray-samplesite.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column SampleSite \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/permanova-jaccard-samplesite.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column SampleSite \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/permanova-unweighted_unifrac-samplesite.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column SampleSite \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/permanova-weighted_unifrac-samplesite.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column SampleSite \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/permanova-rPCA-samplesite.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix-1000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column SampleSite \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/permanova-rPCA-samplesite-1000.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix-4000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column SampleSite \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/permanova-rPCA-samplesite-4000.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix-28000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column SampleSite \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/permanova-rPCA-samplesite-28000.qzv

# Permanova cloacal HospStatus|TurtleSex|AgeRange|DetSex

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-bray-HospStatus.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-jaccard-HospStatus.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-unweighted_unifrac-HospStatus.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-weighted_unifrac-HospStatus.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/permanova-rPCA-HospStatus.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix-28000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/permanova-rPCA-HospStatus-28000.qzv


qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-bray-TurtleSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-jaccard-TurtleSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-unweighted_unifrac-TurtleSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-weighted_unifrac-TurtleSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/permanova-rPCA-TurtleSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix-28000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/permanova-rPCA-TurtleSex-28000.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-bray-AgeRange.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-jaccard-AgeRange.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-unweighted_unifrac-AgeRange.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-weighted_unifrac-AgeRange.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/permanova-rPCA-AgeRange.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix-28000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/permanova-rPCA-AgeRange-28000.qzv


qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-bray-DetSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-jaccard-DetSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-unweighted_unifrac-DetSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-weighted_unifrac-DetSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/permanova-rPCA-DetSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix-28000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/permanova-rPCA-DetSex-28000.qzv


# testing additional categories Age2
qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-bray-Age2.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-jaccard-Age2.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-unweighted_unifrac-Age2.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-cloacal/permanova-weighted_unifrac-Age2.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/permanova-rPCA-Age2.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/rPCA-matrix-28000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-cloacal/permanova-rPCA-Age2-28000.qzv

# Permanova oral HospStatus|TurtleSex|AgeRange

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-bray-HospStatus.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-jaccard-HospStatus.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-unweighted_unifrac-HospStatus.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-weighted_unifrac-HospStatus.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/permanova-rPCA-HospStatus.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix-28000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column HospStatus \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/permanova-rPCA-HospStatus-28000.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-bray-TurtleSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-jaccard-TurtleSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-unweighted_unifrac-TurtleSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-weighted_unifrac-TurtleSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/permanova-rPCA-TurtleSex.qzv


qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix-28000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column TurtleSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/permanova-rPCA-TurtleSex-28000.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-bray-AgeRange.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-jaccard-AgeRange.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-unweighted_unifrac-AgeRange.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-weighted_unifrac-AgeRange.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/permanova-rPCA-AgeRange.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix-28000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column AgeRange \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/permanova-rPCA-AgeRange-28000.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-bray-DetSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-jaccard-DetSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-unweighted_unifrac-DetSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-weighted_unifrac-DetSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/permanova-rPCA-DetSex.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix-28000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column DetSex \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/permanova-rPCA-DetSex-28000.qzv


# Testing additional category Age2
qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-bray-Age2.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-jaccard-Age2.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-unweighted_unifrac-Age2.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000-oral/permanova-weighted_unifrac-Age2.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/permanova-rPCA-Age2.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/rPCA-matrix-28000.qza\
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --m-metadata-column Age2 \
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode-oral/permanova-rPCA-Age2-28000.qzv

#----- ADONIS in Qiime2

qiime diversity adonis \
  --i-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/bray_curtis_distance_matrix.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-formula SampleSite \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/adonis-bray-SampleSite.qzv
qiime diversity adonis \
  --i-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/jaccard_distance_matrix.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-formula SampleSite \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/adonis-jaccard-SampleSite.qzv
qiime diversity adonis \
  --i-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-formula SampleSite \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/adonis-unweighted_unifrac-SampleSite.qzv
qiime diversity adonis \
  --i-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-formula SampleSite \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/core-metrics-results-D28000/adonis-weighted_unifrac-SampleSite.qzv
qiime diversity adonis \
  --i-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix-28000.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-formula SampleSite \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/adonis-rPCA-SampleSite-28000.qzv

qiime diversity adonis \
  --i-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix-28000.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-formula "SampleSite + HospStatus" \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/adonis-rPCA-SampleSiteHosp-28000.qzv

qiime diversity adonis \
  --i-distance-matrix qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-matrix-28000.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-formula "SampleSite + TurtleSex" \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/adonis-rPCA-SampleSiteSex-28000.qzv

#----- Export merged data (table counts, rep seqs, taxonomy)

# original non-filtered reads and samples

qiime feature-table transpose \
  --i-table qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-table-dada2-trimmed.qza \
  --o-transposed-feature-table qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-table-dada2-trimmed-transposed.qza

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-rep-seqs-dada2-trimmed.qza \
  --m-input-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --m-input-file qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-table-dada2-trimmed-transposed.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/merged-all-data-no_filter.qzv

qiime tools export \
  --input-path qiime2_output/16S_endo_v34_2021/merged-all-data-no_filter.qzv \
  --output-path qiime2_output/16S_endo_v34_2021/merged-no_filter

# chloroplasts and mitochondria filtered

qiime feature-table transpose \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --o-transposed-feature-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-transposed.qza

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-rep-seqs-dada2-trimmed.qza \
  --m-input-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --m-input-file qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-transposed.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/merged-all-data-mitochl_filtered.qzv

qiime tools export \
  --input-path qiime2_output/16S_endo_v34_2021/merged-all-data-mitochl_filtered.qzv \
  --output-path qiime2_output/16S_endo_v34_2021/merged-filter


#----- Collapse taxonomy on different levels

mkdir qiime2_output/16S_endo_v34_2021/taxonomy/collapsed

for n in {1..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qza \
  --i-taxonomy qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}.qza; done

for n in {1..7};\
do qiime feature-table relative-frequency \
  --i-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}.qza \
  --o-relative-frequency-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-rel.qza; done

for n in {1..7};\
do qiime tools export \
  --input-path qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-rel.qza \
  --output-path qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-rel; done

for n in {1..7};\
do biom convert \
-i qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-rel/feature-table.biom \
-o qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/feature-table-${n}-rel.tsv \
--to-tsv; done

for n in {1..7};\
do qiime tools export \
  --input-path qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}.qza \
  --output-path qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}; done

for n in {1..7};\
do biom convert \
-i qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-rel/feature-table.biom \
-o qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/feature-table-${n}.tsv \
--to-tsv; done

# Unfiltered table

for n in {1..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-table-dada2-trimmed.qza \
  --i-taxonomy qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/original_table-${n}.qza; done

for n in {1..7};\
do qiime feature-table relative-frequency \
  --i-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/original_table-${n}.qza \
  --o-relative-frequency-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/original_table-${n}-rel.qza; done

for n in {1..7};\
do qiime tools export \
  --input-path qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/original_table-${n}-rel.qza \
  --output-path qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/original_table-${n}-rel; done

for n in {1..7};\
do biom convert \
-i qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/original_table-${n}-rel/feature-table.biom \
-o qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/original_table-${n}-rel.tsv \
--to-tsv; done

for n in {1..7};\
do qiime tools export \
  --input-path qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/original_table-${n}.qza \
  --output-path qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/original_table-${n}; done

for n in {1..7};\
do biom convert \
-i qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/original_table-${n}-rel/feature-table.biom \
-o qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/original_table-${n}.tsv \
--to-tsv; done

#----- Core features

mkdir qiime2_output/16S_endo_v34_2021/core_features

qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-no_mitochl-filter_low_read_samples-core_features.qzv

for n in {1..7};\
do qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/table-${n}-core_features.qzv; done

# investigate core features in SampleSite pairs and individual categories
qiime feature-table filter-samples \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-where "[SampleSite] IN ('CLOACA', 'ORAL')" \
  --o-filtered-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-clo-orl.qza &&
qiime feature-table filter-samples \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-where "[SampleSite] IN ('CLOACA', 'TANK WATER')" \
  --o-filtered-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-clo-wtr.qza &&
qiime feature-table filter-samples \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-where "[SampleSite] IN ('ORAL', 'TANK WATER')" \
  --o-filtered-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-orl-wtr.qza &&
qiime feature-table filter-samples \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl-filter_low_read_samples.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --p-where "[SampleSite] = 'TANK WATER'" \
  --o-filtered-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-water.qza

qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-cloacal-core_asv.qzv &&
qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-oral.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-oral-core_asv.qzv &&
qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-water.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-water-core_asv.qzv &&
qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-clo-orl.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-clo-orl-core_asv.qzv &&
qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-clo-wtr.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-clo-wtr-core_asv.qzv &&
qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-orl-wtr.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-orl-wtr-core_asv.qzv

# Core features in individual tables on collapsed taxonomy - only family and genus and species

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-cloacal.qza \
  --i-taxonomy qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-cloacal.qza &&
qiime taxa collapse \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-oral.qza \
  --i-taxonomy qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-oral.qza &&
qiime taxa collapse \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-water.qza \
  --i-taxonomy qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-water.qza &&
qiime taxa collapse \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-clo-orl.qza \
  --i-taxonomy qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-clo-orl.qza &&
qiime taxa collapse \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-clo-wtr.qza \
  --i-taxonomy qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-clo-wtr.qza &&
qiime taxa collapse \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-orl-wtr.qza \
  --i-taxonomy qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-orl-wtr.qza; done

for n in {5..7};\
do qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-cloacal.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-${n}-cloacal-core_features.qzv &&
qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-oral.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-${n}-oral-core_features.qzv &&
qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-water.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-${n}-water-core_features.qzv &&
qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-clo-orl.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-${n}-clo-orl-core_features.qzv &&
qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-clo-wtr.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-${n}-clo-wtr-core_features.qzv &&
qiime feature-table core-features \
  --i-table qiime2_output/16S_endo_v34_2021/taxonomy/collapsed/table-${n}-orl-wtr.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/core_features/16Sendov34-table-${n}-orl-wtr-core_features.qzv; done

#----- QURRO on rPCA loadings
qiime qurro loading-plot \
  --i-table qiime2_output/16S_endo_v34_2021/filtered/16Sendov34-table-no_mitochl.qza \
  --i-ranks qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/rPCA-ordination-28000.qza \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv \
  --m-feature-metadata-file qiime2_output/16S_endo_v34_2021/taxonomy/16Sendov34_taxonomy_silva_138_v34.qza \
  --o-visualization qiime2_output/16S_endo_v34_2021/filtered/diversity/deicode/qurro-plot-28000.qzv
