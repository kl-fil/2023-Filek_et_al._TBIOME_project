# -----------------------------------------------------------------------------# 
# QIIME2 ITS2 NGS data analyses for TurtleBIOME project (Bosak lab, University
# of Zagreb)
# Paper: TBD holobiont
# By: Borna Branimir Vuković
# Samples: cloacal samples and tank water of loggerhead sea turtles
# Objective: analyze endozoic and environmental ITS2 dataset
# Region: fungal ITS2 region of nucleic ribosomal gene 
# Primers: ITS3 and ITS4 (White et al. 1990)
# Platform: Illumina MiSeq v3 (300x2 bp paired-end)
# Sequences available at European Nucleotide Archive under accession:
# - PRJEB62762
# ZENODO data DOI: 10.5281/zenodo.8054926
# Analyses performed in Qiime2 v. 2021.8
# ---System versions---
# Python version: 3.8.10
# QIIME 2 release: 2021.8
# QIIME 2 version: 2021.8.0
# q2cli version: 2021.8.0
# ---Installed plugins---
# alignment: 2021.8.0
# composition: 2021.8.0
# cutadapt: 2021.8.0
# dada2: 2021.8.0
# deblur: 2021.8.0
# deicode: 0.2.4
# demux: 2021.8.0
# diversity: 2021.8.0
# diversity-lib: 2021.8.0
# emperor: 2021.8.0
# feature-classifier: 2021.8.0
# feature-table: 2021.8.0
# fragment-insertion: 2021.8.0
# gneiss: 2021.8.0
# longitudinal: 2021.8.0
# metadata: 2021.8.0
# phylogeny: 2021.8.0
# quality-control: 2021.8.0
# quality-filter: 2021.8.0
# sample-classifier: 2021.8.0
# taxa: 2021.8.0
# types: 2021.8.0
# vsearch: 2021.8.0
#-----------------------------------------------------------------------------#

# Note: prior to assigning taxonomy a classifier needs to be downloaded and put in classifier folder
# ITS2 classifier: UNITE version 8 dynamic classifier - https://unite.ut.ee/repository.php (or see Mendeley Data)

#----- QIIME 2 environment activation in Conda
conda activate qiime2-2021.8

#----- Change to working directory and set up folders

cd data_analyses
mkdir -p input_data/ITS2_endo_2021
# copy fastq.gz files from source (ENA) to input_data/ITS2_endo_2021
# the sequences were pre-processed by the sequencing company and most non-biological sequences were removed 

# create directory for qiime2 output
mkdir -p qiime2_output/ITS2_analysis-outputs

#----- Import ITS2 cloacal and tank water data in Cassava 1.8 format

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input_data/ITS2_endo_2021 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path qiime2_output/ITS2_analysis-outputs/demux-paired-end_ITS.qza

qiime demux summarize \
  --i-data qiime2_output/ITS2_analysis-outputs/demux-paired-end_ITS.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/demux-paired-end_ITS.qzv

#----- Remove primers in sequences using cutadapt prior to DADA2, avoiding readthrough due to sequence length variability
#for ITS2
#  --p-adapter-f GCATATCAATAAGCGGAGGA\ reverse primer reverse complement
#  --p-front-f GCATCGATGAAGAACGCAGC\ forward primer ITS3
#  --p-adapter-r GCTGCGTTCTTCATCGATGC\ forward primer reverse complement
#  --p-front-r TCCTCCGCTTATTGATATGC\ reverse primer ITS4

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences qiime2_output/ITS2_analysis-outputs/demux-paired-end_ITS.qza \
  --p-adapter-f GCATATCAATAAGCGGAGGA \
  --p-front-f GCATCGATGAAGAACGCAGC \
  --p-adapter-r GCTGCGTTCTTCATCGATGC \
  --p-front-r TCCTCCGCTTATTGATATGC \
  --o-trimmed-sequences qiime2_output/ITS2_analysis-outputs/demux-trimmed_ITS.qza \
  --verbose
  
qiime demux summarize \
  --i-data qiime2_output/ITS2_analysis-outputs/demux-trimmed_ITS.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/demux-trimmed_ITS.qzv
  
#----- Denoising sequences with DADA2

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2_output/ITS2_analysis-outputs/demux-trimmed_ITS.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 257 \
  --p-trunc-len-r 178 \
  --o-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
  --o-representative-sequences qiime2_output/ITS2_analysis-outputs/rep-seqs-trimmed.qza \
  --o-denoising-stats qiime2_output/ITS2_analysis-outputs/denoising-stats-trimmed.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file qiime2_output/ITS2_analysis-outputs/denoising-stats-trimmed.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/denoising-stats-trimmed.qzv

qiime feature-table summarize \
  --i-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/table-trimmed.qzv \
  --m-sample-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv
  
qiime feature-table tabulate-seqs \
  --i-data qiime2_output/ITS2_analysis-outputs/rep-seqs-trimmed.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/rep-seqs-trimmed.qzv
  
#----- Align sequences using mafft and fasttree

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences qiime2_output/ITS2_analysis-outputs/rep-seqs-trimmed.qza \
  --o-alignment qiime2_output/ITS2_analysis-outputs/aligned-rep-seqs.qza \
  --o-masked-alignment qiime2_output/ITS2_analysis-outputs/masked-aligned-rep-seqs.qza \
  --o-tree qiime2_output/ITS2_analysis-outputs/unrooted-tree.qza \
  --o-rooted-tree qiime2_output/ITS2_analysis-outputs/rooted-tree.qza

#----- Assigning taxonomy via UNITE database v.8 dynamic classifier (premade by Qiime 2 team)

qiime feature-classifier classify-sklearn \
  --i-classifier classifiers/unite_classifier_ver8_dynamic.qza \
  --i-reads qiime2_output/ITS2_analysis-outputs/rep-seqs-trimmed.qza \
  --o-classification qiime2_output/ITS2_analysis-outputs/taxonomy.qza

qiime metadata tabulate \
  --m-input-file qiime2_output/ITS2_analysis-outputs/taxonomy.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/taxonomy.qzv

qiime taxa barplot \
  --i-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
  --i-taxonomy qiime2_output/ITS2_analysis-outputs/taxonomy.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --o-visualization qiime2_output/ITS2_analysis-outputs/taxa-bar-plots.qzv

#----- Alpha rarefaction to determine approapriate depth for downstream analyses

qiime diversity alpha-rarefaction \
  --i-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
  --i-phylogeny qiime2_output/ITS2_analysis-outputs/rooted-tree.qza \
  --p-max-depth 69203 \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --o-visualization qiime2_output/ITS2_analysis-outputs/alpha-rarefaction.qzv

#----- Filter tables based on sample sites

qiime feature-table filter-samples \
  --i-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --p-where "[SampleSite]='CLOACA'" \
  --o-filtered-table qiime2_output/ITS2_analysis-outputs/cloacal-samples-filtered-table-trimmed.qza

qiime feature-table summarize \
  --i-table qiime2_output/ITS2_analysis-outputs/cloacal-samples-filtered-table-trimmed.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/cloacal-samples-filtered-table-trimmed.qzv \

qiime feature-table filter-samples \
  --i-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --p-where "[SampleSite]='WATER'" \
  --o-filtered-table qiime2_output/ITS2_analysis-outputs/water-samples-filtered-table-trimmed.qza

# remove negative control
qiime feature-table filter-samples \
  --i-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/samples-to-keep.tsv \
  --o-filtered-table qiime2_output/ITS2_analysis-outputs/no-neg-ctrl-filtered-table-trimmed.qza

#----- Calculate core metrics at depth 25077

# all samples
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny qiime2_output/ITS2_analysis-outputs/rooted-tree.qza \
  --i-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
  --p-sampling-depth 25077 \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --output-dir qiime2_output/ITS2_analysis-outputs/core-metrics-results

# cloacal samples
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny qiime2_output/ITS2_analysis-outputs/rooted-tree.qza \
  --i-table qiime2_output/ITS2_analysis-outputs/cloacal-samples-filtered-table-trimmed.qza \
  --p-sampling-depth 25077 \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --output-dir qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal

mv qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/bray_curtis_emperor.qzv qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/bray_curtis_emperor-cloacal.qzv
mv qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/jaccard_emperor.qzv qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/jaccard_emperor-cloacal.qzv
mv qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/unweighted_unifrac_emperor.qzv qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/unweighted_unifrac_emperor-cloacal.qzv
mv qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/weighted_unifrac_emperor.qzv qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/weighted_unifrac_emperor-cloacal.qzv  

# without negative control cloacal and water tank samples
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny qiime2_output/ITS2_analysis-outputs/rooted-tree.qza \
  --i-table qiime2_output/ITS2_analysis-outputs/no-neg-ctrl-filtered-table-trimmed.qza \
  --p-sampling-depth 25077 \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --output-dir qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl


#----- Core features

qiime feature-table core-features \
  --i-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-features.qzv

qiime feature-table core-features \
  --i-table qiime2_output/ITS2_analysis-outputs/cloacal-samples-filtered-table-trimmed.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/cloacal-core-features.qzv

qiime feature-table core-features \
  --i-table qiime2_output/ITS2_analysis-outputs/water-samples-filtered-table-trimmed.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/water-core-features.qzv


#----- Export merged data (table counts, rep seqs, taxonomy)

qiime feature-table transpose \
  --i-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
  --o-transposed-feature-table qiime2_output/ITS2_analysis-outputs/transposed-table.qza

qiime metadata tabulate \
  --m-input-file qiime2_output/ITS2_analysis-outputs/rep-seqs-trimmed.qza \
  --m-input-file qiime2_output/ITS2_analysis-outputs/taxonomy.qza \
  --m-input-file qiime2_output/ITS2_analysis-outputs/transposed-table.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/merged-table-seq.qzv

qiime tools export \
  --input-path qiime2_output/ITS2_analysis-outputs/merged-table-seq.qzv \
  --output-path qiime2_output/ITS2_analysis-outputs/merged-table-seq


#----- Compositional data analysis with DEICODE: robust Aitchison PCA (rPCA)

# all samples
qiime deicode rpca \
    --i-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
    --o-biplot qiime2_output/ITS2_analysis-outputs/ordination.qza \
    --o-distance-matrix qiime2_output/ITS2_analysis-outputs/distance-matrix.qza

qiime emperor biplot \
    --i-biplot qiime2_output/ITS2_analysis-outputs/ordination.qza \
    --m-sample-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
    --m-feature-metadata-file qiime2_output/ITS2_analysis-outputs/taxonomy.qza \
    --o-visualization qiime2_output/ITS2_analysis-outputs/raitchison-biplot.qzv \
    --p-number-of-features 5

# cloacal samples
qiime deicode rpca \
    --i-table qiime2_output/ITS2_analysis-outputs/cloacal-samples-filtered-table-trimmed.qza \
    --o-biplot qiime2_output/ITS2_analysis-outputs/cloacal-samples-ordination.qza \
    --o-distance-matrix qiime2_output/ITS2_analysis-outputs/cloacal-samples-distance-matrix.qza

qiime emperor biplot \
    --i-biplot qiime2_output/ITS2_analysis-outputs/cloacal-samples-ordination.qza \
    --m-sample-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
    --m-feature-metadata-file qiime2_output/ITS2_analysis-outputs/taxonomy.qza \
    --o-visualization qiime2_output/ITS2_analysis-outputs/cloacal-biplot.qzv \
    --p-number-of-features 5


#----- Calculate alpha diversity group significance

# cloacal samples
qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/faith_pd_vector.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/faith-pd-group-significance-cloacal.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/evenness_vector.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/evenness-group-significance-cloacal.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/shannon_vector.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/shannon-group-significance-cloacal.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/observed_features_vector.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/observed-features-group-significance-cloacal.qzv

# no negative control cloacal and water tank samples
qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/faith_pd_vector.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/faith-pd-group-significance-no-neg-ctrl.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/evenness_vector.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/evenness-group-significance-no-neg-ctrl.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/shannon_vector.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/shannon-group-significance-no-neg-ctrl.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/observed_features_vector.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/observed-features-group-significance-no-neg-ctrl.qzv

#----- PERMANOVA testing

# cloacal samples
qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/cloacal-samples-distance-matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column AgeRange \
  --o-visualization qiime2_output/ITS2_analysis-outputs/robust-aitchison-cloaca-age-range-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/cloacal-samples-distance-matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column HospitalizationLength \
  --o-visualization qiime2_output/ITS2_analysis-outputs/robust-aitchison-cloaca-hospitalization-length-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/bray_curtis_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column AgeRange \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/bray-curtis-age-range-significance-cloacal.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/bray_curtis_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column HospitalizationLength \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/bray-curtis-hospitalization-length-significance-cloacal.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/jaccard_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column AgeRange \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/jaccard-age-range-significance-cloacal.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/jaccard_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column HospitalizationLength \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/jaccard-hospitalization-length-significance-cloacal.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column AgeRange \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/weighted-unifrac-age-range-significance-cloacal.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column HospitalizationLength \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/weighted-unifrac-hospitalization-length-significance-cloacal.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column AgeRange \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/unweighted-unifrac-age-range-significance-cloacal.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column HospitalizationLength \
  --o-visualization qiime2_output/ITS2_analysis-outputs/qiime2_output/ITS2_analysis-outputs/core-metrics-results-cloacal/unweighted-unifrac-hospitalization-length-significance-cloacal.qzv \
  --p-pairwise

# no negative control cloacal and water tank samples
qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/bray_curtis_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column SampleSite \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/bray-curtis-sample-site-significance-no-neg-ctrl.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/jaccard_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column SampleSite \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/jaccard-sample-site-significance-no-neg-ctrl.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column SampleSite \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/weihgted-unifrac-sample-site-significance-no-neg-ctrl.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column SampleSite \
  --o-visualization qiime2_output/ITS2_analysis-outputs/core-metrics-results-no-neg-ctrl/unweighted-unifrac-sample-site-significance-no-neg-ctrl.qzv \
  --p-pairwise

qiime deicode rpca \
    --i-table qiime2_output/ITS2_analysis-outputs/no-neg-ctrl-filtered-table-trimmed.qza \
    --o-biplot qiime2_output/ITS2_analysis-outputs/no-neg-ctrl-ordination.qza \
    --o-distance-matrix qiime2_output/ITS2_analysis-outputs/no-neg-ctrl-distance-matrix.qza

qiime diversity beta-group-significance \
  --i-distance-matrix qiime2_output/ITS2_analysis-outputs/no-neg-ctrl-distance-matrix.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --m-metadata-column SampleSite \
  --o-visualization qiime2_output/ITS2_analysis-outputs/robust-aitchison-sample-site-significance-no-neg-ctrl.qzv \
  --p-pairwise

#----- Prepare data for ANCOM-BC2 in R (by Klara Filek)
mkdir qiime2_output/ITS2_analysis-outputs/ANCOMBC2

qiime feature-table filter-samples \
  --i-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --p-where '[SampleSite] IN ("CLOACA", "WATER")' \
  --o-filtered-table qiime2_output/ITS2_analysis-outputs/ANCOMBC2/filtered-table1.qza

qiime feature-table summarize \
  --i-table qiime2_output/ITS2_analysis-outputs/ANCOMBC2/filtered-table1.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/ANCOMBC2/filtered-table1.qzv

qiime feature-table filter-features \
  --i-table qiime2_output/ITS2_analysis-outputs/ANCOMBC2/filtered-table1.qza \
  --p-min-frequency 10 \
  --p-min-samples 2 \
  --o-filtered-table qiime2_output/ITS2_analysis-outputs/ANCOMBC2/feature-frequency-filtered-table.qza

qiime feature-table summarize \
  --i-table qiime2_output/ITS2_analysis-outputs/ANCOMBC2/feature-frequency-filtered-table.qza \
  --o-visualization qiime2_output/ITS2_analysis-outputs/ANCOMBC2/feature-frequency-filtered-table.qzv

qiime taxa barplot \
  --i-table qiime2_output/ITS2_analysis-outputs/ANCOMBC2/feature-frequency-filtered-table.qza \
  --i-taxonomy qiime2_output/ITS2_analysis-outputs/taxonomy.qza \
  --m-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
  --o-visualization qiime2_output/ITS2_analysis-outputs/ANCOMBC2/taxa-bar-plots-filtered.qzv

for n in {1..7};\
do qiime taxa collapse \
  --i-table qiime2_output/ITS2_analysis-outputs/ANCOMBC2/feature-frequency-filtered-table.qza \
  --i-taxonomy qiime2_output/ITS2_analysis-outputs/taxonomy.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/ITS2_analysis-outputs/ANCOMBC2/table-lvl-${n}.qza; done

# EMPress visualization
qiime empress community-plot \
    --i-tree qiime2_output/ITS2_analysis-outputs/rooted-tree.qza \
    --i-feature-table qiime2_output/ITS2_analysis-outputs/table-trimmed.qza \
    --m-sample-metadata-file qiime2_output/ITS2_analysis-outputs/metadata-ITS-new.tsv \
    --m-feature-metadata-file qiime2_output/ITS2_analysis-outputs/taxonomy.qza \
    --o-visualization qiime2_output/ITS2_analysis-outputs/empress-tree.qzv