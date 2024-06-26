# -----------------------------------------------------------------------------#
# QIIME2 16S NGS data analyses for TurtleBIOME project (Bosak lab, University
# of Zagreb)
# Paper: TBD holobiont
# By: Klara Filek
# Samples: cloacal & oral swabs, tank water, carapace samples of loggerhead sea turtles
# Objective: all sequences of multiple datasets trimmed to V4 region and analyzed together
# Region: V3-V4 & V4 region of 16S rRNA gene
# Primers: V3-V4 - 341F and 805R (Klindworth et al. 2013)
#          V4 - 515F and 806R (Apprill et al. 2015; Parada et al. 2016)       
# Platform: V3-V4 - Illumina MiSeq v3 (300x2 bp paired-end)                   
#           V4 - Illumina MiSeq v2 (250x2 bp paired-end)                      
# Epizoic Sequences available at European Nucleotide Archive under accessions:        
# - Kanjer et al. (2022) 16S V4 epizoic samples: PRJEB51458                
# - Kanjer et al. (2024) 16S V3-V4 epizoic samples: TBD                       
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
# v4 SILVA classifier: SILVA 138 99% OTUs from 515F/806R region of sequences - https://docs.qiime2.org/2023.2/data-resources/ (or see Mendeley Data)

#----- QIIME 2 environment activation in Conda
conda activate qiime2-2023.2

#----- Change to working directory and set up folders
cd data_analyses
mkdir -p input_data/{16S_epi_v4_kanjeretal2022,16S_epi_v34_2022}
# copy fastq.gz files from source (ENA) to input_data/16S_epi_v4_kanjeretal2022 and input_data/16S_epi_v34_2022

mkdir -p qiime2_output/{16S_epi_v4_kanjeretal2022,16S_epi_v34_2022}/{demux-dada2,filtered,taxonomy,phylogeny}

# #----- Rename original files for Cassava format - add "L001_" and remove "_trimmed"

# cd input_data/16S_epi_v34_2022/trimmed && for text in *_R{1,2}_001_trimmed.fastq.gz; do mv "$text" $(echo "$text" | sed 's/_R/_L001_R/g'); done
# for text in *_R{1,2}_001_trimmed.fastq.gz; do mv "$text" $(echo "$text" | sed 's/_trimmed//g'); done
# cd ../../..

#----- Import 16S V34 and V4 carapace data in Cassava 1.8 format

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input_data/16S_epi_v34_2022 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-demux-paired-end.qza

qiime demux summarize \
  --i-data qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-demux-paired-end.qza \
  --o-visualization qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-demux-paired-end.qzv

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path input_data/16S_epi_v4_kanjeretal2022 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-demux-paired-end.qza

qiime demux summarize \
  --i-data qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-demux-paired-end.qza \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-demux-paired-end.qzv


#----- Remove primers in sequences using cutadapt prior to DADA2, avoiding readthrough due to seq length variability
# for V34
#  --p-adapter-f GGATTAGATACCCBDGTAGTC\ reverse primer reverse complement
#  --p-front-f CCTACGGGNGGCWGCAG\ forward primer 341F
#  --p-adapter-r CTGCWGCCNCCCGTAGG\ forward primer reverse complement
#  --p-front-r GACTACHVGGGTATCTAATCC\ reverse primer 805R

# for V4
#  --p-adapter-f ATTAGAWACCCBNGTAGTCC\ reverse primer reverse complement
#  --p-front-f GTGYCAGCMGCCGCGGTAA\ forward primer 515F
#  --p-adapter-r TTACCGCGGCKGCTGRCAC\ forward primer reverse complement
#  --p-front-r GGACTACNVGGGTWTCTAAT\ reverse primer 806R

qiime cutadapt trim-paired \
  --p-cores 4 \
  --i-demultiplexed-sequences qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-demux-paired-end.qza \
  --p-adapter-f GGATTAGATACCCBDGTAGTC\
  --p-front-f CCTACGGGNGGCWGCAG\
  --p-adapter-r CTGCWGCCNCCCGTAGG\
  --p-front-r GACTACHVGGGTATCTAATCC\
  --o-trimmed-sequences qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-demux-trimmed.qza \
  --verbose

qiime demux summarize \
  --i-data qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-demux-trimmed.qza \
  --o-visualization qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-demux-trimmed.qzv

qiime cutadapt trim-paired \
  --p-cores 4 \
  --i-demultiplexed-sequences qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-demux-paired-end.qza \
  --p-adapter-f ATTAGAWACCCBNGTAGTCC\
  --p-front-f GTGYCAGCMGCCGCGGTAA\
  --p-adapter-r TTACCGCGGCKGCTGRCAC\
  --p-front-r GGACTACNVGGGTWTCTAAT\
  --o-trimmed-sequences qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-demux-trimmed.qza \
  --verbose

qiime demux summarize \
  --i-data qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-demux-trimmed.qza \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-demux-trimmed.qzv

#----- Denoising sequences with DADA2

# Non trimmed V4 carapace data
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-demux-paired-end.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-representative-sequences qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-rep-seqs-dada2.qza \
  --o-table qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-table-dada2.qza \
  --o-denoising-stats qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-stats-dada2.qza \
  --p-n-threads 0

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-stats-dada2.qza \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-stats-dada2.qzv

qiime feature-table tabulate-seqs \
  --i-data qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-rep-seqs-dada2.qza \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-rep-seqs-dada2.qzv

qiime feature-table summarize \
  --i-table qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-table-dada2.qza \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-table-dada2.qzv \
  --m-sample-metadata-file master_metadata/16SepiV4.tsv

# V34 carapace data
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-demux-trimmed.qza \
  --p-trunc-len-f 277 \
  --p-trunc-len-r 228 \
  --o-representative-sequences qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-rep-seqs-dada2-trimmed.qza \
  --o-table qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-table-dada2-trimmed.qza \
  --o-denoising-stats qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-stats-dada2-trimmed.qza \
  --p-n-threads 0

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-stats-dada2-trimmed.qza \
  --o-visualization qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-stats-dada2-trimmed.qzv

qiime feature-table tabulate-seqs \
  --i-data qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-rep-seqs-dada2-trimmed.qza \
  --o-visualization qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-rep-seqs-dada2-trimmed.qzv

qiime feature-table summarize \
  --i-table qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-table-dada2-trimmed.qza \
  --o-visualization qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-table-dada2-trimmed.qzv \
  --m-sample-metadata-file master_metadata/16SepiV34.tsv


#----- Assigning taxonomy via SILVA v.138 on full length epizoic sequences

qiime feature-classifier classify-sklearn \
  --i-classifier classifier/silva-138-99-515-806-nb-classifier.qza \
  --p-reads-per-batch 10000 \
  --p-n-jobs -2 \
  --i-reads qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-rep-seqs-dada2.qza \
  --o-classification qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxonomy_silva_138_v4.qza

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxonomy_silva_138_v4.qza \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxonomy_silva_138_v4.qzv

qiime taxa barplot \
  --i-table qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-table-dada2.qza \
  --i-taxonomy qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxonomy_silva_138_v4.qza \
  --m-metadata-file master_metadata/16SepiV4.tsv \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxa-bar-plots.qzv

qiime feature-classifier classify-sklearn \
  --i-classifier classifier/classifier-V34-LK.qza \
  --p-reads-per-batch 10000 \
  --p-n-jobs -2 \
  --i-reads qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-rep-seqs-dada2-trimmed.qza \
  --o-classification qiime2_output/16S_epi_v34_2022/taxonomy/16Sepiv34_taxonomy_silva_138_v34.qza

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_epi_v34_2022/taxonomy/16Sepiv34_taxonomy_silva_138_v34.qza \
  --o-visualization qiime2_output/16S_epi_v34_2022/taxonomy/16Sepiv34_taxonomy_silva_138_v34.qzv

qiime taxa barplot \
  --i-table qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-table-dada2-trimmed.qza \
  --i-taxonomy qiime2_output/16S_epi_v34_2022/taxonomy/16Sepiv34_taxonomy_silva_138_v34.qza \
  --m-metadata-file master_metadata/16SepiV34.tsv \
  --o-visualization qiime2_output/16S_epi_v34_2022/taxonomy/16Sepiv34_taxa-bar-plots.qzv

#----- Align sequences using mafft and fasttree on full length epizoic sequences

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-rep-seqs-dada2.qza \
  --p-n-threads auto \
  --o-alignment qiime2_output/16S_epi_v4_kanjeretal2022/phylogeny/16Sepiv4-aligned-rep-seqs.qza \
  --o-masked-alignment qiime2_output/16S_epi_v4_kanjeretal2022/phylogeny/16Sepiv4-masked-aligned-rep-seqs.qza \
  --o-tree qiime2_output/16S_epi_v4_kanjeretal2022/phylogeny/16Sepiv4-unrooted-tree.qza \
  --o-rooted-tree qiime2_output/16S_epi_v4_kanjeretal2022/phylogeny/16Sepiv4-rooted-tree.qza

qiime empress tree-plot \
    --i-tree qiime2_output/16S_epi_v4_kanjeretal2022/phylogeny/16Sepiv4-rooted-tree.qza \
    --m-feature-metadata-file qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxonomy_silva_138_v4.qza \
    --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/phylogeny/16Sepiv4-rooted-tree-viz-empress.qzv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-rep-seqs-dada2-trimmed.qza \
  --p-n-threads auto \
  --o-alignment qiime2_output/16S_epi_v34_2022/phylogeny/16Sepiv34-aligned-rep-seqs.qza \
  --o-masked-alignment qiime2_output/16S_epi_v34_2022/phylogeny/16Sepiv34-masked-aligned-rep-seqs.qza \
  --o-tree qiime2_output/16S_epi_v34_2022/phylogeny/16Sepiv34-unrooted-tree.qza \
  --o-rooted-tree qiime2_output/16S_epi_v34_2022/phylogeny/16Sepiv34-rooted-tree.qza

qiime empress tree-plot \
    --i-tree qiime2_output/16S_epi_v34_2022/phylogeny/16Sepiv34-rooted-tree.qza \
    --m-feature-metadata-file qiime2_output/16S_epi_v34_2022/taxonomy/16Sepiv34_taxonomy_silva_138_v34.qza \
    --o-visualization qiime2_output/16S_epi_v34_2022/phylogeny/16Sepiv34-rooted-tree-viz-empress.qzv


#----- Cut V34 and V4 data sequences to V4 region using cutadapt (endozoic and epizoic dataset)

mkdir -p qiime2_output/{16S_endo_trimv4_2021,16S_epi_trimv4_2022}/{dada2,taxonomy}

qiime cutadapt trim-paired \
  --p-cores 4 \
  --i-demultiplexed-sequences qiime2_output/16S_endo_v34_2021/demux-dada2/16Sendov34-demux-paired-end.qza \
  --p-adapter-f ATTAGAWACCCBNGTAGTCC\
  --p-front-f GTGYCAGCMGCCGCGGTAA\
  --p-adapter-r TTACCGCGGCKGCTGRCAC\
  --p-front-r GGACTACNVGGGTWTCTAAT\
  --o-trimmed-sequences qiime2_output/16S_endo_trimv4_2021/dada2/16Sendov4-demux-trimmed.qza \
  --verbose

qiime demux summarize \
  --i-data qiime2_output/16S_endo_trimv4_2021/dada2/16Sendov4-demux-trimmed.qza \
  --o-visualization qiime2_output/16S_endo_trimv4_2021/dada2/16Sendov4-demux-trimmed.qzv

qiime cutadapt trim-paired \
  --p-cores 4 \
  --i-demultiplexed-sequences qiime2_output/16S_epi_v34_2022/demux-dada2/16Sepiv34-demux-paired-end.qza \
  --p-adapter-f ATTAGAWACCCBNGTAGTCC\
  --p-front-f GTGYCAGCMGCCGCGGTAA\
  --p-adapter-r TTACCGCGGCKGCTGRCAC\
  --p-front-r GGACTACNVGGGTWTCTAAT\
  --o-trimmed-sequences qiime2_output/16S_epi_trimv4_2022/dada2/16Sepiv42022-demux-trimmed.qza \
  --verbose

qiime demux summarize \
  --i-data qiime2_output/16S_epi_trimv4_2022/dada2/16Sepiv42022-demux-trimmed.qza \
  --o-visualization qiime2_output/16S_epi_trimv4_2022/dada2/16Sepiv42022-demux-trimmed.qzv

#----- Denoising sequences with DADA2
# Trimmed from V34 to V4 cloacal, oral, tank water samples 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2_output/16S_endo_trimv4_2021/dada2/16Sendov4-demux-trimmed.qza \
  --p-trunc-len-f 101 \
  --p-trunc-len-r 228 \
  --o-representative-sequences qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-rep-seqs-dada2.qza \
  --o-table qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-table-dada2.qza \
  --o-denoising-stats qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-stats-dada2.qza \
  --p-n-threads 0

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-stats-dada2.qza \
  --o-visualization qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-stats-dada2.qzv

qiime feature-table tabulate-seqs \
  --i-data qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-rep-seqs-dada2.qza \
  --o-visualization qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-rep-seqs-dada2.qzv

qiime feature-table summarize \
  --i-table qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-table-dada2.qza \
  --o-visualization qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-table-dada2.qzv \
  --m-sample-metadata-file master_metadata/16SendoV34.tsv

# Trimmed to V4 carapace samples - original V4 samples trimmed to same length as V34 samples when trimmed to V4 region

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-demux-trimmed.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-trim-left-r 2 \
  --o-representative-sequences qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-rep-seqs-dada2-trimmed.qza \
  --o-table qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-table-dada2-trimmed.qza \
  --o-denoising-stats qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-stats-dada2-trimmed.qza \
  --p-n-threads 0

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-stats-dada2-trimmed.qza \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-stats-dada2-trimmed.qzv

qiime feature-table tabulate-seqs \
  --i-data qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-rep-seqs-dada2-trimmed.qza \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-rep-seqs-dada2-trimmed.qzv

qiime feature-table summarize \
  --i-table qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-table-dada2-trimmed.qza \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-table-dada2-trimmed.qzv \
  --m-sample-metadata-file master_metadata/16SepiV4.tsv

# Trimmed from V34 to V4 carapace samples 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2_output/16S_epi_trimv4_2022/dada2/16Sepiv42022-demux-trimmed.qza \
  --p-trunc-len-f 90 \
  --p-trunc-len-r 200 \
  --o-representative-sequences qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-rep-seqs-dada2.qza \
  --o-table qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-table-dada2.qza \
  --o-denoising-stats qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-stats-dada2.qza \
  --p-n-threads 0

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-stats-dada2.qza \
  --o-visualization qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-stats-dada2.qzv

qiime feature-table tabulate-seqs \
  --i-data qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-rep-seqs-dada2.qza \
  --o-visualization qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-rep-seqs-dada2.qzv

qiime feature-table summarize \
  --i-table qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-table-dada2.qza \
  --o-visualization qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-table-dada2.qzv \
  --m-sample-metadata-file master_metadata/16SepiV34.tsv


#----- Assigning taxonomy via SILVA v.138 for V4 trimmed sequences
# Cloacal, oral, tank water sequences

qiime feature-classifier classify-sklearn \
  --i-classifier classifier/silva-138-99-515-806-nb-classifier.qza \
  --p-reads-per-batch 10000 \
  --p-n-jobs -2 \
  --i-reads qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-rep-seqs-dada2.qza \
  --o-classification qiime2_output/16S_endo_trimv4_2021/taxonomy/16Sendotrimv4_taxonomy_silva_138_v4.qza &&

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_endo_trimv4_2021/taxonomy/16Sendotrimv4_taxonomy_silva_138_v4.qza \
  --o-visualization qiime2_output/16S_endo_trimv4_2021/taxonomy/16Sendotrimv4_taxonomy_silva_138_v4.qzv &&

qiime taxa barplot \
  --i-table qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-table-dada2.qza \
  --i-taxonomy qiime2_output/16S_endo_trimv4_2021/taxonomy/16Sendotrimv4_taxonomy_silva_138_v4.qza \
  --m-metadata-file master_metadata/16SendoV34.tsv \
  --o-visualization qiime2_output/16S_endo_trimv4_2021/taxonomy/16Sendotrimv4_taxa-bar-plots.qzv

# Carapace sequences

qiime feature-classifier classify-sklearn \
  --i-classifier classifier/silva-138-99-515-806-nb-classifier.qza \
  --p-reads-per-batch 10000 \
  --p-n-jobs -2 \
  --i-reads qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-rep-seqs-dada2.qza \
  --o-classification qiime2_output/16S_epi_trimv4_2022/taxonomy/16Sepitrimv4_taxonomy_silva_138_v4.qza &&

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_epi_trimv4_2022/taxonomy/16Sepitrimv4_taxonomy_silva_138_v4.qza \
  --o-visualization qiime2_output/16S_epi_trimv4_2022/taxonomy/16Sepitrimv4_taxonomy_silva_138_v4.qzv &&

qiime taxa barplot \
  --i-table qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-table-dada2.qza \
  --i-taxonomy qiime2_output/16S_epi_trimv4_2022/taxonomy/16Sepitrimv4_taxonomy_silva_138_v4.qza \
  --m-metadata-file master_metadata/16SepiV34.tsv \
  --o-visualization qiime2_output/16S_epi_trimv4_2022/taxonomy/16Sepitrimv4_taxa-bar-plots.qzv

qiime feature-classifier classify-sklearn \
  --i-classifier classifier/silva-138-99-515-806-nb-classifier.qza \
  --p-reads-per-batch 10000 \
  --p-n-jobs -2 \
  --i-reads qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-rep-seqs-dada2-trimmed.qza \
  --o-classification qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxonomy_silva_138_v4-trimmed.qza &&

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxonomy_silva_138_v4-trimmed.qza \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxonomy_silva_138_v4-trimmed.qzv &&

qiime taxa barplot \
  --i-table qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-table-dada2-trimmed.qza \
  --i-taxonomy qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxonomy_silva_138_v4-trimmed.qza \
  --m-metadata-file master_metadata/16SepiV4.tsv \
  --o-visualization qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxa-bar-plots-trimmed.qzv

#----- Merge trimmed V4 datasets, filter samples with <1000 reads

mkdir -pv qiime2_output/16S_v4_merged

qiime feature-table merge \
  --i-tables qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-table-dada2.qza \
  --i-tables qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-table-dada2.qza \
  --i-tables qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-table-dada2-trimmed.qza \
  --o-merged-table qiime2_output/16S_v4_merged/table-16S_v4_merged.qza

qiime feature-table merge-seqs \
  --i-data qiime2_output/16S_endo_trimv4_2021/dada2/16Sendotrimv4-rep-seqs-dada2.qza \
  --i-data qiime2_output/16S_epi_trimv4_2022/dada2/16Sepitrimv4-rep-seqs-dada2.qza \
  --i-data qiime2_output/16S_epi_v4_kanjeretal2022/demux-dada2/16Sepiv4-rep-seqs-dada2-trimmed.qza \
  --o-merged-data qiime2_output/16S_v4_merged/rep-seqs-16S_v4_merged.qza

qiime feature-table summarize \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged.qza \
  --o-visualization qiime2_output/16S_v4_merged/table-16S_v4_merged.qzv \
  --m-sample-metadata-file master_metadata/16S_all_metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data qiime2_output/16S_v4_merged/rep-seqs-16S_v4_merged.qza \
  --o-visualization qiime2_output/16S_v4_merged/rep-seqs-16S_v4_merged.qzv

qiime feature-table merge-taxa \
  --i-data  qiime2_output/16S_endo_trimv4_2021/taxonomy/16Sendotrimv4_taxonomy_silva_138_v4.qza \
  --i-data  qiime2_output/16S_epi_trimv4_2022/taxonomy/16Sepitrimv4_taxonomy_silva_138_v4.qza \
  --i-data  qiime2_output/16S_epi_v4_kanjeretal2022/taxonomy/16Sepiv4_taxonomy_silva_138_v4-trimmed.qza \
  --o-merged-data qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --o-visualization qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qzv

# Merged taxa bar plot

qiime taxa barplot \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/taxa-bar-plots-16S_v4_merged.qzv

# Merged alignment

mkdir qiime2_output/16S_v4_merged/phylogeny

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences qiime2_output/16S_v4_merged/rep-seqs-16S_v4_merged.qza \
  --p-n-threads auto \
  --o-alignment qiime2_output/16S_v4_merged/phylogeny/16S_v4_merged-aligned-rep-seqs.qza \
  --o-masked-alignment qiime2_output/16S_v4_merged/phylogeny/16S_v4_merged-masked-aligned-rep-seqs.qza \
  --o-tree qiime2_output/16S_v4_merged/phylogeny/16S_v4_merged-unrooted-tree.qza \
  --o-rooted-tree qiime2_output/16S_v4_merged/phylogeny/16S_v4_merged-rooted-tree.qza &&

qiime empress tree-plot \
  --i-tree qiime2_output/16S_v4_merged/phylogeny/16S_v4_merged-rooted-tree.qza \
  --m-feature-metadata-file qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --o-visualization qiime2_output/16S_v4_merged/phylogeny/16S_v4_merged-rooted-tree-viz-empress.qzv

#----- Calculate diversity metrics for trimmed V4 data - alpha, beta diversity

mkdir qiime2_output/16S_v4_merged/diversity

# Filter based on taxonomy

qiime taxa filter-table \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-exclude mitochondria,chloroplast,eukaryota,unassigned \
  --o-filtered-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza

qiime feature-table summarize \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --o-visualization qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qzv \
  --m-sample-metadata-file master_metadata/16S_all_metadata.tsv

# new taxa bar plot
qiime taxa barplot \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/taxa-bar-plots-16S_v4_merged-no-mitochleuk.qzv

# Alpha rarefaction

qiime diversity alpha-rarefaction \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --i-phylogeny qiime2_output/16S_v4_merged/phylogeny/16S_v4_merged-rooted-tree.qza \
  --p-max-depth 100000 \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/diversity/alpha-rarefaction.qzv

# Beta rarefaction na depth 29000 jaccard|braycurtis|unweighted_unifrac|weighted_unifrac

mkdir qiime2_output/16S_v4_merged/diversity/beta_rarefaction

qiime diversity beta-rarefaction \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-metric jaccard \
  --p-clustering-method upgma \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-sampling-depth 29000 \
  --o-visualization qiime2_output/16S_v4_merged/diversity/beta_rarefaction/jaccard-rarefied.qzv &&

qiime diversity beta-rarefaction \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-metric unweighted_unifrac \
  --p-clustering-method upgma \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-sampling-depth 29000 \
  --i-phylogeny qiime2_output/16S_v4_merged/phylogeny/16S_v4_merged-rooted-tree.qza \
  --o-visualization qiime2_output/16S_v4_merged/diversity/beta_rarefaction/unweightedunifrac-rarefied.qzv &&

qiime diversity beta-rarefaction \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-metric braycurtis \
  --p-clustering-method upgma \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-sampling-depth 29000 \
  --o-visualization qiime2_output/16S_v4_merged/diversity/beta_rarefaction/braycurtis-rarefied.qzv &&

qiime diversity beta-rarefaction \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-metric weighted_unifrac \
  --p-clustering-method upgma \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-sampling-depth 29000 \
  --i-phylogeny qiime2_output/16S_v4_merged/phylogeny/16S_v4_merged-rooted-tree.qza \
  --o-visualization qiime2_output/16S_v4_merged/diversity/beta_rarefaction/weightedunifrac-rarefied.qzv

# Calculate core metrics at depth 29000 -> D29000

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny qiime2_output/16S_v4_merged/phylogeny/16S_v4_merged-rooted-tree.qza \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-sampling-depth 29000 \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --output-dir qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000

# Calculate core metrics at depth 100000 -> D100000

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny qiime2_output/16S_v4_merged/phylogeny/16S_v4_merged-rooted-tree.qza \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-sampling-depth 100000 \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --output-dir qiime2_output/16S_v4_merged/diversity/core-metrics-results-D100000

# Alpha diversity group significance

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/shannon_vector.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/shannon_vector.qzv &&

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/faith_pd_vector.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/faith_pd_vector.qzv &&

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/evenness_vector.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/evenness_vector.qzv &&

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/observed_features_vector.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/observed_features_vector.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_v4_merged/diversity/core-metrics-results-D100000/shannon_vector.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D100000/shannon_vector.qzv &&

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_v4_merged/diversity/core-metrics-results-D100000/faith_pd_vector.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D100000/faith_pd_vector.qzv &&

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_v4_merged/diversity/core-metrics-results-D100000/evenness_vector.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D100000/evenness_vector.qzv &&

qiime diversity alpha-group-significance \
  --i-alpha-diversity qiime2_output/16S_v4_merged/diversity/core-metrics-results-D100000/observed_features_vector.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D100000/observed_features_vector.qzv


# DEICODE - robust Aitchison PCA - variable sample filtering and feature filtering

mkdir qiime2_output/16S_v4_merged/diversity/deicode

qiime deicode rpca \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 500 \
  --o-biplot qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination.qza \
  --o-distance-matrix qiime2_output/16S_v4_merged/diversity/deicode/rPCA-matrix.qza

qiime emperor biplot \
    --i-biplot qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination.qza \
    --m-sample-metadata-file master_metadata/16S_all_metadata.tsv \
    --m-feature-metadata-file qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
    --o-visualization qiime2_output/16S_v4_merged/diversity/deicode/rPCA-biplot.qzv \
    --p-number-of-features 8

qiime deicode rpca \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 1000 \
  --o-biplot qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination-1000.qza \
  --o-distance-matrix qiime2_output/16S_v4_merged/diversity/deicode/rPCA-matrix-1000.qza

qiime emperor biplot \
  --i-biplot qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination-1000.qza \
  --m-sample-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-feature-metadata-file qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --o-visualization qiime2_output/16S_v4_merged/diversity/deicode/rPCA-biplot-1000.qzv \
  --p-number-of-features 8 

qiime deicode rpca \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 29000 \
  --o-biplot qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination-29000.qza \
  --o-distance-matrix qiime2_output/16S_v4_merged/diversity/deicode/rPCA-matrix-29000.qza

qiime emperor biplot \
  --i-biplot qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination-29000.qza \
  --m-sample-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-feature-metadata-file qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --o-visualization qiime2_output/16S_v4_merged/diversity/deicode/rPCA-biplot-29000.qzv \
  --p-number-of-features 8

qiime deicode rpca \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 100000 \
  --o-biplot qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination-100000.qza \
  --o-distance-matrix qiime2_output/16S_v4_merged/diversity/deicode/rPCA-matrix-100000.qza

qiime emperor biplot \
  --i-biplot qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination-100000.qza \
  --m-sample-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-feature-metadata-file qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --o-visualization qiime2_output/16S_v4_merged/diversity/deicode/rPCA-biplot-100000.qzv \
  --p-number-of-features 8

#----- PERMANOVA

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/bray_curtis_distance_matrix.qza\
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-metadata-column SampleSite\
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/permanova-bray-samplesite.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/jaccard_distance_matrix.qza\
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-metadata-column SampleSite\
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/permanova-jaccard-samplesite.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/unweighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-metadata-column SampleSite\
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/permanova-unweighted_unifrac-samplesite.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/weighted_unifrac_distance_matrix.qza\
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-metadata-column SampleSite\
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_v4_merged/diversity/core-metrics-results-D29000/permanova-weighted_unifrac-samplesite.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_v4_merged/diversity/deicode/rPCA-matrix-29000.qza\
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-metadata-column SampleSite\
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_v4_merged/diversity/deicode/permanova-29000-samplesite.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_v4_merged/diversity/deicode/rPCA-matrix-1000.qza\
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-metadata-column SampleSite\
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_v4_merged/diversity/deicode/permanova-1000-samplesite.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix  qiime2_output/16S_v4_merged/diversity/deicode/rPCA-matrix-100000.qza\
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-metadata-column SampleSite\
  --p-pairwise TRUE \
  --o-visualization qiime2_output/16S_v4_merged/diversity/deicode/permanova-100000-samplesite.qzv

#----- Export merged data (table counts, rep seqs, taxonomy) for trimmed to V4 sequences

## Not filtered
qiime feature-table transpose \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged.qza \
  --o-transposed-feature-table qiime2_output/16S_v4_merged/table-16S_v4_merged-transposed.qza &&

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_v4_merged/rep-seqs-16S_v4_merged.qza \
  --m-input-file qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --m-input-file qiime2_output/16S_v4_merged/table-16S_v4_merged-transposed.qza \
  --o-visualization qiime2_output/16S_v4_merged/merged-all-data.qzv &&

qiime tools export \
  --input-path qiime2_output/16S_v4_merged/merged-all-data.qzv \
  --output-path qiime2_output/16S_v4_merged/merged-no_filter

# Filtered chloroplasts, mitochondria, eukaryotes and unassigned

qiime feature-table transpose \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --o-transposed-feature-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-transposed.qza

qiime metadata tabulate \
  --m-input-file qiime2_output/16S_v4_merged/rep-seqs-16S_v4_merged.qza \
  --m-input-file qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --m-input-file qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-transposed.qza \
  --o-visualization qiime2_output/16S_v4_merged/merged-all-data-filtered.qzv

qiime tools export \
  --input-path qiime2_output/16S_v4_merged/merged-all-data-filtered.qzv \
  --output-path qiime2_output/16S_v4_merged/merged-filtered


#----- Filter samples with low reads (<1000)

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-min-frequency 1000 \
  --o-filtered-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-filter_low_read_samples.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk.qza \
  --p-min-frequency 29000 \
  --o-filtered-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza

qiime feature-table summarize \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --o-visualization qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qzv \
  --m-sample-metadata-file master_metadata/16S_all_metadata.tsv

# new taxa bar plot

qiime taxa barplot \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --o-visualization qiime2_output/16S_v4_merged/taxa-bar-plots-16S_v4_merged-no-mitochleuk-D29000.qzv


# Collapse taxonomy

mkdir qiime2_output/16S_V4_merged/taxonomy_collapsed

for n in {1..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}.qza; done

for n in {1..7};\
do qiime feature-table relative-frequency \
  --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}.qza \
  --o-relative-frequency-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-rel.qza; done

for n in {1..7};\
do qiime tools export \
  --input-path qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-rel.qza \
  --output-path qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-rel; done

for n in {1..7};\
do biom convert \
-i qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-rel/feature-table.biom \
-o qiime2_output/16S_v4_merged/taxonomy_collapsed/feature-table-${n}-rel.tsv \
--to-tsv; done

for n in {1..7};\
do qiime tools export \
  --input-path qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}.qza \
  --output-path qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}; done

for n in {1..7};\
do biom convert \
-i qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}/feature-table.biom \
-o qiime2_output/16S_v4_merged/taxonomy_collapsed/feature-table-${n}.tsv \
--to-tsv; done


# Filter tables based on SampleSite combinations (3, 2 and 1) and determine core features

mkdir qiime2_output/16S_v4_merged/filtered_tables

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] IN ('CLOACA', 'ORAL', 'CARAPACE')" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-orl-car.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] IN ('CARAPACE', 'ORAL', 'TANK WATER')" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-car-orl-wtr.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] IN ('CLOACA', 'ORAL', 'TANK WATER')" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-orl-wtr.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] IN ('CLOACA', 'CARAPACE', 'TANK WATER')" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-car-wtr.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] = 'CLOACA'" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-cloacal.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] = 'ORAL'" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-oral.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] = 'CARAPACE'" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-carapace.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] = 'TANK WATER'" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-water.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] IN ('CLOACA', 'ORAL')" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-orl.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] IN ('CLOACA', 'CARAPACE')" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-car.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] IN ('CLOACA', 'TANK WATER')" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-wtr.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] IN ('ORAL', 'CARAPACE')" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-orl-car.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] IN ('ORAL', 'TANK WATER')" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-orl-wtr.qza

qiime feature-table filter-samples \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --m-metadata-file master_metadata/16S_all_metadata.tsv \
  --p-where "[SampleSite] IN ('CARAPACE', 'TANK WATER')" \
  --o-filtered-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-car-wtr.qza

qiime feature-table core-features \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --o-visualization qiime2_output/16S_v4_merged/core_features/table-16S_v4_merged-no-mitochleuk-D29000-core_features.qzv

for n in {1..7};\
do qiime feature-table core-features \
  --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}.qza \
  --o-visualization qiime2_output/16S_v4_merged/core_features/table-${n}-core_features.qzv; done

cd qiime2_output/16S_v4_merged/filtered_tables

for file in *.qza;\
do qiime feature-table core-features \
  --i-table ${file} \
  --o-visualization ${file}-core_asv.qzv; done

for text in *.qza-core_asv.qzv; do mv "$text" $(echo "$text" | sed 's/.qza//g'); done

cd ../../..

mv qiime2_output/16S_v4_merged/filtered_tables/*.qzv qiime2_output/16S_v4_merged/core_features

# Core features in SampleSite combinations but on collapsed taxonomy - only family, genus and species

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-car-orl-wtr.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-car-orl-wtr.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-car-wtr.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-car-wtr.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-carapace.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-carapace.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-oral.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-oral.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-water.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-water.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-car-wtr.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-car-wtr.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-car.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-car.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-orl-car.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-orl-car.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-orl.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-orl.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-wtr.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-wtr.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-orl-car.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-orl-car.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-orl-wtr.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-orl-wtr.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-clo-orl-wtr.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-orl-wtr.qza; done

for n in {5..7};\
do qiime taxa collapse \
  --i-table qiime2_output/16S_v4_merged/filtered_tables/16S_v4_merged-table-cloacal.qza \
  --i-taxonomy qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --p-level ${n} \
  --o-collapsed-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-cloacal.qza; done

# core on filtered and collapsed to 5, 6, 7 lvl

for n in {5..7};\
do qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-cloacal.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-cloacal-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-oral.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-oral-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-water.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-water-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-carapace.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-carapace-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-car-orl-wtr.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-car-orl-wtr-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-car-wtr.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-car-wtr-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-car-wtr.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-clo-car-wtr-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-car.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-clo-car-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-orl-car.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-clo-orl-car-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-orl-wtr.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-clo-orl-wtr-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-orl.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-clo-orl-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-clo-wtr.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-clo-wtr-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-orl-car.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-orl-car-core_features.qzv &&
qiime feature-table core-features \
   --i-table qiime2_output/16S_v4_merged/taxonomy_collapsed/table-${n}-orl-wtr.qza \
   --o-visualization qiime2_output/16S_v4_merged/core_features/16S_v4_merged-table-${n}-orl-wtr-core_features.qzv; done


#----- DEICODE on table at 29000 depth and QURRO on rPCA loadings

qiime deicode rpca \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 29000 \
  --o-biplot qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination-29000a.qza \
  --o-distance-matrix qiime2_output/16S_v4_merged/diversity/deicode/rPCA-matrix-29000a.qza

qiime emperor biplot \
  --i-biplot qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination-29000a.qza \
  --m-sample-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-feature-metadata-file qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --o-visualization qiime2_output/16S_v4_merged/diversity/deicode/rPCA-biplot-29000a.qzv \
  --p-number-of-features 8

qiime qurro loading-plot \
  --i-table qiime2_output/16S_v4_merged/table-16S_v4_merged-no-mitochleuk-D29000.qza \
  --i-ranks qiime2_output/16S_v4_merged/diversity/deicode/rPCA-ordination-29000a.qza \
  --m-sample-metadata-file master_metadata/16S_all_metadata.tsv \
  --m-feature-metadata-file qiime2_output/16S_v4_merged/taxonomy-16S_v4_merged.qza \
  --o-visualization qiime2_output/16S_v4_merged/diversity/deicode/qurro-plot-29000a.qzv