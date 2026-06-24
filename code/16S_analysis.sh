#!/bin/bash

cd data/16S

# Import reads
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path import_to_qiime \
  --output-path reads

# Analyze quality scores
qiime demux summarize \
  --p-n 10000 \
  --i-data reads.qza \
  --o-visualization qual_viz

# Denoise reads
mkdir feature_tables

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs reads.qza \
  --o-table feature_tables/unfiltered_table \
  --o-representative-sequences representative_sequences \
  --p-trunc-len-r 200 \
  --p-trunc-len-f 260 \
  --p-trim-left-f 12 \
  --p-trim-left-r 12 \
  --p-n-threads 4 \
  --o-denoising-stats denoise_stats.qza \
  --verbose

  qiime feature-table summarize \
    --i-table feature_tables/unfiltered_table.qza \
    --o-visualization unfiltered_table.qzv

# Taxonomic classification
qiime feature-classifier classify-sklearn \
  --i-classifier /home/rbeiko/links/projects/rrg-rbeiko/shared/Proj-GORP/Scripts/classifiers/16S-18S-silva-138-99-nb-classifier.qza \
  --i-reads representative_sequences.qza \
  --o-classification taxonomy \
  --p-reads-per-batch 1000 \
  --verbose

qiime tools export \
    --input-path taxonomy.qza \
    --output-path .

# Generate phylogeny
mkdir phylogeny

qiime alignment mafft \
  --i-sequences representative_sequences.qza \
  --o-alignment phylogeny/aligned_representative_sequences

qiime alignment mask \
  --i-alignment phylogeny/aligned_representative_sequences.qza \
  --o-masked-alignment phylogeny/masked_aligned_representative_sequences

qiime phylogeny fasttree \
  --i-alignment phylogeny/masked_aligned_representative_sequences.qza \
  --o-tree phylogeny/unrooted_tree

qiime phylogeny midpoint-root \
  --i-tree phylogeny/unrooted_tree.qza \
  --o-rooted-tree phylogeny/rooted_tree

# Filter out low abundance (i.e. possible bleed-through) ASVs
qiime feature-table filter-features \
  --i-table feature_tables/unfiltered_table.qza \
  --p-min-frequency 10 `# 0.1% of mean sample depth` \
  --p-min-samples 1 \
  --o-filtered-table feature_tables/table_bleed.qza

qiime feature-table summarize \
  --i-table feature_tables/table_bleed.qza \
  --o-visualization feature_tables/table_bleed.qzv

# Filter out contaminant (mitochondria, chloroplast and unclassified) ASVs
qiime taxa filter-table \
  --i-table feature_tables/table_bleed.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --p-include p__ \
  --o-filtered-table feature_tables/table_bleed_decon.qza

qiime feature-table summarize \
  --i-table feature_tables/table_bleed_decon.qza \
  --o-visualization feature_tables/table_bleed_decon.qzv

qiime tools export \
    --input-path feature_tables/table_bleed_decon.qza \
    --output-path feature_tables

biom convert \
    -i feature_tables/feature-table.biom \
    -o feature_tables/table_bleed_decon.tsv \
    --to-tsv

# Remove one low-depth sample
qiime feature-table filter-samples \
  --i-table feature_tables/table_bleed_decon.qza \
  --o-filtered-table feature_tables/table_filtered.qza \
  --p-min-frequency 3000

## Additional code by MAF
qiime feature-table filter-samples --i-table table_bleed_decon.qza --o-filtered-table table_filtered.qza --p-min-frequency 2500

# Calculate alpha diversity with rarefaction
for metric in "observed_features" "shannon" "pielou_e" "faith_pd"; do
  qiime boots alpha \
    --i-table feature_tables/table_filtered.qza \
    --i-phylogeny phylogeny/rooted_tree.qza \
    --p-sampling-depth 3000 \
    --p-metric ${metric} \
    --p-n 100 \
    --p-no-replacement \
    --p-average-method median \
    --o-average-alpha-diversity alpha_diversity/${metric}.qza; done

# Export alpha diversity results
for infile in alpha_diversity/*; do
  base="$(basename "$infile" .qza)"
  qiime tools export \
    --input-path ${infile} \
    --output-path alpha_diversity
  mv alpha_diversity/alpha-diversity.tsv alpha_diversity/${base}.tsv
  rm alpha_diversity/alpha-diversity.tsv; done

# Calculate Aitchison distance and Weighted UniFrac with rarefaction
for metric in "aitchison" "weighted_unifrac"; do
qiime boots beta \
  --i-table feature_tables/table_filtered.qza \
  --i-phylogeny phylogeny/rooted_tree.qza \
  --p-sampling-depth 3000 \
  --p-metric aitchison \
  --p-n 100 \
  --p-no-replacement \
  --p-average-method medoid \
  --o-average-distance-matrix beta_diversity/aitchison_distance.qza
# Export distance matrices
qiime tools export \
  --input-path beta_diversity/aitchison_distance.qza \
  --output-path beta_diversity
mv beta_diversity/distance-matrix.tsv beta_diversity/aitchison_distance.tsv

# Calculate phylogenetic RPCA
qiime gemelli phylogenetic-rpca-with-taxonomy \
  --i-table feature_tables/table_filtered.qza \
  --i-phylogeny phylogeny/rooted_tree.qza  \
  --m-taxonomy-file taxonomy.qza \
  --o-biplot beta_diversity/rpca_biplot.qza \
  --o-distance-matrix beta_diversity/rpca_distance.qza \
  --o-counts-by-node-tree beta_diversity//rpca_phylo-tree.qza \
  --o-counts-by-node beta_diversity/rpca_phylo-table.qza \
  --o-t2t-taxonomy beta_diversity/rpca_phylo-taxonomy.qza
# Export distance matrix and ordination results for phylogenetic RPCA
qiime tools export \
  --input-path beta_diversity/rpca_distance.qza \
  --output-path beta_diversity
mv beta_diversity/distance-matrix.tsv beta_diversity/rpca_distance.tsv
qiime tools export \
  --input-path beta_diversity/rpca_biplot.qza \
  --output-path beta_diversity
mv beta_diversity/ordination.txt beta_diversity/rpca_pcoa.tsv

# Pool samples per subject for PERMANOVA
qiime feature-table group \
  --i-table feature_tables/table_filtered.qza \
  --p-axis 'sample' \
  --m-metadata-file ../metadata/sample_history.tsv \
  --m-metadata-column 'SubjectID' \
  --p-mode 'mean-ceiling' \
  --o-grouped-table feature_tables/table_filtered_grouped.qza

# Calculate Aitchison distance and PCoA for pooled data
for metric in "aitchison" "weighted_unifrac"; do
qiime boots beta \
  --i-table feature_tables/table_filtered_grouped.qza \
  --i-phylogeny phylogeny/rooted_tree.qza \
  --p-sampling-depth 3000 \
  --p-metric aitchison \
  --p-n 100 \
  --p-no-replacement \
  --p-average-method medoid \
  --o-average-distance-matrix beta_diversity_grouped/aitchison_distance.qza; done
# Export distance matrices
qiime tools export \
  --input-path beta_diversity_grouped/aitchison_distance.qza \
  --output-path beta_diversity_grouped
mv beta_diversity_grouped/distance-matrix.tsv beta_diversity_grouped/aitchison_distance.tsv
# Calculate PCoA for Aitchison distance and weighted UniFrac
qiime diversity pcoa \
  --i-distance-matrix beta_diversity_grouped/aitchison_distance.qza \
  --o-pcoa beta_diversity_grouped/aitchison_pcoa.qza
# Export ordination results
qiime tools export \
  --input-path beta_diversity_grouped/aitchison_pcoa.qza \
  --output-path beta_diversity_grouped
mv beta_diversity_grouped/ordination.txt beta_diversity_grouped/aitchison_pcoa.tsv

# Phylogenetic RPCA
qiime gemelli phylogenetic-rpca-with-taxonomy \
  --i-table feature_tables/table_filtered_grouped.qza \
  --i-phylogeny phylogeny/rooted_tree.qza  \
  --m-taxonomy-file taxonomy.qza \
  --o-biplot beta_diversity_grouped/rpca_biplot.qza \
  --o-distance-matrix beta_diversity_grouped/rpca_distance.qza \
  --o-counts-by-node-tree beta_diversity_grouped/rpca_phylo-tree.qza \
  --o-counts-by-node beta_diversity_grouped/rpca_phylo-table.qza \
  --o-t2t-taxonomy beta_diversity_grouped/rpca_phylo-taxonomy.qza
qiime tools export \
  --input-path beta_diversity_grouped/rpca_distance.qza \
  --output-path beta_diversity_grouped
mv beta_diversity_grouped/distance-matrix.tsv beta_diversity_grouped/rpca_distance.tsv
qiime tools export \
  --input-path beta_diversity_grouped/rpca_biplot.qza \
  --output-path beta_diversity_grouped
mv beta_diversity_grouped/ordination.txt beta_diversity_grouped/rpca_pcoa.tsv
