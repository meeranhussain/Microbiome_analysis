# Qiime2 Pipeline for 16s rRNA and ITS Sequence Analysis
**This pipeline was used in the following research work:** [Regenerative agriculture augments bacterial community structure for a healthier soil and agriculture](https://doi.org/10.3389/fagro.2023.1134514)

The below pipeline provides codes for processing 16s rRNA data using Qiime2, from initial sequence trimming to taxonomy classification and visualization. Same steps can also be followed for processing ITS data

## 1. Load Sample List & Perform Quality Trimming
First, we load the sample IDs from a text file.

```bash
#!/bin/bash

### Load sample list
input="Sample_id.txt"

while IFS= read -r i
do
  echo $i
  trim_galore -q 20 --paired --fastqc --cores 4 ${i}_R1.fq.gz ${i}_R2.fq.gz -o trimd_files
done < "$input"
```
## 2. Prepare Sample Metadata Sheet
Prepare a metadata sheet based on the study design and save it as soil_analysis_metadata.tsv.

`EXAMPLE SHEET`
| sample_id | Farming_Method | Location   | Kind_of_Vegetation |
|-----------|----------------|------------|--------------------|
| Indira-1  | Regenerative   | Ramnagara  | Beans              |
| Indira-2  | Regenerative   | Ramnagara  | Ragi               |
| Indira-3  | Regenerative   | Magadi     | Tomato             |
| Indira-4  | Conventional   | Ramnagara  | Beans              |
| Indira-5  | Regenerative   | Magadi     | Beans              |
| Indira-6  | Conventional   | Ramnagara  | Beans              |
| Indira-7  | Regenerative   | Ramanagara | Ragi               |
| Indira-8  | Regenerative   | Magadi     | Ragi               |
| Indira-9  | Regenerative   | Hosur      | Tomato/Beans       |

### Visualize the metadata using Qiime2:

```bash
qiime metadata tabulate \
  --m-input-file soil_analysis_metadata.tsv \
  --o-visualization soil_analysis_metadata.qzv
```
## 3. Prepare Sample Load Sheet
Create a load sheet for Qiime2 named soil_samples.tsv:

`EXAMPLE SHEET`
| SampleID | Forward_read_path                   | Reverse_read_path                   |
|----------|-------------------------------------|-------------------------------------|
| Indira-1 | trimd_files/Indira_1_R1.fq.gz       | trimd_files/Indira_1_R2.fq.gz       |
| Indira-2 | trimd_files/Indira_2_R1.fq.gz       | trimd_files/Indira_2_R2.fq.gz       |
| Indira-3 | trimd_files/Indira_3_R1.fq.gz       | trimd_files/Indira_3_R2.fq.gz       |
| Indira-4 | trimd_files/Indira_4_R1.fq.gz       | trimd_files/Indira_4_R2.fq.gz       |
| Indira-5 | trimd_files/Indira_5_R1.fq.gz       | trimd_files/Indira_5_R2.fq.gz       |
| Indira-6 | trimd_files/Indira_6_R1.fq.gz       | trimd_files/Indira_6_R2.fq.gz       |
| Indira-7 | trimd_files/Indira_7_R1.fq.gz       | trimd_files/Indira_7_R2.fq.gz       |
| Indira-8 | trimd_files/Indira_8_R1.fq.gz       | trimd_files/Indira_8_R2.fq.gz       |
| Indira-9 | trimd_files/Indira_9_R1.fq.gz       | trimd_files/Indira_9_R2.fq.gz       |

### Import the sequence data into Qiime2 to create a SampleData[PairedEndSequencesWithQuality] artifact:
```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path soil_samples.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
```
## 4. Visualize Demultiplexed Data
Visualize the paired-end demultiplexed sample data in Qiime2:
```bash
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization paired-end-demux.qzv
```
## 5. Denoise Data with DADA2
Denoise the samples using DADA2 and generate a visualization file:

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trunc-len-f 285 \
  --p-trunc-len-r 250 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table pet-table.qza \
  --o-denoising-stats stats-dada2.qza \
  --p-n-threads 8

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

mv rep-seqs-dada2.qza rep-seqs.qza
mv pet-table.qza table.qza
```

## 6. Create Feature Table and Feature Data Summaries
Summarize the feature table and feature data:
```bash
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file soil_analysis_metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
```
## 7. Taxonomy Classification
Classify the sequences using a pre-trained classifier and visualize the taxonomy:

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file soil_analysis_metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```
## 8. Collapse Taxa and Export Data
Collapse the taxonomy at level 7 and export the data:

```bash
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table collapsed-l2.qza

qiime tools export \
  --input-path collapsed-l2.qza \
  --output-path collapsed-l2-dir
```

