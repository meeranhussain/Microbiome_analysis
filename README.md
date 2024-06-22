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
## 9. Custom taxanomy plotting using CSV file downloaded from "taxa-bar-plots.qzv"
```R
# Install necessary packages
install.packages("dplyr")
install.packages("tibble")

# Load libraries
library(dplyr)
library(tibble)

# Read the CSV file
tax_lvl6 <- read.csv("Taxa_level-6.csv", header = TRUE, sep = ",", row.names = NULL)

# Extract sample IDs
sam_ID <- select(tax_lvl6, c(1))

# Remove specified columns
tax_lvl6 <- select(tax_lvl6, -c(1, 1056:1075))

# Add total species count in each column
tax_lvl6 <- rbind(tax_lvl6, mapply(sum, tax_lvl6[, c(1:1054)]))

# Order and select top 100 rows
list <- order(tax_lvl6[c(15), c(1:1054)], decreasing = TRUE)
top_100 <- tax_lvl6[list]
top_100 <- top_100[-c(15),]

# Calculate percentages
tax_lvl6_pct <- t(apply(top_100, 1, function(x) { x / sum(x) * 100 }))[, 1:100]

# Create a data frame with top 100 species
top_100 <- data.frame(tax_lvl6_pct[, c(1:100)])

# Merge sample ID column
top_100 <- add_column(top_100, sam_ID, .before = 1)
colnames(top_100)[1] <- c("sample.id")

# Selecting only ragi cultivated plants
ragi_100 <- top_100[-c(1, 2, 8, 9, 10, 11, 14),]

# Read metadata table
meta_dta <- read.table("soil_analysis_metadata.tsv", header = TRUE, sep = "\t")

# Remove specified columns
col_to_be_removed <- colnames(ragi_100)[c(4, 22, 34, 45, 49, 54, 63, 64, 81, 86, 88, 93, 98)]
ragi_sel <- ragi_100[, -c(4, 22, 34, 45, 49, 54, 63, 64, 81, 86, 88, 93, 98)]

# Merge data frames
ragi_sel_meta <- merge(ragi_sel, meta_dta)

# Melt data frame
melt_df_2 <- melt(ragi_sel_meta[, -c(1, 89, 90, 92:107)], id.vars = "Farming.Method")

# For 50 species
ragi_50_meta <- ragi_sel_meta[, -c(1, 2, 51:90, 92:107)]
melt_50 <- melt(ragi_50_meta, id.vars = "Farming.Method")

# Rename columns
colnames(melt_50) <- c("Farming.Method", "Species", "percentage")

# Plot using ggplot
ggplot(melt_50, aes(fill = Farming.Method, y = percentage, x = Species)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text.x = element_text(angle = 90))
```
## 9. Qiime Filteration & Alpha Diversity
**Filtering Samples:**
```bash
qiime feature-table filter-samples \
     --i-table ../table.qza \
     --m-metadata-file ../soil_analysis_metadata.tsv \
     --p-where "[Crop_Type] IN ('Ragi', 'Barren')" \
     --o-filtered-table ragi-filtered-table.qza
```
**Summarizing Filtered Table:**
```bash
qiime feature-table summarize \
  --i-table ragi-filtered-table.qza \
  --o-visualization ragi-filtered-table.qzv \
  --m-sample-metadata-file ../soil_analysis_metadata.tsv
```

**Alpha Rarefaction Analysis:**
```bash
qiime diversity alpha-rarefaction \
  --i-table ragi-filtered-table.qza \
  --i-phylogeny ../rooted-tree.qza \
  --p-max-depth 24174 \
  --m-metadata-file ../soil_analysis_metadata.tsv \
  --o-visualization ragi-alpha-rarefaction.qzv
```
**Core Metrics Phylogenetic Analysis:**
```bash
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ../rooted-tree.qza \
  --i-table ragi-filtered-table.qza \
  --p-sampling-depth 14222 \
  --m-metadata-file ../soil_analysis_metadata.tsv \
  --output-dir core-metrics-results
```
**Alpha Group Significance - Faith's PD:**
```bash
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ../soil_analysis_metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv
```

**Alpha Group Significance - Evenness:**

```bash
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file ../soil_analysis_metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv
```

## 10. Custom ploting of Taxonomy, Alpha Rarefaction, Alpha Diversity: (Needs a bit clean) 

```R
install.packages("reshape2")
library(reshape)
library(dplyr)
library(tibble)
install.packages("tidyverse")
library(tidyverse)
install.packages("hrbrthemes")
library(hrbrthemes)
library(viridis)

tax_lvl6 <- read.csv("level-6.csv", header = TRUE, sep = ",", row.names = NULL)
#tax_lvl6 <- tax_lvl6[-c(7),]
sam_ID <- select(tax_lvl6, c(1)) #extract IDs

tax_lvl6 <- select(tax_lvl6,-c(1,750:759)) #removed col1 & 47 to 56
#rownames(tax_lvl6) <- sam_ID[,1]

#tax_lvl6_pct <- t(apply(tax_lvl6[order(rowSums(tax_lvl6), decreasing = T),order(colSums(tax_lvl6), decreasing =T)], 1, function(x) { x / sum(x) * 100}))[,1:45]

tax_lvl6 <- rbind(tax_lvl6, mapply(sum,tax_lvl6[,c(1:748)])) #added the total species count in each col

list <- order(tax_lvl6[c(8),c(1:748)], decreasing = TRUE ) # ordered in increasing - has the index values

top_100 <- tax_lvl6[list] # sorted based on the index values 

top_100 <- top_100[-c(8),] #removed the last row

tax_lvl6_pct <- t(apply(top_100, 1, function(x) { x / sum(x) * 100}))[,1:45]

top_100 <- data.frame(tax_lvl6_pct[,c(1:45)])


top_100 <- add_column(top_100, sam_ID, .before = 1) # merged the sample id col

colnames(top_100)[1] <- c("sample.id") #added 1st col name



#Selecting only raji cultivated plants
#ragi_100 <- top_100[-c(1,2,8,9,10,11,14), ] #selected only ragi plots


meta_dta <- read.table("soil_analysis_metadata.tsv", header = TRUE, sep = "\t") #read the meta table


colnames(top_100) <- gsub(".*p__","",colnames(top_100)) #remove all the strings before ".g__"

#col_to_be_removed <- colnames(ragi_100)[c(4,22,34,45,49,54,63,64,81,86,88,93,98)] #remove the col with unidentified genus

#ragi_sel <- ragi_100[ ,-c(4,22,34,45,49,54,63,64,81,86,88,93,98)]

#ragi_sel_meta <- merge(ragi_sel, meta_dta)

rag_sel_meta <- merge(top_100, meta_dta)



melt_df_2 <- melt(rag_sel_meta[ ,-c(1, 14:49, 51:56 )], id.vars="Farming.Type") #converted to long table

#For 50 species
#ragi_50_meta <- ragi_sel_meta[ ,-c(1,2,51:90,92:107)]
#melt_50 <- melt(ragi_50_meta, id.vars="Farming.Method") #uses reshape2, removed col other than species


colnames(melt_df_2) <- c("Farming.Type","Phylum","Percentage")


#ggplot

melt_df_2$Farming.Type <- factor(melt_df_2$Farming.Type,                                    # Change ordering manually
                                   levels = c("Barren", "Conv", "Reg (1y)", "Reg (>5y)"))
ggplot(melt_df_2, aes(fill=Farming.Type, y=Percentage, x=Phylum)) +  geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(labels = scales::percent)
#plot-2

ggplot(melt_df_2, aes(fill=Phylum, y=Percentage, x=Farming.Type)) +  geom_bar(position="fill", stat="identity") + 
  scale_fill_brewer(palette="Paired") + theme(axis.text.x = element_text(angle = 90))  + scale_y_continuous(labels = scales::percent)

ggplot(melt_df_2, aes(fill=Phylum, y=Percentage, x=Farming.Type)) +  geom_bar(position="stack", stat="summary", fun ="mean" ) + scale_fill_brewer(palette="Paired") + theme(axis.text.x = element_text(angle = 90))
+ scale_y_continuous(labels = scales::percent)


#plot-3 for the faith-PD
ggplot(a, aes(x=as.factor(Farming.Type), y=faith_pd)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("Farming.Type") 

#Alpha_diversity plots
a <- read.csv("Alpha_diversity_metadata.tsv", header = TRUE, sep = "\t", row.names = NULL)

# Jitter plot

ggplot(a, aes(x = Farming.Type, y = faith_pd)) + geom_jitter(position = position_jitter(0.2))

#Pie Chart
#type-1

pie(sapply(split(a$faith_pd, a$Farming.Type), mean))
#type-2
ggplot(a, aes(x="", y=faith_pd, fill=Farming.Type)) +
  geom_bar(stat="summary", fun="mean",  width=1) +
  coord_polar("y", start=0)
#scattered plot
ggplot(a, aes(x=Farming.Type, y=faith_pd, color=Farming.Type)) + 
  geom_point(size=6) +
  theme_ipsum()

#Eveness
evn_data <- read.table("evenness_ragi.tsv", header = T, sep = "\t", row.names = NULL)

# Jitter plot
ggplot(evn_data, aes(x = Farming.Type, y = faith_pd)) + geom_jitter(position = position_jitter(0.2))

#scattered plot
colnames(evn_data)[3] <"Farming Type"
evn_data$Farming.Type <- factor(evn_data$Farming.Type,                                    # Change ordering manually
                                   levels = c("Barren", "Conv", "Reg (1y)", "Reg (>5y)"))
p2 <- ggplot(evn_data, aes(x=Farming.Type, y=`Bacterial Population Evenness`, color=Farming.Type)) + 
  geom_point(size=6) + theme(axis.title = element_text(size = 7))

ggsave("Rplot_evenness(p=0.14).TIFF", plot=p2, height=8, width=14, units=c("cm"), device='tiff', dpi=300)


ggplot(evn_data, aes(x = Farming.Type, y = faith_pd)) + geom_jitter(position = position_jitter(0.2)) + stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.5)


data <- read.table("alpha_rac.tsv", header=T, sep = "\t")
melt_df_alpha <- melt(data[ ,-c(1, 62:101,103:110)], id.vars="Samples") #converted to long table
melt_df_alpha <- melt_df_alpha[order(melt_df_alpha$value),]
# Plot
color <- c("darkolivegreen4", "black", "blue", "darkmagenta", "chocolate1", "darkgoldenrod4", "brown1", "darkviolet", "deeppink")
colnames(data)[2:3] <- c("Number of Reads", "Soil Bacterial Diversity")
p5 <- ggplot(data, aes(x=`Number of Reads`, y=`Soil Bacterial Diversity`, color=Samples)) + geom_line( stat="identity") + geom_point(stat = "identity", size = 1 ) + scale_color_manual(values = color) + theme(axis.text.x = element_text(angle = 90)) + scale_x_comma()


ggsave("Alpha_racfac_ragi.tiff", plot=p5, height=8, width=14, units=c("cm"),device='tiff', dpi=600)


#Selected genus for taxanomy plot

colnames(top_100) <- gsub(".*g__","",colnames(top_100)) #remove all the strings before ".g__"
top_rm_100 <- data.frame(top_100)

sel_genus <- data.frame(top_rm_100$sample.id, top_rm_100$Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium, top_rm_100$Pseudomonas, top_rm_100$Bacillus, top_rm_100$Nocardioides, top_rm_100$Streptomyces, top_rm_100$Nocardia, top_rm_100$Mycobacterium, top_rm_100$Micromonospora )
colnames(sel_genus) <- c("sample.id", "Rhizobium", "Pseudomonas", "Bacillus", "Nocardioides","Streptomyces", "Nocardia", "Mycobacterium", "Micromonospora")
library(tibble)
sel_gen_lng <- merge(sel_genus, meta_dta, by="sample.id")
sel_gen_lng <- melt(sel_gen_lng[ ,-c(1, 10:12, 14:19 )], id.vars="Farming.Type")
colnames(sel_gen_lng) <- c("Farming.Type","Genus","Counts")  
ggplot(sel_gen_lng, aes(fill=Genus, y=Counts, x=Farming.Type)) +  geom_bar(position="stack", stat="summary", fun ="mean")  + scale_fill_brewer(palette="Paired") + theme(axis.text.x = element_text(angle = 90)) + coord_flip()


#Selected taxa
sel_all <- data.frame( top_rm_100$sample.id, top_rm_100$Flavobacterium , top_rm_100$Bacillus , top_rm_100$Streptomyces , top_rm_100$Mesorhizobium , top_rm_100$Achromobacter, top_rm_100$Klebsiella , top_rm_100$Paenibacillus , top_rm_100$Burkholderia.Caballeronia.Paraburkholderia , top_rm_100$Pseudomonas)

colnames(sel_all) <- c("sample.id","Flavobacterium","Bacillus","Streptomyces","Mesorhizobium","Achromobacter","Klebsiella","Paenibacillus","Burkholderia","Pseudomonas")

sel_all_lng <- merge(sel_all, meta_dta, by="sample.id")

sel_all_lng <- melt(sel_all_lng[ ,-c(1, 11:13, 15:20 )], id.vars="Farming.Type")

colnames(sel_all_lng) <- c("Farming_Type","Genus","Percentage") 

sel_all_lng$Farming_Type <- factor(sel_all_lng$Farming_Type,                                    # Change ordering manually
       levels = c("Barren", "Conv", "Reg (1y)", "Reg (>5y)"))

p3 <- ggplot(sel_all_lng, aes(fill=Genus, y=Percentage, x=Farming_Type)) +  geom_bar(position="fill", stat="summary", fun ="mean" )  + scale_fill_brewer(palette="Paired") + theme(axis.text.x = element_text(angle = 90)) + scale_y_continuous(labels = scales::percent) + coord_flip() 
ggsave("sel_taxa_1.TIFF", plot=p3, height=8, width=14, units=c("cm"), device='tiff', dpi=300)


###Richness
richness <- read.table("faith_pd.csv", header = T, sep = ",", row.names = NULL)

#scattered plot
colnames(a)[10] <- "Soil Bacterial Diversity"
a$Farming.Type <- factor(a$Farming.Type,                                    # Change ordering manually
                                levels = c("Barren", "Conv", "Reg (1y)", "Reg (>5y)"))
p2 <- ggplot(a, aes(x=Farming.Type, y=`Soil Bacterial Diversity`, color=Farming.Type)) + 
  geom_point(size=6) + theme(axis.title = element_text(size = 7))

ggsave("Rplot_richness(p=0.169).TIFF", plot=p2, height=8, width=14, units=c("cm"), device='tiff', dpi=300)
```

