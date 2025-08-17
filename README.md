# PSC_IBD-Omics-Preprcossing
this file include the preprocssing steps for RNA-Seq, Metatranscriptomics, 16S rRNA and Bile Acid

RNA-Seq
Preprocessing involved the identification and removal of missing values; 16 genes that were absent across all samples were excluded from further analysis. The remaining 1088 data were then log-transformed using the log1p function, which adds 1 to each value prior to taking the natural logarithm. This transformation helps address the highly skewed and zero-inflated nature of RNA-Seq data. Subsequently, principal component analysis (PCA) was performed to explore the overall data structure and visualize patterns between baseline and week 4 samples.

Metatranscriptomics 
Data preprocessing involved transformation and normalization steps. Metatranscriptomic data were available across all time points, enabling temporal exploration of transcriptional shifts via PCA; however, only baseline and week 4 data were used for downstream predictive modeling.

Bile Acid Metabolites
Subsequent preprocessing included log1p transformation and z-score normalization. PCA was then performed to explore potential patterns or shifts in bile acid profiles across different treatment time points.

16s rRNA 
Due to the compositional and sparse nature of the data, preprocessing steps were adapted from (Hodgkiss and Acharjee, 2025). Specifically, only genus-level taxa were retained, as they were consistently detected across the majority of samples. Low-abundance genera were filtered out, and a pseudocount was added to mitigate the impact of zero values prior to applying the centered log-ratio (CLR) transformation, which is appropriate for compositional data. Following CLR transformation and scaling, Bray–Curtis dissimilarity was calculated, and the data were visualized using Principal Coordinates Analysis (PCoA). After these preprocessing steps out of 177 taxa only 63 remained for downstream analysis.





