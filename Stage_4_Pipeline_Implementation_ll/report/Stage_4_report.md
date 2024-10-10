# Differential Gene Expersion Analysis and Machine Learning of Gliomas

## Authors:
Niguse Kelile Lema (@King),
Ala Khaled (@Alaa043),
Ahmed Abdel-Masqoud (@Abdel-Maqsoud),
Obianujunwa Martha (@obianujunwa),
Raqeebah Rafiu (@Raqeebahh),
Abilashni Arthiswaran (@Abilashni83),
Benedict Orile Ajunku (@orile)

## Introduction
Gliomas are a type of brain tumor originating from glial cells, which support nerve cells in the central nervous system. They are classified by their aggressiveness, with glioblastomas being the most severe. The isocitrate dehydrogenase (IDH) gene mutation is a key factor in glioma classification, as IDH-mutant gliomas have a better prognosis than IDH-wildtype. IDH status plays a critical role in understanding tumor behavior, guiding treatment strategies, and assessing patient survival outcomes in neuro-oncology. The most commonly occurring malignant brain tumors are gliomas (Ceccarelli et al., 2016). 


## Methodology and Results

### Data Acquastion and Preprocessing Steps
Once the data is downloaded, it is prepared for analysis using GDCprepare, which organizes the data into usable formats. 
The raw counts for both LGG and GBM are extracted and stored for further processing.
The next phase involves preparing the metadata associated with the datasets. 
The project collects barcode identifiers and the IDH status (Isocitrate Dehydrogenase status) for each sample, which are essential for downstream analysis. 
Missing values in the metadata are omitted to ensure data integrity.
The final datasets are prepared, consisting of 513 samples for LGG and 154 samples for GBM. These datasets are then merged into a single data frame, facilitating a comprehensive analysis of both cancer types.

### Data Normalization and Filtering
The combined dataset undergoes normalization to adjust for technical variations. 
The normalization methods employed include adjustments for gene length and GC content, which are critical for ensuring that gene expression levels are accurately represented. 
Subsequently, a quantile filtering method is applied to remove lowly expressed genes, retaining only the most informative features for analysis.

### Differential Expression Analysis
Differential expression analysis (DEA) is performed to compare gene expression levels between wild-type (WT) and mutant IDH groups. Using the edgeR package, the project identifies genes that exhibit statistically significant changes in expression between these two groups. 
The results of the DEA are summarized in a table, providing insights into which genes are upregulated or downregulated.

### Data Annotation and Visualization
To enhance the interpretability of the results, the project annotates the differentially expressed genes (DEGs) using the biomaRt package. 
This step associates Ensembl gene IDs with their corresponding HGNC symbols, facilitating biological interpretation.

![Logo](https://example.com/logo.png)
*Figure 1: The volcano plot illustrates the relationship between fold change and false discovery rate (FDR) for the identified DEGs.*

![Logo](https://example.com/logo.png)
*Figure 2: The heatmap visualizes expression patterns of DEGs, employing color coding to differentiate between WT and mutant samples.*

### Machine Learning Analysis
To delve deeper into the data, machine learning techniques are applied. 
The project utilizes the caret package to preprocess the gene expression data, focusing on variability and correlation among genes. 
The top 3,000 genes with the highest variability are selected for clustering analysis.
-Means clustering is performed to identify potential subgroups within the glioma data, with the results visualized using fviz_cluster. 
This clustering analysis helps uncover underlying patterns in gene expression and assess the impact of IDH status on gene expression profiles.

## Conclusion
The project successfully integrates cancer genomics with machine learning techniques to advance the understanding of gliomas. 

## Rscript:
## Report:

## References
1. Ceccarelli, M., Barthel, F. P., Malta, T. M., Sabedot, T. S., Salama, S. R., Murray, B. A., Morozova, O., Newton, Y., Radenbaugh, A., Pagnotta, S. M., Anjum, S., Wang, J., Manyam, G., Zoppoli, P., Ling, S., Rao, A. A., Grifford, M., Cherniack, A. D., Zhang, H., … Verhaak, R. G. W. (2016). Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma. Cell, 164(3), 550–563. https://doi.org/10.1016/j.cell.2015.12.028
2. Neil, Z. D., Pierzchajlo, N., Boyett, C., Little, O., Kuo, C. C., Brown, N. J., & Gendreau, J. (2023). Assessing Metabolic Markers in Glioblastoma Using Machine Learning: A Systematic Review. Metabolites, 13(2). https://doi.org/10.3390/metabo13020161
3. Shields, L. B. E., & Choucair, A. K. (2014). Management of Low-Grade Gliomas: A Review of Patient-Perceived Quality of Life and Neurocognitive Outcome. World Neurosurgery, 82(1–2), e299–e309. https://doi.org/10.1016/j.wneu.2014.02.033
4. Vrazhnov, D., Mankova, A., Stupak, E., Kistenev, Y., Shkurinov, A., & Cherkasova, O. (2023). Discovering Glioma Tissue through Its Biomarkers’ Detection in Blood by Raman Spectroscopy and Machine Learning. Pharmaceutics, 15(1). https://doi.org/10.3390/pharmaceutics15010203

