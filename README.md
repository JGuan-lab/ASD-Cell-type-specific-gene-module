# Cell type-specific gene network-based analysis

Cell type-specific gene network-based analysis can be applied to analyze single cell/nucleus gene expression data for identifying cell type-specific gene modules. Further it can be used to identify cell type-specific disease-associated gene modules.

Specifically, for each cell type, a gene co-expression network is constructed. Then disease genes were used to identify disease-associated gene modules. Next, cell type specificity of each gene is calculated, and cell type-specific genes were used to assess the cell type enrichment for each disease-associated module. The disease-associated gene modules built from a cell type significantly enriched with the specific genes of the cell type but not enriched with other kinds of cell type-specific genes are considered as candidate cell type-specific disease-associated gene modules. Then the function of modulePreservation in WGCNA is applied to calculate module preservation statistics between gene modules build from each cell type and from other cell types. The candidate cell type-specific disease-associated gene modules whose medians of Zsummary are smaller than two are reported as cell type-specific disease-associated gene modules. 

To characterize the cell type heterogeneity of autism spectrum disorder (ASD), we applied cell type-specific gene network-based analysis and identified several cell type-specific ASD-associated gene modules. 

The codes correspond to the following paper, where further details can be found:
Guan J, Lin Y and Ji G (2020) Cell Type-Specific Gene Network-Based Analysis Depicts the Heterogeneity of Autism Spectrum Disorder. Front. Cell. Neurosci. 14:59. doi: 10.3389/fncel.2020.00059
