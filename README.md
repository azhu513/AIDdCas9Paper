# AIDdCas9paper analysis code

* `functions` - helper functions throughout the analysis
    * `utils.R` - data processing and cleaning functions
    * `model_helper.R` - helper functions to apply DESeq2 and apeglm 
    * `ggplot_themes.R` - ggplot themes
* `Tiling` - analysis code for the EGFR/BRAF tiling experiment
    * `EGFRBRAF_tiling_bulkcount.R` - model the bulk count in the tiling experiment
* `MissionBio` - analysis code for the Mission Bio validation experiment
    * `data_process.R` - process the txt files of single cell amplicon-seq from crisprscan workflow
    * `annotate_raw.sh` - annotate the single cell amplicon-seq alleles with amino acids change
    * `single_amino_process.R` - process the single cell data to have single amino level information
    * `allele_count_pseudobulk_model.R` - model the pseudobulk allele count
    * `sc_cell_count_model.R` - model the cell count with single amino level data
* `Plot` - code for generating the figures
    * `Fig1_SuppFig1_SuppFig2.R` - Figure 1 & Supplementary Figure 1 & 2
    * `Fig2_SuppFig4.R` - Figure 2 and Supplementary Figure 4
    * `SuppFig3.R` - Supplementary Figure 3

Open source software/packages used and their license:
    * R version 4.3.0 (GPL-2 | GPL-3)
    * ggplot2 3.4.4 (MIT License)
    * dplyr 1.1.3 (MIT License)
    * tidyr 1.3.0 (MIT License)
    * DESeq2 1.41.12 (LGPL (>= 3))
    * apeglm 1.23.1 (GPL-2)
    * data.table 1.14.8 (MIT license)
    * stringr 1.5.0 (MIT license)
    * ggpubr 0.6.0 (GPL (>= 2))
    * reshape 0.8.9 (MIT license)
    * reshape2 1.4.4 (MIT license)
    * readr 2.1.4 (MIT license)
    
The software/packages can be obtained from the [official R website](https://www.r-project.org/), [The Comprehensive R Archive Network](https://cran.r-project.org/) or [BioConductor](https://www.bioconductor.org/) websites. 