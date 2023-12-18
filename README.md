# Cognition-associated long noncoding RNAs are dysregulated upon severe COVID-19
This repository provides instructions, code, and supplementary data files to reproduce results and figures from our lncRNA COVID-19 brain paper 
> JD Lee, IH Solomon, FJ Slack, M Mavrikaki. Cognition-associated long noncoding RNAs are dysregulated upon severe COVID-19, Front Immunol (In revision 2023). 

## Abstract

Severe COVID-19 leads to widespread transcriptomic changes in the human brain, mimicking diminished cognitive performance. As long noncoding RNAs (lncRNAs) play crucial roles in the regulation of gene expression, identification of the lncRNAs differentially expressed upon COVID-19 may nominate key regulatory nodes underpinning cognitive changes. Here we identify hundreds of lncRNAs differentially expressed in the brains of COVID-19 patients relative to uninfected age/sex-matched controls, many of which are associated with decreased cognitive performance and inflammatory cytokine response. Our analyses reveal pervasive transcriptomic changes in lncRNA expression upon severe COVID-19, which may serve as key regulators of neurocognitive changes in the brain.

## Code information

The available [R markdown file](https://github.com/jonathandlee12/covid19-brain-lnc/blob/main/covid19_lncrna_rscripts_final.Rmd) contains R scripts used to generate all figures and data tables available with the current iteration of the manuscript. The code was last run 12/17/2023. 

## Resources for data preprocessing
- [Ensembl v104 reference transcriptome, protein-coding genes](http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz)
- [Ensembl v104 reference transcriptome, noncoding RNA genes](http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz)
- [Frontal cortex transcriptome profiles from this manuscript (.fastq)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188847)

NOTE: while similar .sf files are available at GEO (from Mavrikaki et al. Nature Aging 2022), the alignments/quantifications do not include noncoding RNAs.

## Resources for data analysis
- [Cohort characteristics - COVID-19](https://github.com/jonathandlee12/covid19-brain-lnc/blob/main/TableS1_covid.csv) may also be found in Supplementary Table 1
- [Transcript-to-gene mappings (ensembl v96, from biomaRt)](https://github.com/jonathandlee12/covid19-brain-lnc/blob/main/tx2gene.txt)
- [ROSMAP preprocessed RNA-seq dataset](https://www.synapse.org/#!Synapse:syn8691134)
- [ROSMAP clinical annotation](https://www.synapse.org/Portal.html#!Synapse:syn3157322)
- [ROSMAP biospecimen metadata](https://www.synapse.org/Portal.html#!Synapse:syn3157322)

TCGA RNA-seq files are available at: https://www.cancer.gov/tcga


## Session Info from R

R version 4.0.4 (2021-02-15)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 14.1.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] viridis_0.6.0               viridisLite_0.4.0           Rtsne_0.15                  DESeq2_1.30.1              
 [5] SummarizedExperiment_1.20.0 MatrixGenerics_1.2.1        matrixStats_0.58.0          GenomicRanges_1.42.0       
 [9] GenomeInfoDb_1.26.7         tximport_1.18.0             pheatmap_1.0.12             reshape2_1.4.4             
[13] dplyr_1.0.5                 org.Hs.eg.db_3.12.0         biomaRt_2.46.3              fgsea_1.16.0               
[17] GO.db_3.12.1                AnnotationDbi_1.52.0        IRanges_2.24.1              S4Vectors_0.28.1           
[21] Biobase_2.50.0              BiocGenerics_0.36.1         ggrepel_0.9.1               ggplot2_3.3.3              
[25] stringr_1.4.0              

loaded via a namespace (and not attached):
 [1] colorspace_2.0-0       ellipsis_0.3.2         XVector_0.30.0         rstudioapi_0.13        farver_2.1.0          
 [6] bit64_4.0.5            fansi_0.4.2            mvtnorm_1.1-1          apeglm_1.12.0          xml2_1.3.2            
[11] splines_4.0.4          cachem_1.0.4           geneplotter_1.68.0     knitr_1.33             jsonlite_1.7.2        
[16] annotate_1.68.0        dbplyr_2.1.1           readr_1.4.0            compiler_4.0.4         httr_1.4.2            
[21] assertthat_0.2.1       Matrix_1.3-4           fastmap_1.1.0          cli_2.5.0              htmltools_0.5.2       
[26] prettyunits_1.1.1      tools_4.0.4            coda_0.19-4            gtable_0.3.0           glue_1.4.2            
[31] GenomeInfoDbData_1.2.4 rappdirs_0.3.3         fastmatch_1.1-0        Rcpp_1.0.7             bbmle_1.0.23.1        
[36] vctrs_0.3.8            xfun_0.31              lifecycle_1.0.0        XML_3.99-0.6           zlibbioc_1.36.0       
[41] MASS_7.3-54            scales_1.1.1           hms_1.0.0              RColorBrewer_1.1-2     yaml_2.2.1            
[46] curl_4.3.1             memoise_2.0.0          gridExtra_2.3          emdbook_1.3.12         bdsmatrix_1.3-4       
[51] stringi_1.5.3          RSQLite_2.2.7          genefilter_1.72.1      BiocParallel_1.24.1    rlang_0.4.11          
[56] pkgconfig_2.0.3        bitops_1.0-7           evaluate_0.14          lattice_0.20-44        purrr_0.3.4           
[61] labeling_0.4.2         bit_4.0.4              tidyselect_1.1.1       plyr_1.8.7             magrittr_2.0.1        
[66] R6_2.5.0               generics_0.1.0         DelayedArray_0.16.3    DBI_1.1.1              pillar_1.6.0          
[71] withr_2.4.2            survival_3.2-11        RCurl_1.98-1.3         tibble_3.1.1           crayon_1.4.1          
[76] utf8_1.2.1             BiocFileCache_1.14.0   rmarkdown_2.14         progress_1.2.2         locfit_1.5-9.4        
[81] grid_4.0.4             data.table_1.14.0      blob_1.2.1             digest_0.6.27          xtable_1.8-4          
[86] numDeriv_2016.8-1.1    openssl_1.4.4          munsell_0.5.0          askpass_1.1
