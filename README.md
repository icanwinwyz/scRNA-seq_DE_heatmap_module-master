# scRNA-seq_DE_heatmap_module
This is one of the modules for 10X scRNA-seq data analysis pipeline - scRNA-seq DE analysis and heatmap Module


### Input
* library required: Seurat,filesstrings,gplots,ggplots2,patchwork,dplyr
* Raw expression count matrix
* Cluster (predifined cluster used as labels of x axist)
* name for file prefix (folder and violin plots with all genes)
* min_pct (optional) minimum % of cells express such gene (default:0.1)
* logfc (optional) minimum logfc for such gene (default:0.1)
* cluster_column (optional) column used in cluster file (default:1)


### Demo data:
* demo_expression_matrix - raw expression count matrix
* demo_cluster.csv - predifined cluster used as labels of x axis
* test_DE_heatmap_plots - DE analysis results and heatmap

Please follow the format of "demo_cluster.csv" to generate the corresponding files

### Usage(on titan): 
```r
source activate my_python_3_7

Rscript Rscript scRNA_DE_heatmap_plot.R raw_expression.csv cluster.csv name <min_pct> <logfc> <column used in cluster file>
```


### Outputs:
* A folder with prefix given by user will be created and DE analysis resulst and heatmap will be generated there.

```r
> sessionInfo()
R version 3.6.2 (2019-12-12)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: Ubuntu 18.04.2 LTS

Matrix products: default
BLAS/LAPACK: /home/genomics/anaconda3/envs/my_python_3_7/lib/libopenblasp-r0.3.7.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] filesstrings_3.1.5 stringr_1.4.0      Seurat_3.1.5

loaded via a namespace (and not attached):
 [1] httr_1.4.1          tidyr_1.0.2         jsonlite_1.6.1
 [4] viridisLite_0.3.0   splines_3.6.2       lsei_1.2-0
 [7] leiden_0.3.3        gtools_3.8.1        RcppParallel_4.4.4
[10] assertthat_0.2.1    ggrepel_0.8.1       globals_0.12.5
[13] pillar_1.4.3        lattice_0.20-38     glue_1.3.1
[16] reticulate_1.14     digest_0.6.24       RColorBrewer_1.1-2
[19] colorspace_1.4-1    cowplot_1.0.0       htmltools_0.4.0
[22] Matrix_1.2-18       plyr_1.8.5          pkgconfig_2.0.3
[25] tsne_0.1-3          listenv_0.8.0       purrr_0.3.3
[28] patchwork_1.0.0     scales_1.1.0        processx_3.4.2
[31] RANN_2.6.1          gdata_2.18.0        Rtsne_0.15
[34] tibble_2.1.3        ggplot2_3.2.1       ROCR_1.0-7
[37] pbapply_1.4-2       lazyeval_0.2.2      survival_3.1-8
[40] magrittr_1.5        crayon_1.3.4        ps_1.3.2
[43] strex_1.2.0         future_1.16.0       nlme_3.1-144
[46] MASS_7.3-51.5       gplots_3.0.1.2      ica_1.0-2
[49] tools_3.6.2         fitdistrplus_1.0-14 data.table_1.12.8
[52] lifecycle_0.1.0     plotly_4.9.2        munsell_0.5.0
[55] cluster_2.1.0       irlba_2.3.3         compiler_3.6.2
[58] rsvd_1.0.2          caTools_1.18.0      rlang_0.4.4
[61] grid_3.6.2          ggridges_0.5.2      RcppAnnoy_0.0.14
[64] rappdirs_0.3.1      htmlwidgets_1.5.1   igraph_1.2.4.2
[67] bitops_1.0-6        npsurv_0.4-0        gtable_0.3.0
[70] codetools_0.2-16    reshape2_1.4.3      R6_2.4.1
[73] gridExtra_2.3       zoo_1.8-7           dplyr_0.8.4
[76] uwot_0.1.5          future.apply_1.4.0  KernSmooth_2.23-16
[79] ape_5.3             stringi_1.4.5       parallel_3.6.2
[82] Rcpp_1.0.3          vctrs_0.2.2         sctransform_0.2.1
[85] png_0.1-7           tidyselect_1.0.0    lmtest_0.9-37
```
