# Parameters files

## Main parameters file: configfile
It is the main parameters file that contains non-optionnal settings to run the pipeline.

Be careful: in a yaml file, the indentation is important.

This file is organized in 2 parts:

### 1. steps: choose the steps to run

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `steps` | steps to run | ["Alignment_countTable_GE","Droplets_QC_GE","Filtering_GE"] | No default value | "Alignment_countTable_GE", "Alignment_countTable_ADT", "Alignment_annotations_TCR_BCR", "Droplets_QC_GE", "Filtering_GE", "Norm_DimRed_Eval_GE", "Clust_Markers_Annot_GE", "Adding_ADT", "Adding_TCR", "Adding_BCR", "Cerebro" |

Note: to have more details on steps, see Pipeline details page of the wiki.

### 2. parameters for each step

* Alignment_countTable_GE:

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `sample.name.ge` | list of samples names of genes expression | ["sample1_GE", "sample2_GE"] | No default value | NA |
| `input.dir.ge` | absolute path to fastq files of genes expression | "/mnt/beegfs/userdata/m_aglave/fastq/" | No default value | NA
| `output.dir.ge` | absolute path to output folder | "/mnt/beegfs/userdata/m_aglave/pipeline/output/" | No default value | NA |
| `sctech` | technology of 10X used to generate fastq files | "10xv2" | "10xv3" | "10xv2","10xv3" |
| `kindex.ge` | absolute path to index file for the aligment of genes expression | "/mnt/beegfs/database/bioinfo/single-cell/INDEX/KB-python_KALLISTO/0.24.4_0.46.2/homo_sapiens/GRCh38/Ensembl/r99/cDNA_LINCs_MIRs/GRCH38_r99_cDNA_linc_mir.kidx" | "/mnt/beegfs/database/bioinfo/single-cell/INDEX/KB-python_KALLISTO/0.24.4_0.46.2/homo_sapiens/GRCh38/Ensembl/r99/cDNA_LINCs_MIRs/GRCH38_r99_cDNA_linc_mir.kidx" | NA |
| `tr2g.file.ge` | absolute path to tr2g file for the aligment of genes expression | "/mnt/beegfs/database/bioinfo/single-cell/INDEX/KB-python_KALLISTO/0.24.4_0.46.2/homo_sapiens/GRCh38/Ensembl/r99/cDNA_LINCs_MIRs/GRCH38_r99_cDNA_linc_mir_tr2gs.txt" | "/mnt/beegfs/database/bioinfo/single-cell/INDEX/KB-python_KALLISTO/0.24.4_0.46.2/homo_sapiens/GRCh38/Ensembl/r99/cDNA_LINCs_MIRs/GRCH38_r99_cDNA_linc_mir_tr2gs.txt" | NA |
| `reference.txt` | text for the aligment of genes expression in Materials and Methods | "Ensembl reference transcriptome v99 corresponding to the homo sapiens GRCH38 build" | "<insert_you_reference_here>" | NA |
| `fastqscreen_index` | absolute path to the configuration file of references for fastq-screen alignment | "/mnt/beegfs/database/bioinfo/single-cell/INDEX/FASTQ_SCREEN/0.14.0/fastq_screen.conf" | "/mnt/beegfs/database/bioinfo/single-cell/INDEX/FASTQ_SCREEN/0.14.0/fastq_screen.conf" | NA |

* Droplets_QC_GE:

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `sample.name.ge` | list of samples names of genes expression | ["sample1_GE", "sample2_GE"] | determined from `sample.name.ge` of `Alignment_countTable_GE` if it exists | NA |
| `input.dir.ge` | absolute path to the aligment results folder | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/KALLISTOBUS/","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/KALLISTOBUS/"] | determined from `output.dir.ge` of `Alignment_countTable_GE` if it exists | NA |
| `output.dir.ge` | absolute path to output folder | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/"] | determined from `output.dir.ge` of `Alignment_countTable_GE` if it exists | NA |
| `species` | species of genes expression | "homo_sapiens" | "homo_sapiens" | "homo_sapiens","mus_musculus" |
| `author.name` | name of the analysis author | "marine_aglave" | No default value | NA |
| `author.mail` | mail of the analysis author | "monmail@gustaveroussy.fr" | No default value | NA |
| `emptydrops.fdr` | FDR threshold for emptydrops tool | "5E-02" | "1E-03" | NA |
| `droplets.limit` | number min of droplets to run emptydrops | "1E+04" | "1E+05" | NA |
| `emptydrops.retain` | all droplets with a number of UMI above this value is considered as a cell | 1000 | No default value | NA |
| `translation` | bool to translate ENSG into Gene Symbol | TRUE | TRUE | TRUE/FALSE |
| `pcmito.min` | threshold min for percentage of mitochondrial RNA (below this threshold the cells are eliminated) | 0 | 0 | NA |
| `pcmito.max` | threshold max for percentage of mitochondrial RNA (above this threshold the cells are eliminated) | 0.1 | 0.2 | NA |
| `pcribo.min` | threshold min for percentage of ribosomal RNA (below this threshold the cells are eliminated) | 0.1 | 0 | NA |
| `pcribo.max` | threshold max for percentage of ribosomal RNA (above this threshold the cells are eliminated) | 0.9 | 1 | NA |
| `min.features` | threshold min for number of genes (below this threshold the cells are eliminated) | 150 | 200 | NA |
| `min.counts` | threshold min for number of UMI (below this threshold the cells are eliminated) | 1500 | 1000 | NA |
| `min.cells` | include genes expressed in at least this many cells (minimum cells covering) | 10 | 5 | NA |
| `mt.genes.file` | RDS file with list of mitochondrial genes | "/mnt/beegfs/pipelines/single-cell/resources/GENELISTS/homo_sapiens_mito_symbols_20191001.rds" | determined from `species` parameter of `Droplets_QC_GE` | NA |
| `crb.genes.file` | RDS file with list of ribosomal genes | "/mnt/beegfs/pipelines/single-cell/resources/GENELISTS/homo_sapiens_cribo_symbols_20191015.rds" | determined from `species` parameter of `Droplets_QC_GE` | NA |
| `str.genes.file` | RDS file with list of mecanic stress genes | "/mnt/beegfs/pipelines/single-cell/resources/GENELISTS/homo_sapiens_stress_symbols_20200224.rds" | determined from `species` parameter of `Droplets_QC_GE` | NA |
| `translation.file` | file of translation between ENSG into Gene Symbol | "/mnt/beegfs/pipelines/single-cell/resources/GENE_CONVERT/EnsemblToGeneSymbol_Homo_sapiens.GRCh38.txt" | determined from `species` parameter of `Droplets_QC_GE` | NA |

* Filtering_GE:

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `sample.name.ge` | list of samples names | ["sample1_GE", "sample2_GE"] | determined from `sample.name.ge` of `Droplets_QC_GE` if it exists | NA |
| `input.rda.ge` | absolute path to the file.rda containing the seurat R object | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/QC_droplets/sample1_GE_QC_NON-NORMALIZED.rda","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/QC_droplets/sample2_GE_QC_NON-NORMALIZED.rda"] | determined from `output.dir.ge` of `Droplets_QC_GE` if it exists | NA |
| `output.dir.ge` | absolute path to output folder | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/"] | determined from `output.dir.ge` of `Droplets_QC_GE` if it exists | NA |
| `author.name` | name of the analysis author | "marine_aglave" | No default value | NA |
| `author.mail` | mail of the analysis author | "monmail@gustaveroussy.fr" | No default value | NA |
| `pcmito.min` | threshold min for percentage of mitochondrial RNA (below this threshold the cells are eliminated) | 0 | 0 | NA |
| `pcmito.max` | threshold max for percentage of mitochondrial RNA (above this threshold the cells are eliminated) | 0.1 | 0.2 | NA |
| `pcribo.min` | threshold min for percentage of ribosomal RNA (below this threshold the cells are eliminated) | 0.1 | 0 | NA |
| `pcribo.max` | threshold max for percentage of ribosomal RNA (above this threshold the cells are eliminated) | 0.9 | 1 | NA |
| `min.features` | threshold min for number of genes (below this threshold the cells are eliminated) | 150 | 200 | NA |
| `min.counts` | threshold min for number of UMI (below this threshold the cells are eliminated) | 1500 | 1000 | NA |
| `min.cells` | include genes expressed in at least this many cells (minimum cells covering) | 10 | 5 | NA |
| `doublets.filter.method` | method used to filter doublets. To not filter set this parameter to "none" | "all" | "all" | "all","scDblFinder","scds","none" |
| `cc.seurat.file` | RDS file with list of cell cycle genes for seurat | "/mnt/beegfs/pipelines/single-cell/resources/GENELISTS/homo_sapiens_cyclone_pairs_symbols_20191001.rds" | determined from `species` into seurat object | NA |
| `cc.cyclone.file` | RDS file with list of cell cycle genes for cyclone | "/mnt/beegfs/pipelines/single-cell/resources/GENELISTS/homo_sapiens_cyclone_pairs_symbols_20191001.rds" | determined from `species` into seurat object | NA |

* Norm_DimRed_Eval_GE:

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `sample.name.ge` | list of samples names | ["sample1_GE", "sample2_GE"] | determined from `sample.name.ge` of `Filtering_GE` if it exists | NA |
| `input.rda.ge` | absolute path to the file.rda containing the seurat R object | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/sample1_GE_FILTERED_NON-NORMALIZED.rda","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/sample2_GE_FILTERED_NON-NORMALIZED.rda"] | determined from `output.dir.ge` of `Filtering_GE` if it exists | NA |
| `output.dir.ge` | absolute path to output folder | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/"] | determined from `output.dir.ge` of `Filtering_GE` if it exists | NA |
| `author.name` | name of the analysis author | "marine_aglave" | No default value | NA |
| `author.mail` | mail of the analysis author | "monmail@gustaveroussy.fr" | No default value | NA |
| `eval.markers` | list of genes to evaluate normalization and dimension reduction | "GAPDH,CD4,CD8A,CD24,CTLA4" | No default value | NA |
| `features.n` | number of High Variable Genes to consider | 3000 | 2000 | NA |
| `norm.method` | name of normalization method | "LogNormalize" | "SCTransform" | "LogNormalize","SCTransform" |
| `dimred.method` | name of dimension reduction method | "scbfa" | "pca" | "scbfa","bpca","pca","ica","mds" |
| `vtr` | list of biases to regress | "nFeature_RNA,percent_mt" | No default value | percent_mt, percent_rb, nFeature_RNA, percent_st, Cyclone.Phase, and all other column name in metadata |
| `vtr.scale` | bool to center biaises to regress (for scbfa and bpca only) | FALSE | FALSE | TRUE,FALSE |
| `dims.max` | number max of dimensions to compute | 100 | 50 | NA |
| `dims.min` | number min of dimensions to compute | 10 | 3 | NA |
| `dims.steps` | steps for dimensions to compute for evaluation | 3 | 2 | NA |
| `res.max` | number max of resolution to compute for evaluation | 3 | 1.2 | NA |
| `res.min` | number min of resolution to compute for evaluation | 0.1 | 0.1 | NA |
| `res.steps` | steps for resolution to compute for evaluation | 0.2 | 0.1 | NA |

* Clust_Markers_Annot_GE:

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `sample.name.ge` | list of samples names | ["sample1_GE", "sample2_GE"] | determined from `sample.name.ge` of `Norm_DimRed_Eval_GE` if it exists | NA |
| `input.rda.ge` | absolute path to the normalized and reduced seurat object (in .rda format) | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/sample1_GE_SCTransform_pca.rda","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/sample2_GE_SCTransform_pca.rda"] | determined from `output.dir.ge` of `Norm_DimRed_Eval_GE` if it exists | NA |
| `output.dir.ge` | absolute path to output folder | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/"] | determined from `output.dir.ge` of `Norm_DimRed_Eval_GE` if it exists | NA |
| `author.name` | name of the analysis author | "marine_aglave" | No default value | NA |
| `author.mail` | mail of the analysis author | "monmail@gustaveroussy.fr" | No default value | NA |
| `markfile` | genes to plot on umap | see makfile section of Parameters files | No default value | NA |
| `keep.dims` | number of dimension to keep for clustering (from 0 to keep.dims) | 25 | No default value | NA |
| `keep.res` | resolution value for clustering | 0.5 | No default value | NA |
| `cfr.minscore` | minimum correlation score for clustifyr to consider | 0.40 | 0.35 | NA |
| `sr.minscore` | minimum correlation score for SingleR to consider | 0.20 | 0.25 | NA |

* Cerebro:

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `input.rda` | absolute path to the seurat object (in .rda format) to convert in cerebro object | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/dims25_res0.5/sample1_GE_SCTransform_pca_25_0.5.rda","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/dims25_res0.5/sample2_GE_SCTransform_pca_25_0.5.rda"] | determined from seurat object output of `Adding_BCR` or `Adding_TCR` or `Adding_ADT` or `Clust_Markers_Annot_GE` | NA |
| `author.name` | name of the analysis author | "marine_aglave" | No default value | NA |
| `author.mail` | mail of the analysis author | "monmail@gustaveroussy.fr" | No default value | NA |
| `version` | version of cerebro to use | "v1.2" | "v1.3" | "v1.2","v1.3" |
| `groups` | column name (in meta.data) to define clusters/comparisons for Cerebro object (usefull for TCR/BCR part). | "conditions" | the last RNA clustering and samples information are already included | all column name in metadata |
| `remove.other.reductions` | remove all other reductions present in seurat object (keep only final umap) | FALSE | FALSE | TRUE,FALSE |
| `remove.other.idents` | remove all other clustering present in seurat object (keep only the last clustering) | TRUE | FALSE | TRUE,FALSE |
| `remove.mt.genes` | remove mitochondrial genes (to see better the other genes) | FALSE | FALSE | TRUE,FALSE |
| `remove.crb.genes` | remove ribosomal genes (to see better the other genes) | FALSE | FALSE | TRUE,FALSE |
| `remove.str.genes` | remove stress genes (to see better the other genes) | FALSE | FALSE | TRUE,FALSE |
| `only.pos.DE` | keep only positive DE genes from customized differential expression analysis (for genes markers identification is always only positive) | FALSE | FALSE | TRUE,FALSE |
| `remove.custom.DE` | remove results from customized differential expression analysis | FALSE | FALSE | TRUE,FALSE |
| `gmt.file` | GMT file for cerebro | "/mnt/beegfs/pipelines/single-cell/resources/DATABASE/MSIGDB/7.1/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt" | "/mnt/beegfs/pipelines/single-cell/resources/DATABASE/MSIGDB/7.1/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt" | NA |

* Alignment_countTable_ADT:

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `sample.name.adt` | list of samples names of cell surface proteins | ["sample1_ADT", "sample2_ADT"] | No default value | NA |
| `input.dir.adt` | absolute path to cell surface proteins fastq files | "/mnt/beegfs/userdata/m_aglave/fastq/" | No default value | NA |
| `output.dir.adt` | absolute path to output folder | "/mnt/beegfs/userdata/m_aglave/pipeline/output/" | No default value | NA |
| `sctech` | technology of 10X used to generate fastq files | "10xv2" | "10xv3" | "10xv2","10xv3" |
| `kindex.adt` | absolute path to index file for aligment | "/mnt/beegfs/userdata/m_aglave/ADT/kallisto_index/project_CITEseq_kallisto_index" | No default value | NA |
| `tr2g.file.adt` | absolute path to tr2g file for aligment | "/mnt/beegfs/userdata/m_aglave/ADT/kallisto_index/project_CITEseq_tr2gs.txt" | No default value | NA |

* Adding_ADT:

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `input.rda.ge` | absolute path to the seurat object (in .rda format) | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/dims25_res0.5/sample1_GE_SCTransform_pca_25_0.5.rda","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/dims25_res0.5/sample2_GE_SCTransform_pca_25_0.5.rda"] | determined from seurat object output of  `Clust_Markers_Annot_GE` | NA |
| `input.dir.adt` | absolute path to the aligment results folder of cell surface proteins | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_ADT/KALLISTOBUS/","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_ADT/KALLISTOBUS/"] | determined from `output.dir.adt` of `Alignment_countTable_ADT` | NA |
| `author.name` | name of the analysis author | "marine_aglave" | No default value | NA |
| `author.mail` | mail of the analysis author | "monmail@gustaveroussy.fr" | No default value | NA |
| `gene.names` | list of gene names wich correspond to the ADT proteins | "CD3G,CD4,CTLA4" | No default value | NA |
| `ADT.max.cutoff` | list of quantile max to cutoff protein expression for plot | "q70,q95,q85" | "q95" * number of gene in "gene.names" parameter | NA |
| `ADT.min.cutoff` | list of quantile min to cutoff protein expression for plot | "q30,q30,q55" | "q30" * number of gene in "gene.names" parameter | NA |

* Alignment_annotations_TCR_BCR:

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `sample.name.tcr` | list of samples names of TCR | ["sample1_TCR", "sample2_TCR"] | No default value | NA |
| `input.dir.tcr` | absolute path to TCR fastq files | "/mnt/beegfs/userdata/m_aglave/fastq/" | No default value | NA |
| `sample.name.bcr` | list of samples names of BCR | ["sample1_BCR", "sample2_BCR"] | No default value | NA |
| `input.dir.bcr` | absolute path to BCR fastq files | "/mnt/beegfs/userdata/m_aglave/fastq/" | No default value | NA |
| `output.dir.tcr_bcr` | absolute path to output folder | "/mnt/beegfs/userdata/m_aglave/pipeline/" | No default value | NA |
| `crindex.tcr_bcr` | CellRanger index for vdj analysis | "/mnt/beegfs/database/bioinfo/single-cell/TCR_REFERENCES/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0" | "/mnt/beegfs/database/bioinfo/single-cell/TCR_REFERENCES/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0" | NA |
| `fastqscreen_index` | absolute path to the configuration file of references for fastq-screen alignment | "/mnt/beegfs/database/bioinfo/single-cell/INDEX/FASTQ_SCREEN/0.14.0/fastq_screen.conf" | "/mnt/beegfs/database/bioinfo/single-cell/INDEX/FASTQ_SCREEN/0.14.0/fastq_screen.conf" | NA |

* Adding_TCR:

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `input.rda` | absolute path to the seurat object (in .rda format) | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/dims25_res0.5/sample1_GE_SCTransform_pca_25_0.5_ADT.rda","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/dims25_res0.5/sample2_GE_SCTransform_pca_25_0.5_ADT.rda"] | determined from seurat object output of `Adding_ADT` or `Clust_Markers_Annot_GE` | NA |
| `vdj.input.file.tcr` | file filtered_contig_annotations.csv from CellRanger aligment pipeline | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_TCR/sample1_TCR_CellRanger/outs/filtered_contig_annotations.csv","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_TCR/sample2_TCR_CellRanger/outs/filtered_contig_annotations.csv"] | determined from `output.dir.tcr_bcr` of Alignment_annotations_TCR_BCR | NA |
| `author.name` | name of the analysis author | "marine_aglave" | No default value | NA |
| `author.mail` | mail of the analysis author | "monmail@gustaveroussy.fr" | No default value | NA |

* Adding_BCR:

| Name | Description | Example | Default value | Possible values |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `input.rda` | absolute path to the seurat object (in .rda format) | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/dims25_res0.5/sample1_GE_SCTransform_pca_25_0.5_ADT_TCR.rda","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_GE/F200_C1000_M0-0.2_R0-1_G5/DOUBLETSFILTER_all/SCTransform/pca/dims25_res0.5/sample2_GE_SCTransform_pca_25_0.5_ADT_TCR.rda"] | determined from seurat object output of `Adding_TCR` or `Adding_ADT` or `Clust_Markers_Annot_GE` | NA |
| `vdj.input.file.bcr` | file filtered_contig_annotations.csv from CellRanger aligment pipeline | ["/mnt/beegfs/userdata/m_aglave/pipeline/output/sample1_BCR/sample1_BCR_CellRanger/outs/filtered_contig_annotations.csv","/mnt/beegfs/userdata/m_aglave/pipeline/output/sample2_BCR/sample2_BCR_CellRanger/outs/filtered_contig_annotations.csv"] | determined from `output.dir.tcr_bcr` of `Alignment_annotations_TCR_BCR` | NA |
| `author.name` | name of the analysis author | "marine_aglave" | No default value | NA |
| `author.mail` | mail of the analysis author | "monmail@gustaveroussy.fr" | No default value | NA |


### Notes:

* If there is not default value, the parameter is obligatory.
* sctech is a common parameter for GE and ADT
* `str.genes.file` parameter of `Droplets_QC_GE` step corresponds to a list of mecanic stress genes from the thesis of Léo Machado entitled « From skeletal muscle stem cells to tissue atlases: new tools to investigate and circumvent dissociation-induced stress », 2019.
* The dimension and resolution parameters depend on sample complexity and number of cells.
* The index for adt alignment must be specific of antibodies used. It can be made thanks to kb-python tool usable by conda.
* The list of genes name in `gene.names` parameter, of Alignment_countTable_ADT step, must be in the same ordre than proteins name in the index, to keep the correspondance.
* The `Droplets_QC_GE`, `Filtering_GE`, `Norm_DimRed_Eval_GE`, `Clust_Markers_Annot_GE`, `Adding_ADT`, `Adding_TCR`, `Adding_BCR` and `Cerebro` steps are proceded in this order. The input/output parameters are automatically detemined thanks to the step before. For exemple, if there is no `Adding_ADT` step, parameters of `Adding_TCR` step will be determined thanks to the `Clust_Markers_Annot_GE` step; if there are no `Adding_ADT`, `Adding_TCR` and `Adding_BCR` steps, parameters of `Cerebro` step will be determined thanks to the `Clust_Markers_Annot_GE` step.

### Example of parameters file:

* Example of minimal parameter file: see /mnt/beegfs/pipelines/single-cell/examples/Param_min.yaml
* Empty maximal parameter file:  see /mnt/beegfs/pipelines/single-cell/examples/Param_max.yaml
