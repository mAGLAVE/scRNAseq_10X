## INTEGRATION PROTOCOL (Seurat or scbfa)
## To Do: test Harmony and LIGER integration

#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda.int", help="Input filtred, normalized and reducted seurat object (in .rda format)."),
  make_option("--output.dir.int", help="Output path"),
  make_option("--markfile", help="Genes to plot on umap (# )format: 2 columns named Genes and Signature"),
  make_option("--author.name", help="Name of auhtor of the analysis"),
  make_option("--author.mail", help="Email of auhtor of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Analysis Parameters
  # Clustering
  make_option("--keep.dims", help="Number of dimension to keep for clustering (from 0 to keep.dims)"),
  make_option("--keep.res", help="Resolution value for clustering"),
  # Annotation
  make_option("--cfr.minscore", help="Minimum correlation score for clustifyr to consider"),
  make_option("--sr.minscore", help="Minimum correlation score for SingleR to consider"),
  ### Yaml parameters file to remplace all parameters before (usefull to use R script without snakemake)
  make_option("--yaml", help="Patho to yaml file with all parameters")
)
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)


#### Formatting Parameters ####
#convert "NULL"/"FALSE"/"TRUE" (in character) into NULL/FALSE/TRUE
for (i in names(args$options)){
  if ((length(args$options[i]) == 0) || (length(args$options[i]) == 1 && toupper(args$options[i]) == "NULL")) { args$options[i] <- NULL
  } else if ((length(args$options[i]) == 1) && (toupper(args$options[i]) == "FALSE")) { args$options[i] <- FALSE
  } else if ((length(args$options[i]) == 1) && (toupper(args$options[i]) == "TRUE")) { args$options[i] <- TRUE
  }
}

#### Get Paramaters ####
### Project
input.rda.int <- args$options$input.rda.int
output.dir.int <- args$options$output.dir.int
markfile <- if (!is.null(args$options$markfile)) unlist(stringr::str_split(args$options$markfile, ","))
author.name <- args$options$author.name
author.mail <- args$options$author.mail
### Computational Parameters
nthreads <-  if (!is.null(args$options$nthreads)) as.numeric(args$options$nthreads)
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
# Clustering
keep.dims <- if (!is.null(args$options$keep.dims)) as.numeric(args$options$keep.dims)
keep.res <- if (!is.null(args$options$keep.res)) as.numeric(args$options$keep.res)
# Annotation
cfr.minscore <- if (!is.null(args$options$cfr.minscore)) as.numeric(args$options$cfr.minscore)
sr.minscore <- if (!is.null(args$options$sr.minscore)) as.numeric(args$options$sr.minscore)
### Yaml parameters file to remplace all parameters before (usefull to use R script without snakemake)
if (!is.null(args$options$yaml)){
  yaml_options <- yaml::yaml.load_file(args$options$yaml)
  for(i in names(yaml_options)) {
    #convert "NULL"/"FALSE"/"TRUE" (in character) into NULL/FALSE/TRUE
    if ((length(yaml_options[[i]]) == 0) || (length(yaml_options[[i]]) == 1 && toupper(yaml_options[[i]]) == "NULL")) { yaml_options[[i]] <- NULL
    } else if ((length(yaml_options[[i]]) == 1) && (toupper(yaml_options[[i]]) == "FALSE")) { yaml_options[[i]] <- FALSE
    } else if ((length(yaml_options[[i]]) == 1) && (toupper(yaml_options[[i]]) == "TRUE")) { yaml_options[[i]] <- TRUE
    }
    #assign values
    assign(i, yaml_options[[i]])
    if(i %in% c("nthreads","keep.dims","keep.res","cfr.minscore","sr.minscore")) assign(i, as.numeric(yaml_options[[i]]))else assign(i, yaml_options[[i]])
    
  }
  rm(yaml_options, i)
}
### Clean
rm(option_list,parser,args)

#### Get path if snakemake/singularity/local ####
if(is.null(pipeline.path)) stop("--pipeline.path parameter must be set!")

#### Check non-optional parameters ####
if (is.null(input.rda.int)) stop("input.rda.int parameter can't be empty!")
if (is.null(output.dir.int)) stop("output.dir.int parameter can't be empty!")
if (is.null(keep.dims)) stop("keep.dims parameter can't be empty!")
if (is.null(keep.res)) stop("keep.res parameter can't be empty!")

### Load data
load(input.rda.int)

### Save project parameters
if (!is.null(author.name) && !tolower(author.name) %in% tolower(sobj@misc$params$author.name)) sobj@misc$params$author.name <- c(sobj@misc$params$author.name, author.name)
if (!is.null(author.mail) && !tolower(author.mail) %in% tolower(sobj@misc$params$author.mail)) sobj@misc$params$author.mail <- c(sobj@misc$params$author.mail, author.mail)

#### Get Missing Paramaters ####
### Project
sample.name.int <- sobj@misc$params$sample.name.int
species <- sobj@misc$params$species
### Computational Parameters
if (is.null(nthreads)) nthreads <- 4
### Analysis Parameters
# Integration
integration.method <- sobj@misc$params$integration$method
# Normalization and dimension reduction
assay <- sobj@misc$params$integration$out.assay
norm.method <- sobj@assays[[assay]]@misc$params$normalization$normalization.method
norm_vtr <- paste0(c(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr[1])) paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = '_') else NULL), collapse = '_')
dimred.method <- sobj@assays[[assay]]@misc$params$reductions$method
if(integration.method == "Seurat") red.name <- paste(c("integrated", dimred.method, integration.method), collapse = '_')
if(integration.method %in% c('scbfa', 'bpca','Liger')) red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
if(integration.method == "Harmony") red.name <- paste(c(assay, dimred.method), collapse = '_')
dimred_vtr <- paste0(c(dimred.method, if(!is.na(sobj@reductions[[red.name]]@misc$vtr[1])) paste(sobj@reductions[[red.name]]@misc$vtr, collapse = '_') else NULL), collapse = '_')
if(integration.method == "Harmony") red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
# Annotation
if (is.null(cfr.minscore)) cfr.minscore <- 0.35
if (is.null(sr.minscore)) sr.minscore <- 0.25

#### Fixed parameters ####
# Annotation
if (species == "homo_sapiens") {
  singler.setnames <- c("HumanPrimaryCellAtlasData", "NovershternHematopoieticData", "DatabaseImmuneCellExpressionData", "MonacoImmuneData")
  clustifyr.setnames <- c("pbmc_avg", "hema_microarray_matrix", "gtex_bulk_matrix")
}
if (species == "mus_musculus") {
  singler.setnames <- c("MouseRNAseqData", "ImmGenData")
  clustifyr.setnames <- c("ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs", "ref_moca_main", "ref_immgen", "ref_mouse.rnaseq")
}
if (species == "rattus_norvegicus") {
  singler.setnames <- c("MouseRNAseqData", "ImmGenData")
  clustifyr.setnames <- c("ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs", "ref_moca_main", "ref_immgen", "ref_mouse.rnaseq")
}

#### Fixed parameters ####
# Plots
solo.pt.size <- 3
multi.pt.size <- 2
gradient.cols <- c("gold", "blue")

#### Get genes markers ####
if (is.null(markfile)){
  markers <- NULL
}else{
  markers <- c()
  for (i in markfile){
    mark.xl <- openxlsx::read.xlsx(i, sheet = 1, startRow = 1, fillMergedCells = TRUE, colNames = TRUE)
    mark.xl <- mark.xl[order(mark.xl[,2]),]
    markers_tmp <- setNames(mark.xl[,1], mark.xl[,2])
    markers <- c(markers,markers_tmp)
  }
}

#### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

## RUN
######

print("#####################################")
print(paste0("Sample: ", sample.name.int))
print(paste0("RDA file: ", input.rda.int))
print(paste0("Dimension: ", keep.dims))
print(paste0("Resolution: ", keep.res))
print("#####################################")

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### Building clustered output directory
clust.dir <- paste(output.dir.int, paste0("dims", keep.dims, "_res", keep.res), sep = '/')
dir.create(path = clust.dir, recursive = TRUE, showWarnings = FALSE)

### Replotting final clusters
cat("\nClustering...\n")
sobj <- louvain.cluster(sobj = sobj, reduction = red.name, max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size, algorithm = 1)

## Setting ident name and RNA.reduction
ident.name <- paste0(paste0(red.name, ".", keep.dims), '_res.', stringr::str_replace(keep.res, pattern = ",", replacement = "."))
INT.reduction <- paste(c(red.name, keep.dims, 'umap'), collapse = '_')
#sobj@reductions[[red.name]]@misc$from.assay <- assay

### uMAP plot by sample
cat("\nuMAP plot by sample...\n")
blockpix = 600
png(filename = paste0(clust.dir, '/', paste(c(sample.name.int, red.name, 'uMAP.png'), collapse = '_')), width = 1000, height = 1000)
print(Seurat::DimPlot(object = sobj, reduction = paste(c(red.name, keep.dims, 'umap'), collapse = '_'), order = sample(x = 1:ncol(sobj), size = ncol(sobj), replace = FALSE), group.by = 'orig.ident', pt.size = solo.pt.size) + ggplot2::ggtitle("uMAP for all samples ") + Seurat::DarkTheme())
dev.off()
grid.xy <- grid.scalers(length(unique(sobj@meta.data$orig.ident)))
png(filename = paste0(clust.dir, '/', paste(c(sample.name.int, red.name, 'split', 'uMAP.png'), collapse = '_')), width = grid.xy[1]*blockpix, height = grid.xy[2]*blockpix)
print(Seurat::DimPlot(object = sobj, reduction = paste(c(red.name, keep.dims, 'umap'), collapse = '_'), group.by = ident.name, split.by = 'orig.ident', pt.size = solo.pt.size, ncol = grid.xy[1]) + ggplot2::ggtitle(paste0("uMAP split on samples")) + Seurat::DarkTheme())
dev.off()

### Technical plots
cat("\nSaving technical plots...\n")
technical.plot(sobj = sobj, ident = ident.name, out.dir = clust.dir, multi.pt.size = multi.pt.size)

### Finding markers
cat("\nFinding markers...\n")
sobj <- find.markers.quick(sobj = sobj, ident = ident.name, test.use = 'wilcox', min.pct = .75, logfc.threshold = .5, only.pos = TRUE, adjp.p.max = 5E-02, topn = 10, heatmap.cols = gradient.cols, out.dir = clust.dir)

### Automatic cell type annotation
cat("\nAutomatic cell type annotation...\n")
sobj <- cells.annot(sobj = sobj, ident = ident.name, singler.setnames = singler.setnames, clustifyr.setnames = clustifyr.setnames, sr.minscore = .25, cfr.minscore = .35, out.dir = clust.dir, solo.pt.size = solo.pt.size)

### Assessing clusters : Plotting provided marker genes
cat("\nPlotting provided marker genes...\n")
if(!is.null(markers)) sobj <- markers.umap.plot(sobj = sobj, markers = markers, ident = ident.name, out.dir = clust.dir, dimplot.cols = gradient.cols, multi.pt.size = 2)

### Materials and Methods
sobj@misc$parameters$Materials_and_Methods$Integration_Clust_Markers_Annot <- paste0(
  "An automatic annotation of cell types was perfom by SingleR (version ",sobj@misc$technical_info$SingleR,") (with fine-tuning step) and ClustifyR (version ",sobj@misc$technical_info$clustifyr,"), using packages built-in references. It labels clusters (or cells) from a dataset based on similarity (Spearman correlation score) to a reference dataset with known labels. The labels with a correlation score greater than ",sr.minscore," for SingleR or greater than ",cfr.minscore," for ClustifyR were kept.",
  "Marker genes for Louvain clusters and samples, were identified through a «one versus others» differential anaylisis using the Wilcoxon test through the FindAllMarkers() function from Seurat, considering only genes with a minimum log fold-change of 0.5 in at least 75% of cells from one of the groups compared, and FDR-adjusted p-values <0.05 (Benjaminin-Hochberg method)."
)
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
write_MandM(sobj=sobj, output.dir=clust.dir)

### Saving final object
cat("\nSaving object...\n")
GE_file=paste0(clust.dir, '/', paste(c(sample.name.int, norm_vtr, dimred_vtr, keep.dims, keep.res), collapse = "_"))
save(sobj, file = paste0(GE_file, '.rda'), compress = "bzip2")

# #TCR view
# #Fusion of *highlight_aa_top10_freq TCR_highlight_aa_top11to20_freq *highlight_aa_top11to20_freq (only TCR)
# col.names.top <- grep("highlight_aa_top", names(sobj@meta.data), value=TRUE)
# col.name.all <- grep("highlight_aa_all", names(sobj@meta.data), value=TRUE)
# for (i in 1:length(sobj@meta.data[[col.name.all]])){
#   if (is.na(sobj@meta.data[[col.names.top[1]]][i])){
#     if (is.na(sobj@meta.data[[col.names.top[2]]][i])){
#       sobj@meta.data$TCR_highlight_aa_top20_freq[i] <- NA
#     }else{
#       sobj@meta.data$TCR_highlight_aa_top20_freq[i] <- sobj@meta.data[[col.names.top[2]]][i]
#     }
#   }else if (is.na(sobj@meta.data[[col.names.top[2]]][i])){
#     sobj@meta.data$TCR_highlight_aa_top20_freq[i] <- sobj@meta.data[[col.names.top[1]]][i]
#   }else{
#       stop(paste0("Error :", col.names.top[1]," AND ", col.names.top[2]," are full in ",i,"! "))
#     }
# }
