# contains helper functions for analysis notebooks

# this used to collect summary stats for every dataset
count_all <- function(step, object){
  cat(count_stats(paste("counts per", step), object@meta.data$nCount_RNA),
      count_stats(paste("features per", step), object@meta.data$nFeature_RNA),
      count_stats(paste("mitochondrial proportion per", step), object@meta.data$Mito_proportion),
      count_stats(paste("ribosomal proportion per", step), object@meta.data$Ribo_proportion), sep = "\n")
}

count_all_atac <- function(step, object){
  cat(count_stats(paste("fragments per", step), object@meta.data$passed_filters),
      count_stats(paste("fragments in peaks per", step), object@meta.data$peak_region_fragments),
      count_stats(paste("fraction reads in peaks per", step), object@meta.data$fraction_reads_in_peaks),
      count_stats(paste("nucleosome signal per", step), object@meta.data$nucleosome_signal),
      count_stats(paste("TSS enrichment per", step), object@meta.data$TSS.enrichment), sep = "\n")
}

count_stats <- function(x, y){
  cat(paste("Min", x, round(min(y), 2)),
      paste("Max", x, round(max(y), 2)),
      paste("Mean", x, round(mean(y), 2)),
      paste("Median", x, round(median(y), 2)),
      paste("Std. dev.", x, round(sd(y), 2)),
      paste("Std. error", x, round(sd(y)/sqrt(length(y)), 2)),
      paste(),
      sep = "\n")
}

# Find markers
# Taken from Warren Lab package, does FindAllMarkers but with conserved variables
FindAllConservedMarkers <- function(
  seurat,
  grouping.var,
  test.use = 'wilcox',
  verbosity = TRUE
) {
  Reduce(function(df, cluster) {
    # FindConservedMarkers throws an error if it cannot find any
    # conserved markers for a cluster, so catch those errors and
    # deal with them by just not adding any lines to the dataframe
    # for that cluster and outputting a message.
    tryCatch({
      df2 <- FindConservedMarkers(
        seurat,
        assay = "SCT",
        ident.1 = cluster,
        grouping.var = grouping.var,
        only.pos = TRUE,
        logfc.threshold = 0.585,
        test.use = test.use,
        verbose = verbosity
      )
      df2$cluster <- cluster
      df2$feature <- rownames(df2)
      return(plyr::rbind.fill(df, df2))
    }, error = function(e) {
      print(paste0('Could not find biomarkers for cluster ', cluster))
      return(df)
    })
  }, levels(Idents(seurat)), data.frame())
}

## Cluster markers (wilcoxon, single-cell)
markers_sc_wilcox <- function(x, label, verbosity = TRUE){
  Idents(x) <- label
  marker_list <- FindAllConservedMarkers(seurat = x, grouping.var = "Region")
  marker_list <- cbind.data.frame(marker_list$feature, marker_list$minimump_p_val, marker_list$cluster)
  colnames(marker_list) <- c("gene", "combined_p_val", "cluster")
  marker_list <- split(marker_list, marker_list$cluster)
  marker_list <- lapply(marker_list, function(x) as.data.frame(x[,1:2]))
  return(marker_list)
}
## Region markers (edgeR-LRT, pseudobulk)
markers_pb_edgeRLRT <- function(x, label){
  marker_lists <- vector("list", length(levels(as.factor(x$Region))))
  names(marker_lists) <- levels(as.factor(x$Region))
  for (i in 1:length(levels(as.factor(x$Region)))){
    metadata_col <- ifelse(x$Region == levels(as.factor(x$Region))[i], levels(as.factor(x$Region))[i], "others")
    x <- AddMetaData(x, metadata = metadata_col, col.name = "test_regions")
    marker_list <- run_de(x, cell_type_col = label, replicate_col = "orig.ident", label_col = "test_regions")
    marker_list <- marker_list[,1:5]
    marker_list <- split(marker_list, marker_list$cell_type)
    marker_list <- lapply(marker_list, function(x) as.data.frame(x[,2:5]))
    for (j in 1:length(marker_list)){
      rownames(marker_list[[j]]) <- marker_list[[j]]$gene
      marker_list[[j]] <- marker_list[[j]][,2:4]
    }
    marker_lists[[i]] <- marker_list
  }
  return(marker_lists)
}

## Region markers (edgeR-LRT, pseudobulk)
markers_pb_edgeRLRT_pairwise <- function(x, label){
  Idents(x) <- label
  deg_list <- vector(mode = "list", length = 4)
  names(deg_list) <- levels(as.factor(x$Region))

  for (i in 1:4){
    deg_list[[i]] <- vector(mode = "list", length = 4)
    names(deg_list[[i]]) <- levels(as.factor(x$Region))
    for (j in 1:4){
      if (i == j){
        marker_list <- cbind.data.frame(0, 0, 0)
        names(marker_list) <- c("avg_logFC", "p_val", "p_val_adj")
        deg_list[[i]][[j]] <- replicate(length(as.factor(levels(Idents(x)))), marker_list, simplify = FALSE)
        names(deg_list[[i]][[j]]) <- as.factor(levels(Idents(x)))
      }
      else {
        array_subset <- x[, as.factor(x$Region) %in% c(levels(as.factor(x$Region))[i], levels(as.factor(x$Region))[j])]
        marker_list <- run_de(array_subset, cell_type_col = label, replicate_col = "orig.ident", label_col = "Region")
        marker_list <- marker_list[,1:5]
        marker_list <- split(marker_list, marker_list$cell_type)
        marker_list <- lapply(marker_list, function(x) as.data.frame(x[,2:5]))
        for (k in 1:length(marker_list)){
          rownames(marker_list[[k]]) <- marker_list[[k]]$gene
          marker_list[[k]] <- marker_list[[k]][,2:4]
        }
        deg_list[[i]][[j]] <- marker_list
      }
    }
  }

  return(deg_list)
}

targeted_comp_markers <- function(x, i1, i2, gene_vector){

  gene_list <- FindMarkers(x, test.use = "t", only.pos = TRUE, logfc.threshold = 0.5, verbose = FALSE,
                           ident.1 = i1, ident.2 = i2)
  gene_list <- gene_list[gene_list$p_val_adj < 0.05,]
  gene_list <- gene_list[gene_list$pct.1 > 2*gene_list$pct.2,]
  gene_list <- gene_list[!rownames(gene_list) %in% gene_vector,]
}

pc_list <- c(1, 2, 3, 50)

# Dimensionality reduction IDs
dr_params <- cbind.data.frame(c("PCA", "tSNE", "UMAP"),
                                  c("pca", "tsne", "umap"))
colnames(dr_params) <- c("name", "type")


vln_metrics <- cbind.data.frame(c("nCount_RNA", "nFeature_RNA", "Mito_proportion", "Ribo_proportion"),
                                    c("UMIs", "Genes", "% Mito", "% Ribo"))
colnames(vln_metrics) <- c("id", "label")

batch_region <- c("Batch", "Region")
identity_region <- c("orig.ident", "Region")
tissue_class <- c("Tissue", "Class")# contains helper functions for analysis notebooks

# this used to collect summary stats for every dataset
count_all <- function(step, object){
  cat(count_stats(paste("counts per", step), object@meta.data$nCount_RNA),
      count_stats(paste("features per", step), object@meta.data$nFeature_RNA),
      count_stats(paste("mitochondrial proportion per", step), object@meta.data$Mito_proportion),
      count_stats(paste("ribosomal proportion per", step), object@meta.data$Ribo_proportion), sep = "\n")
}

count_stats <- function(x, y){
  cat(paste("Min", x, round(min(y), 2)),
      paste("Max", x, round(max(y), 2)),
      paste("Mean", x, round(mean(y), 2)),
      paste("Median", x, round(median(y), 2)),
      paste("Std. dev.", x, round(sd(y), 2)),
      paste("Std. error", x, round(sd(y)/sqrt(length(y)), 2)),
      paste(),
      sep = "\n")
}

create_dendrogram <- function(array_id, label_id, color_id){

  # figure out cluster medians
  medians <- do.call("cbind",
                     tapply(names(label_id),
                            label_id,
                            function(x){matrixStats::rowMedians(as.matrix(GetAssayData(object = array_id, slot = "data")[,x]))
                            }))

  pvclust_matrix <- pvclust(data = medians, method.dist = "cor", method.hclust = "average", nboot = 100, parallel = TRUE)

  dend <- as.dendrogram(pvclust_matrix$hclust)

  color_labels <- setNames(color_id, levels(label_id))

  dend <- dend %>% set("labels_cex", 0.7)
  dend <- dend %>% set("leaves_pch", 19) %>% set("leaves_cex", 0.5)
  dend <- dend %>% set("labels_col", color_labels[labels(dend)])
  dend <- dend %>% set("leaves_col", color_labels[labels(dend)])

  output <- list(medians, dend)
  names(output) <- c("Medians", "Dendrogram")

  return(output)
}

test_region_specific_expression <- function(x, gene){

  gene_expression <- FetchData(x, vars = c(gene, "orig.ident", "Region"))
  gene_expression <- aggregate(gene_expression[[gene]] ~ Region + orig.ident, data = gene_expression, FUN = "mean")

  colnames(gene_expression) <- c("Region", "Batch", gene)

  anova_gene <- aov(gene_expression[[gene]] ~ Region, data = gene_expression)

  gene_test_list <- list(gene_expression,
                         summary(anova_gene),
                         TukeyHSD(anova_gene))
  names(gene_test_list) <- c("Summary", "ANOVA", "Tukey test")

  return(gene_test_list)
}

pseudobulk_pca <- function(x, identifier, label, colors){

  pseudobulk_array <- AggregateExpression(x, group.by = identifier, assays = "RNA", slot = "counts", return.seurat = TRUE)
  pseudobulk_array <- AddMetaData(pseudobulk_array, metadata = as.factor(names(pseudobulk_array$orig.ident)), col.name = label)
  pseudobulk_array <- NormalizeData(pseudobulk_array, verbose = FALSE) # results are slightly different
  pseudobulk_array <- FindVariableFeatures(pseudobulk_array, nfeatures = 3000)
  pseudobulk_array <- ScaleData(pseudobulk_array, verbose = FALSE)
  pseudobulk_array <- RunPCA(pseudobulk_array, verbose = FALSE, npcs = 2)

  var_explained <- pseudobulk_array@reductions$pca@stdev ^ 2 /
    sum(matrixStats::rowVars(GetAssayData(pseudobulk_array, assay = "RNA", slot = "scale.data")))
  var_explained <- 100 * round(var_explained, 3)

  plot <- DimPlot_scCustom(pseudobulk_array, reduction = 'pca', group.by = label, pt.size = 5) +
    xlab(paste0("PC 1 (", var_explained[1], "%)")) + ylab(paste0("PC 2 (", var_explained[2], "%)")) +
    scale_colour_manual(values = colors) +
    theme_classic() +
    theme(plot.title = element_blank(), legend.text=element_text(size=8))
  return(plot)
}

d1d2_coexpression_stats <- function(x, region){

  subset_array <- x[, as.factor(x$Region) %in% region]
  subset_array <- FetchData(subset_array, slot = "counts", vars = c("Drd1", "Drd2"))

  print(table(subset_array$Drd1, subset_array$Drd2))
  print(table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1))

  print(paste(region, "D1+/D2-:",
              table(subset_array$Drd1 >= 1, subset_array$Drd2 == 0)[2,2],
              round(table(subset_array$Drd1 >= 1, subset_array$Drd2 == 0)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
  print(paste(region, "D1-/D2+:",
              table(subset_array$Drd1 == 0, subset_array$Drd2 >= 1)[2,2],
              round(table(subset_array$Drd1 == 0, subset_array$Drd2 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
  print(paste(region, "D1+/D2+:",
              table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1)[2,2],
              round(table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
}

d1d2_coexpression_stats_batch <- function(x, batch){

  subset_array <- x[, as.factor(x$orig.ident) %in% batch]
  subset_array <- FetchData(subset_array, slot = "counts", vars = c("Drd1", "Drd2"))

  print(table(subset_array$Drd1, subset_array$Drd2))
  print(table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1))

  print(paste(batch, "D1+/D2-:",
              table(subset_array$Drd1 >= 1, subset_array$Drd2 == 0)[2,2],
              round(table(subset_array$Drd1 >= 1, subset_array$Drd2 == 0)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
  print(paste(batch, "D1-/D2+:",
              table(subset_array$Drd1 == 0, subset_array$Drd2 >= 1)[2,2],
              round(table(subset_array$Drd1 == 0, subset_array$Drd2 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
  print(paste(batch, "D1+/D2+:",
              table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1)[2,2],
              round(table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
}

d1d2_coexpression_stats_celltype <- function(x, region, grouping, cell_type){

  subset_array <- x[, as.factor(x$Region) %in% region]
  subset_array <- FetchData(subset_array, slot = "counts", vars = c("Drd1", "Drd2", grouping))
  colnames(subset_array) <- c("Drd1", "Drd2", "group")
  subset_array <- subset_array[subset_array$group %in% cell_type,]
  subset_array <- subset_array[,1:2]

  print(table(subset_array$Drd1, subset_array$Drd2))
  print(table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1))

  print(paste(region, cell_type, "D1+/D2-:",
              table(subset_array$Drd1 >= 1, subset_array$Drd2 == 0)[2,2],
              round(table(subset_array$Drd1 >= 1, subset_array$Drd2 == 0)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
  print(paste(region, cell_type, "D1-/D2+:",
              table(subset_array$Drd1 == 0, subset_array$Drd2 >= 1)[2,2],
              round(table(subset_array$Drd1 == 0, subset_array$Drd2 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
  print(paste(region, cell_type, "D1+/D2+:",
              table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1)[2,2],
              round(table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
}

d1d2_coexpression_stats_4group <- function(x, region_id){

  subset_array <- x[, as.factor(x$Region) %in% region_id]
  subset_array <- FetchData(subset_array, slot = "counts", vars = c("Drd1", "Drd2", "Class"))
  colnames(subset_array) <- c("Drd1", "Drd2", "group")
  subset_array <- subset_array[subset_array$group %in% c("GABA.D1", "GABA.D2"),]

  d1d2_vector <- NULL
  for (i in 1:length(rownames(subset_array))){
    if (subset_array[i,1] >= 1) {
      if (subset_array[i,2] >= 1){
        d1d2_vector[i] <- "dual D1/D2"
      } else {
        d1d2_vector[i] <- "D1"
      }
    } else {
      if (subset_array[i,2] >= 1){
        d1d2_vector[i] <- "D2"
      } else {
        d1d2_vector[i] <- "no D1/D2"
      }
    }
  }
  stat_table <- table(d1d2_vector, subset_array$group)
  stat_table <- stat_table[,1:2]

return(stat_table)
}

pc_list <- c(1, 2, 3, 30)

# Dimensionality reduction IDs
dr_params <- cbind.data.frame(c("PCA", "tSNE", "UMAP"),
                              c("pca", "tsne", "umap"))
colnames(dr_params) <- c("name", "type")

vln_metrics <- cbind.data.frame(c("nCount_RNA", "nFeature_RNA", "Mito_proportion", "Ribo_proportion"),
                                c("UMIs", "Genes", "% Mito", "% Ribo"))
colnames(vln_metrics) <- c("id", "label")

batch_region <- c("Batch", "Region")
identity_region <- c("orig.ident", "Region")

marker_features <- c("Syp", "Eno2", # pan-neuronal
                     "Slc17a7", "Slc17a6", # glutamatergic neuron
                     "Gad1", "Gad2", # GABAergic neuron
                     "Drd1", "Drd2", "Adora2a", # dopaminergic neuron markers
                     "Aldh1l1", "Aqp4", # astrocyte
                     "Tmem119", "Ptprc", # microglia
                     "Pdgfra", "Cspg4", #OPC
                     "Bmp4", "Enpp6", # NFOL
                     "Mog", "Mal", # oligodendrocyte/OLG
                     "Cldn5", "Adgrl4", # endothelial cell
                     "Kcnj8", "Abcc9", # mural cell
                     "Ccdc153", "Tmem212") # ependymal cell

level2_marker_features <- c(
                            # subclusters
                            "Gad2", # confirm GABAergic
                            "Drd1", # generalized D1
                            "Drd2", # generalized D2
                            "Slc17a7", "Slc17a6", # generalized glutamatergic
                            "Kremen1", "Sema5b", "Tshz1", "Tac1", "Pdyn", "Oprm1", # generalized patch
                            "Id4", "Epha4", "Sgk1", "Wnt2", "Gng2", # generalized matrix
                            "Necab1", # D1 patch
                            "Asic4", # D2 patch
                            "Chat", #ChAT neurons
                            "Npy", # NPY neurons
                            "Sst", "Chodl", # NPY-SST neurons
                            "Mia", "Car4", # NPY-Car4 neurons
                            "Pthlh", # Pthlh neurons
                            "Th", # Th neurons
                            "Otof", "Cacng5", "Pcdh8", "Casz1", "Col11a1", "Adarb2",  # eSPN
                            # within region
                            "Cdh6", # ASt general
                            "Itga6", "Arhgef25", # ASt D1
                            "Crtac1", "Tesk2", # ASt D2
                            "Gabrg1", # CeA D1
                            "Pou3f1", # DS D1
                            #NONE FOR TS D1
                            "Dpyd", # CeA D2
                            "C79798", # DS D2
                            #NONE FOR TS D2
                            "Tnfaip8", # CeA Astro
                            "Grip1") # DS Astro

level2_marker_features_subset <- c(
  # subclusters
  "Gad2", # confirm GABAergic
  "Drd1", # generalized D1
  "Drd2", # generalized D2
  "Slc17a7", "Slc17a6", # generalized glutamatergic
  "Kremen1", "Sema5b", "Tshz1", "Tac1", "Pdyn", "Oprm1", # generalized patch
  "Id4", "Epha4", "Sgk1", "Wnt2", "Gng2", # generalized matrix
  "Necab1", # D1 patch
  "Asic4", # D2 patch
  "Chat", #ChAT neurons
  "Npy", # NPY neurons
  "Sst", "Chodl", # NPY-SST neurons
  "Pthlh", # Pthlh neurons
  "Otof", "Cacng5", "Pcdh8", "Casz1", "Col11a1", "Adarb2",  # eSPN
  # within region
  "Cdh6", # ASt general
  "Itga6", "Arhgef25", # ASt D1
  "Crtac1", "Tesk2", # ASt D2
  "Gabrg1", # CeA D1
  "Pou3f1", # DS D1
  #NONE FOR TS D1
  "Dpyd", # CeA D2
  "C79798", # DS D2
  #NONE FOR TS D2
  "Tnfaip8", # CeA Astro
  "Grip1") # DS Astro

astro_markers <- c("Aqp4", "Gfap", "S100b", "Aldh1l1", "Agt", "Slc1a3", # general markers
                   "Hapln1", "Fam163a", "Lypd6", # main astrocytes
                   "Thbs4", "Epha3", "Adgrv1", # GFAP/reactive astrocytes
                   "Gria1", "Slc6a11", "Lrig1") # Agt astrocytes

micro_markers <- c("Tmem119", "Ptprc", # general markers
                   "Cd83", "Csf1", # dendritic cells
                   "F13a1", "Aoah") # macrophages
