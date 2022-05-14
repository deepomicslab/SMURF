library(splatter)
library(pheatmap)
library(Seurat)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(viridis)

plot_doheatmap <- function(m, cluster, markers=NULL){
  s_object <- CreateSeuratObject(m)
  # s_object <- NormalizeData(s_object, display.progress = FALSE)
  s_object <- FindVariableFeatures(s_object)
  s_object <- ScaleData(s_object, display.progress = FALSE)
  s_object <- RunPCA(s_object, do.print = FALSE)
  s_object <- FindNeighbors(s_object)
  s_object <- FindClusters(s_object)
  s_object@meta.data$seurat_clusters <- cluster
  # Set cell identity to sample identity
  s_object <- SetIdent(object = s_object, value = "seurat_clusters")
  s_object <- RunUMAP(object = s_object, dims = 1:5)
  s_object <- RunTSNE(object = s_object, dims = 1:5, check_duplicates = FALSE)
  if(is.null(markers)){
    markers <- FindAllMarkers(s_object, only.pos = TRUE)  
  }
  doheatmap_fig <- DoHeatmap(s_object, markers$gene, slot='scale.data', label=F, draw.lines = F, 
                             group.colors = get_palette('nejm', 5)
  ) + theme(axis.text.y = element_blank()) + 
    scale_fill_gradientn(colours = rev(magma(5)))
  umap_fig <- DimPlot(s_object, reduction='umap', cols = get_palette('nejm', length(unique(cluster)))) +
    theme_bw() + 
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank()
    )
  tsne_fig <- DimPlot(s_object, reduction='tsne', cols = get_palette('nejm', length(unique(cluster)))) +
    theme_bw() + 
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank()
    )  
  return (list(
    doheatmap_fig=doheatmap_fig,
    umap_fig=umap_fig,
    tsne_fig=tsne_fig,
    markers=markers,
    s_object=s_object
  ))
}



fig_list <- list()
metric_df <- NULL

for (dataset in c(
  # 'p3cl'
  # 'sc_10x_5cl'
   # 'sc_10x'
  # 'sc_celseq2_5cl_p1'
  # 'sc_celseq2_5cl_p2'
  # 'sc_celseq2_5cl_p3'
   # 'sc_celseq2'
  # 'sc_dropseq'
  "3Line-qPCR"
)){
  # prefix <- paste0('cellbench/', dataset)
  prefix <- paste0('cellcircle/', dataset)
  prefix
  
  cell_line <- 'cell_line'
  cell_line <- 'Label'
  if (dataset %in% c('sc_celseq2_5cl_p1', 'sc_celseq2_5cl_p2', 'sc_celseq2_5cl_p3')){
    cell_line <- 'cell_line_demuxlet'
  }
  cluster <- read.csv(paste0(prefix, '_metadata.csv'))[[cell_line]]
  
  fig_list[[prefix]] <- list()
  fig_list[[prefix]]$umap <- list()
  fig_list[[prefix]]$doheatmap <- list()
  for (tool in c(
    # 'SMURF_CV'
     # 'SMURF_F'
     # 'SMURF_V'
    "SC_SMURF_F"
  )){
    # m_fn <- paste0(prefix, '_rightMatrix_', tool, '.csv')  
    m_fn <- paste0(prefix, "_",tool, '_H.csv')  
    m_fn
    m <- read.csv(m_fn, row.names = 1)
    # m <- m + abs(min(m))
    
    tool <- paste0(tool, '_H')
    
    res <- plot_doheatmap(m, cluster)  
    fig_list[[prefix]]$doheatmap[[tool]] <- res$doheatmap_fig
    fig_list[[prefix]]$umap[[tool]] <- res$umap_fig
    fig_list[[prefix]]$tsne[[tool]] <- res$tsne_fig
  }
}

for (dataset in c(
   # "p3cl"
  # 'sc_10x_5cl'
    # 'sc_10x'
    # 'sc_celseq2_5cl_p1'
  # 'sc_celseq2_5cl_p2'
   # 'sc_celseq2_5cl_p3'
   # 'sc_celseq2'
   # 'sc_dropseq'
  "3Line-qPCR"
)){
  # prefix <- paste0('cellbench/', dataset)
  prefix <- paste0('cellcircle/', dataset)
  prefix
  
  # cell_line <- 'cell_line'
  cell_line <- 'Label'
  if (dataset %in% c('sc_celseq2_5cl_p1', 'sc_celseq2_5cl_p2', 'sc_celseq2_5cl_p3')){
    cell_line <- 'cell_line_demuxlet'
  }
  cluster <- read.csv(paste0(prefix, '_metadata.csv'))[[cell_line]]
  
  
  for (tool in c(# 'count',
    # 'samp',
    #'C-Impute',
    # 'I-Impute',
    # 'SAVER',
    # 'MAGIC',
    # 'kNN-smoothing'
    # 'SMURF_CV_H'
    # 'SMURF_F_H'
    # 'SMURF_V_H'
    "SC_SMURF_F_H"

  )){
    umap_data <- fig_list[[prefix]]$umap[[tool]]$data[, c('UMAP_1', 'UMAP_2')]
    write.csv(umap_data, paste0(prefix, '_', tool, '_umap.csv'), quote = F)
    #eva_cluster(umap_data, cluster, length(unique(cluster)))
  }
}
