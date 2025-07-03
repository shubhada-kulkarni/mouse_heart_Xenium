library(monocle3)
library(RColorBrewer)
library(Seurat)

xenium.obj.female_mi$celltype <- Idents(xenium.obj.female_mi)
data <- as(as.matrix(xenium.obj.female_mi@assays$SCT@counts), 'sparseMatrix')
cell_metadata <- data.frame(xenium.obj.female_mi@meta.data)
gene_annotation <- data.frame(rownames(data))
colnames(gene_annotation) <- "gene_short_name"
rownames(gene_annotation) <- gene_annotation$gene_short_name
cds <- new_cell_data_set(data,cell_metadata  = cell_metadata,gene_metadata  = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds,preprocess_method = 'PCA',reduction_method = c("UMAP"),cores = 0)
cds@int_colData@listData$reducedDims$UMAP <- xenium.obj.female_mi@reductions$umap@cell.embeddings
cds@clusters$UMAP$clusters <- xenium.obj.female_mi@active.ident
cds <- cluster_cells(cds)

#color by seurat cluster
plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "celltype")

#color by monocle cluster
plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster=FALSE)
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE, 
           label_leaves=F, label_branch_points=F, cell_size = 1.5)
cds <- order_cells(cds, reduction_method = "UMAP")
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE, 
           label_leaves=F, label_branch_points=F)

# Plotting genes in pseudotime
plot_genes_in_pseudotime(cds["Ankrd1-FSJ2",], color_cells_by="celltype", 
                              min_expr=0.5, ncol = 2)