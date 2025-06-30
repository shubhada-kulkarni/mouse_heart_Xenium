library(monocle3)

gene_annotation = as.data.frame(rownames(xenium.obj.female_mi_niche@assays$SCT@counts))
colnames(gene_annotation) = "gene_short_name"
rownames(gene_annotation) = gene_annotation$gene_short_name


cds <- new_cell_data_set(
  xenium.obj.female_mi_niche@assays$SCT@counts,
  cell_metadata = xenium.obj.female_mi_niche@meta.data,
  gene_metadata = gene_annotation
)

cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by="celltype", trajectory_graph_color = "black", group_label_size = 5,
           graph_label_size = 5)
cds <- order_cells(cds) # this step requires interactive selection for ordering the cells.