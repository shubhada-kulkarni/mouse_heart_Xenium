# R script to re-process the qupath_mod data which is the same data by removing two bad outlier cells.
# The script was written in order to get cellnames for Xenium explorer.

xenium.obj.female_reprocessed <- LoadXenium(data.dir = "/prj/XeniumProbeDesign/mouse_heart_03012025/output-XETG00046__0046388__Region_1__20241219__110827/qupath_seg_mod/outs/")
# filter to avoid SCtransform error
xenium.obj.female_reprocessed <- subset(xenium.obj.female_reprocessed, subset = nCount_Xenium > 0) # this filtering needs to be done or there is an error

options(future.globals.maxSize = 8000 * 1024^2) # to avoid SCTrasnform getGlobalsAndPackages error
xenium.obj.female_reprocessed <- SCTransform(xenium.obj.female_reprocessed, assay = "Xenium")
xenium.obj.female_reprocessed <- RunPCA(xenium.obj.female_reprocessed, npcs = 30, features = rownames(xenium.obj.female_reprocessed))
xenium.obj.female_reprocessed <- RunUMAP(xenium.obj.female_reprocessed, dims = 1:30)
xenium.obj.female_reprocessed <- FindNeighbors(xenium.obj.female_reprocessed, reduction = "pca", dims = 1:30)
xenium.obj.female_reprocessed <- FindClusters(xenium.obj.female_reprocessed, resolution = 0.5)

# Remove the cells that have low feature counts
xenium.obj.female_reprocessed <- subset(xenium.obj.female_reprocessed, subset = nFeature_Xenium >= 25)

## Marker gene analysis to assign clusters
#Cardiomyocytes
VlnPlot(xenium.obj.female_reprocessed, c("Ttn"), pt.size = 0)
# Fibroblasts - 1, 4, 10
FeaturePlot(xenium.obj.female_reprocessed, c("Dpt", "Lum"), label = T)
# Scahwann cells - 15
FeaturePlot(xenium.obj.female_reprocessed, "Plp1", label = T)
 # EC cells - Coronary artery EC cells - 2, 11 - using reverse engineering with Tabula Muris expression database
FeaturePlot(xenium.obj.female_reprocessed, c("Eng1", "Tie1"), label = T)
# Lymphatic ECs - 14
FeaturePlot(xenium.obj.female_reprocessed, c("Lyve1"), label = T)
# Erythrocytes & pericytes - reverse engineering - 9
FeaturePlot(xenium.obj.female_reprocessed, c("Higd1b", "Rsad2"), label = T)
# Endocardial - 11
FeaturePlot(xenium.obj.female_reprocessed, c("Vwf", "Plvap"), label = T)
# Epicardium
FeaturePlot(xenium.obj.female_reprocessed, c("Upk3b"), label = T)
# SMCs
FeaturePlot(xenium.obj.female_reprocessed, c("Myh11"), label = T)
# Muscle contraction
FeaturePlot(xenium.obj.female_reprocessed, c("Cnn1"), label = T)

backup <- xenium.obj.female_reprocessed
xenium.obj.female_reprocessed <- backup

## Rename clusters
new.cluster.ids.1 <- c("ECs", "Fibroblasts", "CMs", "CMs", "Fibroblasts", "Fibroblasts", "ECs", "ECs",
                       "MuscleContraction", "Erythrocytes", "Fibroblasts", "Endocardial", "CMs",
                       "SMCc", "LymphaticEC", "Neural", "CMs", "Epicardial")

names(new.cluster.ids.1) <- levels(xenium.obj.female_reprocessed)
xenium.obj.female_reprocessed <- RenameIdents(xenium.obj.female_reprocessed, new.cluster.ids.1)

## Save the cell ID information to upload in Xenium Explorer

clusters <- data.frame(names(Idents(xenium.obj.female_reprocessed)), as.vector(Idents(xenium.obj.female_reprocessed)))
colnames(clusters) <- c("cell_id", "group")
write.csv(clusters, "../r_objects/celltype_xenium_female_reprocessed.csv", quote = F, row.names = F)


