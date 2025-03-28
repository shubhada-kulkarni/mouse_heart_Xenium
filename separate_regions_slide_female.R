# R script used to read in the readxenium object for combined regions with female heart samples and 
# separate these regions into four different Seurat objects

female <- readRDS("/prj/XeniumProbeDesign/mouse_heart_03012025/r_objects/mouse_female_fourregions.rds")

# wild-type cell-id extraction
count <- 0
count_na <- 0
cellids_wt <- c()
for (i in 1:nrow(female$centroids)){
  # WT
  if ( (female$centroids[i,"x"] > 500) & (female$centroids[i,"x"] < 5600) & 
       (female$centroids[i,"y"] > 600) & (female$centroids[i,"y"] < 6600)){
    #print(paste("This is wild-type sample", female$centroids[i,"x"], female$centroids[i,"y"], female$centroids[i,"cell"]))
    count <- count + 1
    cellids_wt <- c(cellids_wt, female$centroids[i,"cell"])
  }
  # MI
  else if ( (female$centroids[i,"x"] > 5500) & (female$centroids[i,"x"] < 11000) &
            (female$centroids[i,"y"] > 6100) & (female$centroids[i,"y"] < 12100)){ 
    #print(paste("This is MI sample", female$centroids[i,"x"], female$centroids[i,"y"], female$centroids[i,"cell"])) 
    cell_IDs_samples[i,1:2] <-  c(female$centroids[i,"cell"], "MI") 
    count <- count + 1 
    #cellids_wt <- c(cellids_wt, female$centroids[i,"cell"]) 
  } 
  # ORAB1
  else if ( (female$centroids[i,"x"] > 400) & (female$centroids[i,"x"] < 5200) &
            (female$centroids[i,"y"] > 14500) & (female$centroids[i,"y"] < 20200)){ 
    #print(paste("This is ORAB1 sample", female$centroids[i,"x"], female$centroids[i,"y"], female$centroids[i,"cell"])) 
    cell_IDs_samples[i,1:2] <-  c(female$centroids[i,"cell"], "ORAB1") 
    count <- count + 1 
    #cellids_wt <- c(cellids_wt, female$centroids[i,"cell"]) 
  } 
  # ORAB2 
  else if ( (female$centroids[i,"x"] > 6200) & (female$centroids[i,"x"] < 11000) &
            (female$centroids[i,"y"] > 15800) & (female$centroids[i,"y"] < 21350)){ 
    #print(paste("This is ORAB2 sample", female$centroids[i,"x"], female$centroids[i,"y"], female$centroids[i,"cell"])) 
    cell_IDs_samples[i,1:2] <-  c(female$centroids[i,"cell"], "ORAB2") 
    count <- count + 1 
    #cellids_wt <- c(cellids_wt, female$centroids[i,"cell"]) 
  } 
  else {
    print(paste("Cell not assigned - ", female$centroids[i,"x"], female$centroids[i,"y"], female$centroids[i,"cell"])) 
    count_na <- count_na + 1
  }
}
colnames(cell_IDs_samples) <- c("cellid", "sample") 
print(table(cell_IDs_samples)) 
colnames(cell_IDs_samples) <- c("cellid", "condition") 
cell_IDs_samples <- cell_IDs_samples[!is.na(cell_IDs_samples$cellid),] 

print(count)

cm_female_wt <- female$matrix$`Gene Expression`[,cellids_wt]