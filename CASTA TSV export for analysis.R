##Outputs TSVs of the whole ST dataset, as well as subsets for the epithelium and stroma
##Also makes a TSV for the alignment values
##These datasets need to be merged in python

#libraries
install.packages("Seurat")
#install.packages("SeuratData")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("dplyr")
install.packages("devtools")
install.packages("hdf5r")
install.packages('BiocManager')
BiocManager::install('glmGamPoi')
BiocManager::install('ggseurat')
remotes::install_github('mojaveazure/ggseurat')
install.packages("ggplot2")
install.packages("raster")

library(ggseurat)
library(glmGamPoi)
library(hdf5r)
library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(raster)

setwd("/Users/samellis/Library/CloudStorage/OneDrive-VUMC/Collagen_Alignment_Project/Refined analysis pipeline/TR11_18105/TR11_18105_ST")

#load data
tumor <- Load10X_Spatial(
  "outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  bin.size = NULL,
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL
)

#perform SCT transformation
tumor <- SCTransform(tumor, assay = "Spatial", verbose = FALSE)

#optionally visualize gene expression
SpatialFeaturePlot(tumor, features = c("FAP", "ACTA2", "VIM", "DDR1"), pt.size = 3 )

#run PCA and generate UMAP. This is needed to subset the spots
tumor <- RunPCA(tumor, assay = "SCT", verbose = FALSE)
tumor <- FindNeighbors(tumor, reduction = "pca", dims = 1:30)
tumor <- FindClusters(tumor, verbose = FALSE)
tumor <- RunUMAP(tumor, reduction = "pca", dims = 1:30)

#optional: add in UMAP plotting

#make a dataframe with the expression data
tumor.df <- as.data.frame(tumor[["SCT"]]$counts)
tumor.df <- t(tumor.df)

#save it
setwd("/Users/samellis/Library/CloudStorage/OneDrive-VUMC/Collagen_Alignment_Project/Refined analysis pipeline/TR11_18105/TR11_18105_ST")
write.table(tumor.df,  
            file= paste("TR11_18105_full_dataframe.tsv"), 
            sep='\t', 
            col.names = T
)

##Now subset the area to stroma and epithelium

SpatialDimPlot(tumor, interactive = TRUE)
#subset out just the stroma
#change "idents" to the numbers of the regions that are stromal
#for TR11_206: 2,3,4,5,8
#for TR11_16184: 2,4,5,7
#for TR11_18105: 0,1,2,3,4,6
#for TR11_21723: 0,1,4,6,8,9
#for TR11_23542: 0,1,2,4,7,8
#need to hit done before proceeding
clusters1 <- c(0,1,2,3,4,6)
stroma <- subset(tumor, idents = clusters1)

#make a stroma dataframe
stroma.df <- as.data.frame(stroma[["SCT"]]$counts)
stroma.df <- as.data.frame(t(stroma.df))
stroma.df$barcodes <- rownames(stroma.df)

#save it
setwd("/Users/samellis/Library/CloudStorage/OneDrive-VUMC/Collagen_Alignment_Project/Refined analysis pipeline/TR11_18105/TR11_18105_ST")
write.table(stroma.df,  
            file= paste("TR11_18105_stroma.tsv"), 
            sep='\t', 
            col.names = T
)

#make a df with the epithelium
#change "idents" to the stromal clusters
#for TR11_206: 0,7
#for TR11_16184: 0,8,9,11,13
#for TR11_18105: 5, 9
#for TR11_21723: 3,2,7
#for TR11_23542: 5,3,6
epithelium <- subset(tumor, idents = c(5, 9))
epithelium.df <- as.data.frame(epithelium[["SCT"]]$counts)
epithelium.df <- as.data.frame(t(epithelium.df))
epithelium.df$barcodes <- rownames(epithelium.df)

#save it
setwd("/Users/samellis/Library/CloudStorage/OneDrive-VUMC/Collagen_Alignment_Project/Refined analysis pipeline/TR11_18105/TR11_18105_ST")
write.table(epithelium.df,  
            file= paste("TR11_18105_epithelium.tsv"), 
            sep='\t', 
            col.names = T
)


##Do the same for the heatmap

setwd("/Users/samellis/Library/CloudStorage/OneDrive-VUMC/Collagen_Alignment_Project/Refined analysis pipeline/TR11_18105/Transformed_Heatmaps")

#open image
heatmap_raster <- raster("heat_aligned.tif")
#heatmap_matrix <- rasterToPoints(heatmap_raster, progress="text")

#load in the locations of the barcodes
#this is the "tissue positions" csv
setwd("/Users/samellis/Library/CloudStorage/OneDrive-VUMC/Collagen_Alignment_Project/Refined analysis pipeline/TR11_18105/TR11_18105_ST/outs/spatial")
bar_pos <- as.data.frame(read.csv("tissue_positions.csv"))

#make pixel locations integers
bar_pos <- cbind(bar_pos, as.integer(bar_pos$pxl_row_in_fullres))
bar_pos <- cbind(bar_pos, as.integer(bar_pos$pxl_col_in_fullres))

#rename the columns to reflect the change
colnames(bar_pos)[7] <-  "pxl_row_in_fullres_int"
colnames(bar_pos)[8] <-  "pxl_col_in_fullres_int"

#make a new dataframe to hold the alignment values
heatmap_df <- bar_pos
heatmap_df$alignment <- NA

#extract pixel values from the image at the locations of the spots and add to the matrix
#get scale factor from outs/spatial/scale_factors.json
#for TR11_2172 and TR16_23542 : 0.08794
#for TR11_18105: 0.459155 (not sure why. That's not what's reported in the scale factor file)

scale_factor <- 0.459155

for (i in 1:length(heatmap_df[,1])){
  row_pix <- as.integer((heatmap_df$pxl_row_in_fullres_int[i] * scale_factor))
  col_pix <- as.integer((heatmap_df$pxl_col_in_fullres_int[i] * scale_factor) )
  align <- extract(heatmap_raster, SpatialPoints(cbind(col_pix, row_pix)))
  #if (align != 0){
  heatmap_df$alignment[i] <- align
  #}
  #else {
  # heatmap_df$alignment[i] <- NA
  #}
}

colnames(heatmap_df)[1] <- "barcodes"

#plot the alignment locations to check that they look right
align_plot <- ggplot(data = heatmap_df, 
                     aes(y = pxl_row_in_fullres, 
                         x = pxl_col_in_fullres, 
                         color = alignment
                     )
)
align_plot <- align_plot + geom_point()
align_plot

#save it
setwd("/Users/samellis/Library/CloudStorage/OneDrive-VUMC/Collagen_Alignment_Project/Refined analysis pipeline/TR11_18105/TR11_18105_ST")
write.table(heatmap_df,  
            file= paste("TR11_18105_alignment_vals.tsv"), 
            sep='\t', 
            col.names = T
)

