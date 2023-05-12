library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)

coch.p8 =  LoadH5Seurat("data/processed/cochlear_p8_Annotated.h5Seurat")
coch.p15 = LoadH5Seurat("data/processed/cochlear_p15_Annotated.h5Seurat")

DimPlot(coch.p8, reduction = 'umap',group.by = "Annotation", cols = c("red", "green", "blue", "purple", 
                                                                           "cyan", "orange","gray"))
ggsave(
 "fig/coch-p8-umap.png",
 width = 7.5,
 height = 5,
 dpi = 1200
)


DimPlot(coch.p15, reduction = 'umap',group.by = "Annotation", cols = c("red", "green", "blue", "purple", 
                                                                            "cyan", "orange","gray"))
ggsave(
  "fig/coch-p15-umap.png",
  width = 7.5,
  height = 5,
  dpi = 1200
)