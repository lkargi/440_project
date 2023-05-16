library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(icesTAF)

coch.p8 =  LoadH5Seurat("data/processed/cochlear_p8_Annotated.h5Seurat")
coch.p15 = LoadH5Seurat("data/processed/cochlear_p15_Annotated.h5Seurat")

mkdir('fig/Figure3')
mkdir('fig/FigureS4')

# P8
## Fig 3
pdf(NULL)
DimPlot(coch.p8, reduction = 'umap',group.by = "Annotation", cols = c("red", "green", "blue", "purple", 
                                                                           "cyan", "orange","gray"))
ggsave(
 "fig/Figure3/umap_p8.png",
 width = 7.5,
 height = 5,
 dpi = 1200
)
dev.off()

## Fig S4
pdf(NULL)
FeaturePlot(coch.p8, features = "SLC1A3") + ggtitle("GLAST")
ggsave(
  "fig/FigureS4/GLAST_p8.png",
  width = 6,
  height = 4,
  dpi = 1200
)
dev.off()


# P15
## Fig 3
pdf(NULL)
DimPlot(coch.p15, reduction = 'umap',group.by = "Annotation", cols = c("red", "green", "blue", "purple", 
                                                                            "cyan", "orange","gray"))
ggsave(
  "fig/Figure3/umap_p15.png",
  width = 7.5,
  height = 5,
  dpi = 1200
)
dev.off()

## Fig S4
pdf(NULL)
FeaturePlot(coch.p15, features = "SLC1A3") + ggtitle("GLAST")
ggsave(
  "fig/FigureS4/GLAST_p15.png",
  width = 6,
  height = 4,
  dpi = 1200
)
dev.off()

