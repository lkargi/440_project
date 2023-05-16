library(Seurat)
library(SeuratDisk)
library(SeuratData)

coch.p8 =  LoadH5Seurat("data/processed/cochlear_p8.h5Seurat")
coch.p15 = LoadH5Seurat("data/processed/cochlear_p15.h5Seurat")

#ANNOTATE CLUSTERS
ct_labels = as.vector(coch.p8$seurat_clusters)
ct_labels[] = "Unclassified"
ct_labels[(coch.p8$seurat_clusters == 0) | (coch.p8$seurat_clusters == 1) | (coch.p8$seurat_clusters == 2) | (coch.p8$seurat_clusters == 3) ] = "Deiters/Pillar/Phalangeal Cells"
ct_labels[(coch.p8$seurat_clusters == 8) | (coch.p8$seurat_clusters == 10)] = "Roof Cells"
ct_labels[(coch.p8$seurat_clusters == 6) | (coch.p8$seurat_clusters == 14)] = "Hair Cells"
ct_labels[(coch.p8$seurat_clusters == 4) | (coch.p8$seurat_clusters == 9)] = "Interdental Cells"
ct_labels[(coch.p8$seurat_clusters == 17) | (coch.p8$seurat_clusters == 19)] = "Glia Cells"
ct_labels[(coch.p8$seurat_clusters == 5) | (coch.p8$seurat_clusters == 15) | (coch.p8$seurat_clusters == 16)] = "Lateral Greater Epithelial Range (LGER)"
coch.p8[["Annotation"]] = ct_labels

ct_labels = as.vector(coch.p15$seurat_clusters)
ct_labels[] = "Unclassified"
ct_labels[(coch.p15$seurat_clusters == 7) | (coch.p15$seurat_clusters == 4) | (coch.p15$seurat_clusters == 2)] = "Deiters/Pillar/Phalangeal Cells"
ct_labels[(coch.p15$seurat_clusters == 14)] = "Roof Cells"
ct_labels[(coch.p15$seurat_clusters == 11)] = "Hair Cells"
ct_labels[(coch.p15$seurat_clusters == 0)] = "Interdental Cells"
ct_labels[(coch.p15$seurat_clusters == 6)] = "Glia Cells"
ct_labels[(coch.p15$seurat_clusters == 5) | (coch.p15$seurat_clusters == 8) | (coch.p15$seurat_clusters == 17)] = "Lateral Greater Epithelial Range (LGER)"
coch.p15[["Annotation"]] = ct_labels

#RESAVE ANNOTATED CLUSTERS
SaveH5Seurat(coch.p8, "data/processed/cochlear_p8_Annotated" ,overwrite = TRUE)
SaveH5Seurat(coch.p15, "data/processed/cochlear_p15_Annotated" ,overwrite = TRUE)