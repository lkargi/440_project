ct_labels = as.vector(coch.p8$seurat_clusters)
ct_labels[] = "Unclassified"
ct_labels[(coch.p8$seurat_clusters == 0) | (coch.p8$seurat_clusters == 1) | (coch.p8$seurat_clusters == 2) | (coch.p8$seurat_clusters == 3) ] = "Deiters/Pillar/Phalangeal Cells"

ct_labels[(coch.p8$seurat_clusters == 8) | (coch.p8$seurat_clusters == 10)] = "Roof Cells"
ct_labels[(coch.p8$seurat_clusters == 6) | (coch.p8$seurat_clusters == 14)] = "Hair Cells"
ct_labels[(coch.p8$seurat_clusters == 4) | (coch.p8$seurat_clusters == 9)] = "Interdental Cells"
#ct_labels[(coch.p8$seurat_clusters == 7)] = "IPhC"
ct_labels[(coch.p8$seurat_clusters == 17) | (coch.p8$seurat_clusters == 19)] = "Glia Cells"
ct_labels[(coch.p8$seurat_clusters == 5) | (coch.p8$seurat_clusters == 15) | (coch.p8$seurat_clusters == 16)] = "Lateral Greater Epithelial Range (LGER)"
#ct_labels[(coch.p8$seurat_clusters == 18)] = "Claudius"


coch.p8[["Annotation"]] = ct_labels
DimPlot(coch.p8, reduction = 'umap',group.by = "Annotation", cols = c("red", "green", "blue", "purple", 
                                                                           "cyan", "orange","gray"))
ggsave(
 "coch-p8-umap.png",
 width = 7.5,
 height = 5,
 dpi = 1200
)

ct_labels = as.vector(coch.p15$seurat_clusters)
ct_labels[] = "Unclassified"
ct_labels[(coch.p15$seurat_clusters == 7) | (coch.p15$seurat_clusters == 4) | (coch.p15$seurat_clusters == 2)] = "Deiters/Pillar/Phalangeal Cells"
ct_labels[(coch.p15$seurat_clusters == 14)] = "Roof Cells"
ct_labels[(coch.p15$seurat_clusters == 11)] = "Hair Cells"
ct_labels[(coch.p15$seurat_clusters == 0)] = "Interdental Cells"
ct_labels[(coch.p15$seurat_clusters == 6)] = "Glia Cells"
ct_labels[(coch.p15$seurat_clusters == 5) | (coch.p15$seurat_clusters == 8) | (coch.p15$seurat_clusters == 17)] = "LGER"



coch.p15[["Annotation"]] = ct_labels
DimPlot(coch.p15, reduction = 'umap',group.by = "Annotation", cols = c("red", "green", "blue", "purple", 
                                                                            "cyan", "orange","gray"))
ggsave(
  "coch-p15-umap.png",
  width = 7.5,
  height = 5,
  dpi = 1200
)