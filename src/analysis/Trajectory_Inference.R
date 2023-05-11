library(slingshot)


## P8
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

umap = Embeddings(coch.p8[["umap"]])
for(start in c("Deiters/Pillar/Phalangeal Cells", "Glia Cells", "Interdental Cells", "LGER", "Roof Cells")){
  filt = (ct_labels == start) | (ct_labels == "Hair Cells")
  lineages <- getLineages(data = Embeddings(coch.p8, reduction = "pca")[filt,1:50],
                          clusterLabels = ct_labels[filt],
                          end.clus = "Hair Cells",
                          start.clus = start)
  
  crv = getCurves(lineages)
  psd_t = slingPseudotime(crv)
  
  colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
  plotcol <- colors[cut(psd_t, breaks=100)]
  
  #png(paste0(paste0('../fig/Figure4_P8', substr(start,1,5)), '.png'), width = 1200, height = 1100)
  #layout(matrix(1:2,ncol=2), width = c(2,.75),height = c(1,1))
  legend_image <- as.raster(matrix(rev(colors), ncol=1))
  par(
    mar      = c(6, 6, 3, 3),
    cex.axis = 2.5,
    cex.lab  = 2.5
  )
  plot(umap[filt, 1:2], col = plotcol, cex = 2.5, pch = 16)
  #title(paste0(start," Development"))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Legend')
  text(x=1.5, y = seq(0,1,l=3), labels = seq(0,1,l=3))
  rasterImage(legend_image, 0, 0, 1,1)
  dev.off()
  
}



## P15
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

umap = Embeddings(coch.p15[["umap"]])
for(start in c("Deiters/Pillar/Phalangeal Cells", "Glia Cells", "Interdental Cells", "LGER", "OC90+ Cells")){
  filt = (ct_labels == start) | (ct_labels == "Hair Cells")
  lineages <- getLineages(data = Embeddings(coch.p15, reduction = "pca")[filt,1:50],
                          clusterLabels = ct_labels[filt],
                          end.clus = "Hair Cells",
                          start.clus = start)

  crv = getCurves(lineages)
  psd_t = slingPseudotime(crv)
  
  colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
  plotcol <- colors[cut(psd_t, breaks=100)]
  
  #png(paste0(paste0('../fig/Figure4_P15', substr(start,1,5)), '.png'), width = 1200, height = 1100)
  #layout(matrix(1:2,ncol=2), width = c(2,.75),height = c(1,1))
  legend_image <- as.raster(matrix(rev(colors), ncol=1))
  par(
    mar      = c(6, 6, 3, 3),
    cex.axis = 2.5,
    cex.lab  = 2.5,
    cex.main = 4
  )
  plot(umap[filt, 1:2], col = plotcol, cex = 2.5, pch = 16)
  title(paste0(start," Development"))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pseudotime')
  text(x=1.5, y = seq(0,1,l=3), labels = seq(0,1,l=3), cex = 2.5)
  rasterImage(legend_image, 0, 0, 1,1)
  dev.off()
}