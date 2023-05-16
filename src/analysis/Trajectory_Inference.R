library(slingshot)
library(Seurat)
library(SeuratDisk)
library(icesTAF)
library(viridis)
library(RColorBrewer)

coch.p8 =  LoadH5Seurat("data/processed/cochlear_p8_Annotated.h5Seurat")
coch.p15 = LoadH5Seurat("data/processed/cochlear_p15_Annotated.h5Seurat")

mkdir('fig/FigureS2')
mkdir('fig/Figure4-S1')

vec_pairwise_sub <- function(c,u){
  return(abs(u-c))
}

# P8 Traj
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
ct_labels = coch.p8[["Annotation"]]
umap = Embeddings(coch.p8[["umap"]])
for(start in c("Deiters/Pillar/Phalangeal Cells", "Glial Cells", "Interdental Cells", "Lateral Greater Epithelial Range (LGER)", "Roof Cells")){
  filt = as.vector(((ct_labels == start) | (ct_labels == "Hair Cells")) & (as.factor(coch.p8$groups) == 'rosa26GAP_p8'))
  lineages <- getLineages(data = Embeddings(coch.p8, reduction = "pca")[filt,1:50],
                          clusterLabels = ct_labels[filt,],
                          end.clus = "Hair Cells",
                          start.clus = start)
  
  
  crv = getCurves(lineages)
  psd_t = slingPseudotime(crv)
  
  ## Fig S2
  gene = "SOX9"
  kern_avg = ksmooth(psd_t, coch.p8[["RNA"]][gene ,as.vector(filt)], kernel = "normal", bandwidth = 10)
  tt <- apply(psd_t,1,vec_pairwise_sub,kern_avg$x)
  colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(400))
  colors <- rev(viridis(400))
  log_col = log(rowSums(tt<2.5)+1)
  
  plotcol <- colors[cut(log_col, breaks=400)]
  png(paste0(paste0('fig/FigureS2/SOX9_', substring(start,1,5)), '_p8.png'), width = 1200, height = 1100)
  par(
    mar      = c(7, 7, 4, 4),
    cex.axis = 3,
    cex.lab  = 3,
    cex.main = 4.5
  )
  plot(kern_avg$x/max(kern_avg$x),  kern_avg$y, col = plotcol, xlab = 'Pseudotime', ylab = paste0("Average ",paste0(gene, ' Expression')), cex = 3, pch = 19)
  legend_image <- as.raster(matrix(rev(colors), ncol=1))
  if( start == "Lateral Greater Epithelial Range (LGER)"){
    title("LGER")
  }
  else{ title(start) }
  dev.off()
  ### Fig S2 Legends
  png(paste0(paste0('fig/FigureS2/SOX9_',substring(start,1,5)),'_p8_legend.png'), width = 600, height =1200)
  par(
    mar      = c(6, 6, 3, 3),
    cex.axis = 2.5,
    cex.lab  = 2.5,
    cex.main = 4
  )
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Log Cell Density')
  text(x=1.5, y = seq(0,1,l=3), labels = seq(round(min(log_col)),round(max(log_col)),l=3), cex = 4)
  rasterImage(legend_image, 0, 0, 1,1)
  dev.off()
  
  ## Fig 4/S1
  colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
  plotcol <- colors[cut(psd_t, breaks=100)]
  png(paste0(paste0('fig/Figure4-S1/Traj_', substr(start,1,5)), '_p8.png'), width = 900, height = 850)
  legend_image <- as.raster(matrix(rev(colors), ncol=1))
  par(
    mar      = c(8, 8, 5, 5),
    cex.axis = 2.5,
    cex.lab  = 2.5,
    cex.main = 3.5
  )
  plot(umap[filt, 1:2], col = plotcol, cex = 2.5, pch = 16)
  title(start)
  
  dev.off()
}

### Fig 4/S1 legend
png('fig/Figure4-S1/Traj_legend.png', width = 900, height = 1200)
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Legend', cex = 3, cex.main =4)
text(x=1.5, y = seq(0,1,l=3), labels = seq(0,1,l=3), cex = 3)
rasterImage(legend_image, 0, 0, 1,1)


# P15 Traj
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
ct_labels = coch.p15[["Annotation"]]
umap = Embeddings(coch.p15[["umap"]])
for(start in c("Deiters/Pillar/Phalangeal Cells", "Glial Cells", "Interdental Cells", "Lateral Greater Epithelial Range (LGER)", "Roof Cells")){
  filt = as.vector(((ct_labels == start) | (ct_labels == "Hair Cells")) & (as.factor(coch.p15$groups) == 'rosa26GAP_p15'))
  lineages <- getLineages(data = Embeddings(coch.p15, reduction = "pca")[filt,1:50],
                          clusterLabels = ct_labels[filt,],
                          end.clus = "Hair Cells",
                          start.clus = start)
  
  
  crv = getCurves(lineages)
  psd_t = slingPseudotime(crv)
  
  ## Fig S2
  gene = "SOX9"
  kern_avg = ksmooth(psd_t, coch.p15[["RNA"]][gene ,as.vector(filt)], kernel = "normal", bandwidth = 10)
  tt <- apply(psd_t,1,vec_pairwise_sub,kern_avg$x)
  colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(400))
  colors <- rev(viridis(400))
  log_col = log(rowSums(tt<2.5)+1)
  
  plotcol <- colors[cut(log_col, breaks=400)]
  png(paste0(paste0('fig/FigureS2/SOX9_',substring(start,1,5)),'_p15.png'), width = 1200, height = 1100)
  par(
    mar      = c(7, 7, 4, 4),
    cex.axis = 3,
    cex.lab  = 3,
    cex.main = 4.5
  )
  plot(kern_avg$x/max(kern_avg$x),  kern_avg$y, col = plotcol, xlab = 'Pseudotime', ylab = paste0("Average ",paste0(gene, ' Expression')), cex = 3, pch = 19)
  legend_image <- as.raster(matrix(rev(colors), ncol=1))
  if( start == "Lateral Greater Epithelial Range (LGER)"){
    title("LGER")
  }
  else{ title(start) }
  dev.off()
  ### Fig S2 Legends
  png(paste0(paste0('fig/FigureS2/SOX9_',substring(start,1,5)),'_p15_legend.png'), width = 600, height =1200)
  par(
    mar      = c(6, 6, 3, 3),
    cex.axis = 2.5,
    cex.lab  = 2.5,
    cex.main = 4
  )
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Log Cell Density')
  text(x=1.5, y = seq(0,1,l=3), labels = seq(round(min(log_col)),round(max(log_col)),l=3), cex = 4)
  rasterImage(legend_image, 0, 0, 1,1)
  dev.off()
  
  ## Fig 4/S1
  colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
  plotcol <- colors[cut(psd_t, breaks=100)]
  png(paste0(paste0('fig/Figure4-S1/Traj_', substr(start,1,5)), '_p15.png'), width = 900, height = 850)
  legend_image <- as.raster(matrix(rev(colors), ncol=1))
  par(
    mar      = c(8, 8, 5, 5),
    cex.axis = 2.5,
    cex.lab  = 2.5,
    cex.main = 3.5
  )
  plot(umap[filt, 1:2], col = plotcol, cex = 2.5, pch = 16)
  title(start)
  dev.off()
}


