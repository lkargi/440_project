vec_pairwise_sub <- function(c,u){
  return(abs(u-c))
}

## P8
start = "Deiters/Pillar/Phalangeal Cells"
ct_labels = coch.p8[["Annotation"]]
filt = (ct_labels == start) | (ct_labels == "Hair Cells")
lineages <- getLineages(data = Embeddings(coch.p8, reduction = "pca")[filt,1:50],
                        clusterLabels = ct_labels[filt],
                        end.clus = "Hair Cells",
                        start.clus = start)

crv = getCurves(lineages)
psd_t = slingPseudotime(crv)
psd_t <- as.matrix(psd_t)

kern_avg = ksmooth(psd_t, coch.p8[["RNA"]]["NOTCH1",as.vector(filt)], kernel = "normal", bandwidth = 5)
tt <- apply(psd_t,1,vec_pairwise_sub,kern_avg$x)
colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(400))
colors <- rev(viridis(400))
log_col = log(rowSums(tt<1.25))

plotcol <- colors[cut(log_col, breaks=400)]
layout(matrix(1:2,ncol=2), width = c(2,.75),height = c(1,1))
plot(kern_avg$x/max(kern_avg$x),  kern_avg$y, col = plotcol, xlab = 'Pseudotime', ylab = 'SOX9 Expression')
legend_image <- as.raster(matrix(rev(colors), ncol=1))
title("Moving Average Expression")
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Log Density', cex.main=.75)
text(x=1.5, y = seq(0,1,l=3), labels = seq(round(min(log_col)),round(max(log_col)),l=3))
rasterImage(legend_image, 0, 0, 1,1)



## P15
start = "Deiters/Pillar/Phalangeal Cells"
ct_labels = coch.p15[["Annotation"]]
filt = (ct_labels == start) | (ct_labels == "Hair Cells")
lineages <- getLineages(data = Embeddings(coch.p15, reduction = "pca")[filt,1:50],
                        clusterLabels = ct_labels[filt],
                        end.clus = "Hair Cells",
                        start.clus = start)

crv = getCurves(lineages)
psd_t = slingPseudotime(crv)
psd_t <- as.matrix(psd_t)

kern_avg = ksmooth(psd_t, coch.p15[["RNA"]]["NOTCH1",as.vector(filt)], kernel = "normal", bandwidth = 5)
tt <- apply(psd_t,1,vec_pairwise_sub,kern_avg$x)
colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(400))
colors <- rev(viridis(400))
log_col = log(rowSums(tt<1.25))

plotcol <- colors[cut(log_col, breaks=400)]
layout(matrix(1:2,ncol=2), width = c(2,.75),height = c(1,1))
plot(kern_avg$x/max(kern_avg$x),  kern_avg$y, col = plotcol, xlab = 'Pseudotime', ylab = 'SOX9 Expression')
legend_image <- as.raster(matrix(rev(colors), ncol=1))
title("Moving Average Expression")
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Log Density', cex.main=.75)
text(x=1.5, y = seq(0,1,l=3), labels = seq(round(min(log_col)),round(max(log_col)),l=3))
rasterImage(legend_image, 0, 0, 1,1)
