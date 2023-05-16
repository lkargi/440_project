library(slingshot)
library(Seurat)
library(SeuratDisk)
library(stats)
library(icesTAF)

coch.p8 =  LoadH5Seurat("data/processed/cochlear_p8_Annotated.h5Seurat")
coch.p15 = LoadH5Seurat("data/processed/cochlear_p15_Annotated.h5Seurat")

mkdir('fig/Figure5')
mkdir('fig/FigureS3')

vec_pairwise_sub <- function(c,u){
  return(abs(u-c))
}

# P8
start = "Deiters/Pillar/Phalangeal Cells"
ct_labels = coch.p8[["Annotation"]]
filt = as.vector(((ct_labels == start) | (ct_labels == "Hair Cells")) &(as.factor(coch.p8$groups) == 'rosa26A_p8'))
lineages <- getLineages(data = Embeddings(coch.p8, reduction = "pca")[filt,1:50],
                        clusterLabels = ct_labels[filt,],
                        end.clus = "Hair Cells",
                        start.clus = start)

crv = getCurves(lineages)
psd_t = slingPseudotime(crv)
psd_t <- as.matrix(psd_t)

## Figure 5A
png('fig/Figure5/Figure5A_p8.png', width = 900, height = 850)
par(
  mar      = c(8, 8, 5, 5),
  cex.axis = 2.5,
  cex.lab  = 2.5,
  cex.main = 3.5
)
plot(NULL, NULL,xlim = c(0,1), ylim = c(0,2.5), ylab = "Expression", xlab = "Pseudotime" )
colors = c("black", "red", "blue")
color_num = 1
for(gene in c("ATOH1", "GFI1", "POU4F3")){
  kern_avg = ksmooth(psd_t, coch.p8[["RNA"]][gene ,as.vector(filt)], kernel = "normal", bandwidth = 10)
  lines(kern_avg$x/max(kern_avg$x),  kern_avg$y, xlab = 'Pseudotime', ylab = paste0(gene, ' Expression'), cex = 3, col = colors[color_num], lwd = 4.0) 
  color_num= color_num+1
}
legend('topleft', legend = c("ATOH1", "GFI1", "POU4F3"), col = colors, lty = 1, cex = 2.5, lwd = 4)
dev.off()

## Figure 5B
png('fig/Figure5/Figure5B_p8.png', width = 900, height = 850)
par(
  mar      = c(8, 8, 5, 5),
  cex.axis = 2.5,
  cex.lab  = 2.5,
  cex.main = 3.5
)
plot(NULL, NULL,xlim = c(0,1), ylim = c(0,1.5), ylab = "NOTCH1 Expression", xlab = "Pseudotime" )
colors = c("black")
color_num = 1
for(gene in c("NOTCH1")){
  kern_avg = ksmooth(psd_t, coch.p8[["RNA"]][gene ,as.vector(filt)], kernel = "normal", bandwidth = 10)
  lines(kern_avg$x/max(kern_avg$x),  kern_avg$y, xlab = 'Pseudotime', ylab = paste0(gene, ' Expression'), cex = 3, col = colors[color_num],lwd = 4.0) 
  color_num= color_num+1
}
dev.off()

## Figure S3
png('fig/FigureS3/Density_p8.png', width = 900, height = 850)
par(
  mar      = c(6, 6, 3, 3),
  cex.axis = 2.5,
  cex.lab  = 2.5,
  cex.main = 4
)
tt <- apply(psd_t,1,vec_pairwise_sub,kern_avg$x)
log_col = log(rowSums(tt<2.5))
plot(kern_avg$x/max(kern_avg$x), log_col, xlab = "Pseudotime", ylab = "Log Cell Density", type = "l", lwd = 2, ylim = c(0,7.5))


#P15
start = "Deiters/Pillar/Phalangeal Cells"
ct_labels = coch.p15[["Annotation"]]
filt = as.vector(((ct_labels == start) | (ct_labels == "Hair Cells")) &(as.factor(coch.p15$groups) == 'rosa26GAP_p15'))
lineages <- getLineages(data = Embeddings(coch.p15, reduction = "pca")[filt,1:50],
                        clusterLabels = ct_labels[filt,],
                        end.clus = "Hair Cells",
                        start.clus = start)
crv = getCurves(lineages)
psd_t = slingPseudotime(crv)
psd_t <- as.matrix(psd_t)


## Fig 5A
png('fig/Figure5/Figure5A_p15.png', width = 900, height = 850)
par(
  mar      = c(8, 8, 5, 5),
  cex.axis = 2.5,
  cex.lab  = 2.5,
  cex.main = 3.5
)
plot(NULL, NULL,xlim = c(0,1), ylim = c(0,2.5), ylab = "Expression", xlab = "Pseudotime" )
colors = c("black", "red", "blue")
color_num = 1
for(gene in c("ATOH1", "GFI1", "POU4F3")){
  kern_avg = ksmooth(psd_t, coch.p15[["RNA"]][gene ,as.vector(filt)], kernel = "normal", bandwidth = 10)
  lines(kern_avg$x/max(kern_avg$x),  kern_avg$y, xlab = 'Pseudotime', ylab = paste0(gene, ' Expression'), cex = 3, col = colors[color_num], lwd = 4.0) 
  color_num= color_num+1
}
legend('topleft', legend = c("ATOH1", "GFI1", "POU4F3"), col = colors, lty = 1, cex = 2.5, lwd = 4)
dev.off()

## Fig 5B
png('fig/Figure5/Figure5B_p15.png', width = 900, height = 850)
par(
  mar      = c(8, 8, 5, 5),
  cex.axis = 2.5,
  cex.lab  = 2.5,
  cex.main = 3.5
)
plot(NULL, NULL,xlim = c(0,1), ylim = c(0,1.5), ylab = "NOTCH1 Expression", xlab = "Pseudotime" )
colors = c("black")
color_num = 1
for(gene in c("NOTCH1")){
  kern_avg = ksmooth(psd_t, coch.p15[["RNA"]][gene ,as.vector(filt)], kernel = "normal", bandwidth = 10)
  lines(kern_avg$x/max(kern_avg$x),  kern_avg$y, xlab = 'Pseudotime', ylab = paste0(gene, ' Expression'), cex = 3, col = colors[color_num], lwd = 4.0) 
  color_num= color_num+1
}
dev.off()

## Fig S3
png('fig/FigureS3/Density_p15.png', width = 900, height = 850)
par(
  mar      = c(6, 6, 3, 3),
  cex.axis = 2.5,
  cex.lab  = 2.5,
  cex.main = 4
)
tt <- apply(psd_t,1,vec_pairwise_sub,kern_avg$x)
log_col = log(rowSums(tt<2.5))
plot(kern_avg$x/max(kern_avg$x), log_col, xlab = "Pseudotime", ylab = "Log Cell Density", type = "l", lwd = 2)

