library(Seurat)
library(SeuratDisk)
library(icesTAF)

coch.p8 =  LoadH5Seurat("data/processed/cochlear_p8.h5Seurat")
coch.p15 = LoadH5Seurat("data/processed/cochlear_p15.h5Seurat")

mkdir('data/processed/Celltype_Annotation/')

markers.p8 = FindAllMarkers(object = coch.p8)
markers.p15 = FindAllMarkers(object = coch.p15)
write.csv(markers.p8, 'data/processed/Celltype_Annotation/markers_p8.csv')
write.csv(markers.p15, 'data/processed/Celltype_Annotation/markers_p15.csv')


#Label Clusters
char_genes = read.csv('data/processed/KollaCellClusterGenes.csv')

#METHOD 1
num_clust = length(levels(coch.p8$seurat_clusters))-1
hits_med.p8 = NULL
for(clust in 0:num_clust){
  hits = c()
  for(cell_type in levels(as.factor(char_genes[,2]))){
    clust_mark = markers.p8[markers.p8[,'cluster'] == clust,'gene']
    ct_filt = char_genes[,2] == cell_type
    ct_mark = char_genes[ct_filt,1]
    matches = match(toupper(ct_mark),toupper(clust_mark))
    hits = c(hits, median(matches[!is.na(matches)]))
    
    #hits = c(hits,length(intersect(toupper(clust_mark), toupper(ct_mark))))
  }
  if(is.null(hits_med.p8)){
    hits_med.p8 = as.matrix(hits)
  }
  else{
    hits_med.p8 = cbind(hits_med.p8, hits)
  }
}
rownames(hits_med.p8) = levels(as.factor(char_genes[,2]))
colnames(hits_med.p8) = 0:num_clust


num_clust = length(levels(coch.p15$seurat_clusters))-1
hits_med.p15 = NULL
for(clust in 0:num_clust){
  hits = c()
  for(cell_type in levels(as.factor(char_genes[,2]))){
    clust_mark = markers.p15[markers.p15[,'cluster'] == clust,'gene']
    ct_filt = char_genes[,2] == cell_type
    ct_mark = char_genes[ct_filt,1]
    matches = match(toupper(ct_mark),toupper(clust_mark))
    hits = c(hits, median(matches[!is.na(matches)]))
    
    #hits = c(hits,length(intersect(toupper(clust_mark), toupper(ct_mark))))
  }
  if(is.null(hits_med.p15)){
    hits_med.p15 = as.matrix(hits)
  }
  else{
    hits_med.p15 = cbind(hits_med.p15, hits)
  }
}
rownames(hits_med.p15) = levels(as.factor(char_genes[,2]))
colnames(hits_med.p15) = 0:num_clust


#METHOD 2
num_clust = length(levels(coch.p8$seurat_clusters))-1
hits_num.p8 = NULL
for(clust in 0:num_clust){
  hits = c()
  for(cell_type in levels(as.factor(char_genes[,2]))){
    clust_mark = markers.p8[markers.p8[,'cluster'] == clust,'gene'][1:150]
    ct_filt = char_genes[,2] == cell_type
    ct_mark = char_genes[ct_filt,1][1:150]
    hits = c(hits,length(intersect(toupper(clust_mark), toupper(ct_mark))))
  }
  if(is.null(hits_num.p8)){
    hits_num.p8 = as.matrix(hits)
  }
  else{
    hits_num.p8 = cbind(hits_num.p8, hits)
  }
}
rownames(hits_num.p8) = levels(as.factor(char_genes[,2]))
colnames(hits_num.p8) = 0:num_clust

num_clust = length(levels(coch.p15$seurat_clusters))-1
hits_num.p15 = NULL
for(clust in 0:num_clust){
  hits = c()
  for(cell_type in levels(as.factor(char_genes[,2]))){
    clust_mark = markers.p15[markers.p15[,'cluster'] == clust,'gene'][1:150]
    ct_filt = char_genes[,2] == cell_type
    ct_mark = char_genes[ct_filt,1][1:150]
    hits = c(hits,length(intersect(toupper(clust_mark), toupper(ct_mark))))
  }
  if(is.null(hits_num.p15)){
    hits_num.p15 = as.matrix(hits)
  }
  else{
    hits_num.p15 = cbind(hits_num.p15, hits)
  }
}
rownames(hits_num.p15) = levels(as.factor(char_genes[,2]))
colnames(hits_num.p15) = 0:num_clust


#SAVE HITS
write.csv(hits_med.p8, 'data/processed/Celltype_Annotation/hits_med_p8.csv')
write.csv(hits_med.p15, 'data/processed/Celltype_Annotation/hits_med_p15.csv')
write.csv(hits_num.p8, 'data/processed/Celltype_Annotation/hits_num_p8.csv')
write.csv(hits_num.p15, 'data/processed/Celltype_Annotation/hits_num_p15.csv')
