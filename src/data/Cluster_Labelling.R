markers.p8 = FindAllMarkers(object = coch.p8)
markers.p15 = FindAllMarkers(object = coch.p15)

#Label Clusters
char_genes = read.csv('../data/processed/KollaCellClusterGenes.csv')

#METHOD 1
num_clust = length(levels(coch.p8$seurat_clusters))-1
hits_mat = NULL
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
  if(is.null(hits_mat)){
    hits_mat = as.matrix(hits)
  }
  else{
    hits_mat = cbind(hits_mat, hits)
  }
}
rownames(hits_mat) = levels(as.factor(char_genes[,2]))
colnames(hits_mat) = 0:num_clust

hits_mat


#METHOD 2
num_clust = length(levels(coch.p8$seurat_clusters))-1
hits_mat = NULL
for(clust in 0:num_clust){
  hits = c()
  for(cell_type in levels(as.factor(char_genes[,2]))){
    clust_mark = markers.p8[markers.p8[,'cluster'] == clust,'gene'][1:150]
    ct_filt = char_genes[,2] == cell_type
    ct_mark = char_genes[ct_filt,1][1:150]
    hits = c(hits,length(intersect(toupper(clust_mark), toupper(ct_mark))))
  }
  if(is.null(hits_mat)){
    hits_mat = as.matrix(hits)
  }
  else{
    hits_mat = cbind(hits_mat, hits)
  }
}
rownames(hits_mat) = levels(as.factor(char_genes[,2]))
colnames(hits_mat) = 0:num_clust

hits_mat


#SAVE HITS