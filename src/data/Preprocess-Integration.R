library(Seurat)
library(RColorBrewer)
library(SeuratDisk)
library(SeuratData)
library(Matrix)
library(class)
library(cccd)
library(igraph)
library(umap)

#Initializing
set.seed(4202023)
coch.p8.list = list()
coch.p15.list = list()

#Loading Data
for(dir in list.dirs("data/raw", recursive = FALSE)){
  print(dir)
  p8 = grepl("p8", dir)
  
  count_path = paste0(dir,'/matrix.mtx')
  barcodes_path = paste0(dir,'/barcodes.tsv')
  features_path = paste0(dir,'/features.tsv')
  
  count = Matrix::readMM(count_path)
  barcodes = read.table(barcodes_path)
  features = read.table(features_path)
  group = rep(substring(dir, 13), dim(count)[2])
  
  if(mean(grepl("_", features[,1])) == 1){
    count = count[grepl("mm10___", features[,1]),]
    features = features[grepl("mm10___", features[,1]),]
    features[,1] = substring(features[,1], 8)
    features[,2] = substring(features[,2], 8)
  }
  colnames(count) = barcodes[,1]
  rownames(count) = toupper(features[,2])
  
  rownames(count)[which(duplicated(features[,2]))] = features[which(duplicated(features[,2])),1]
  
  rownames(count)[grepl("^Mt", features[,2])] = features[grepl("^Mt", features[,2]),2]
  
  
  
  
  coch = CreateSeuratObject(count)
  group = rep(substring(dir, 13), dim(count)[2])
  coch[["groups"]] = group
  
  cell_filt = (coch$nFeature_RNA >=200)
  gene_filt  = rowSums(count != 0) >=5
  coch = coch[gene_filt, cell_filt]  
  
  
  mt = PercentageFeatureSet(coch, pattern = "^MT-")
  coch[["percent.mt"]] = mt
  print(VlnPlot(coch, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  
  c = GetAssayData(object = coch, slot = "counts")
  
  
  cell_filt = (coch$nFeature_RNA <=5000) & (coch$nFeature_RNA >=750) & (coch$percent.mt <= 30)
  coch = coch[, cell_filt]  
  
  
  coch = NormalizeData(coch)
  coch = FindVariableFeatures(coch, selection.method = "vst", nfeatures = 2000)
  
  coch = list(dir = coch)
  if(p8){
    coch.p8.list = append(coch.p8.list, coch)
    
  }
  else{
    coch.p15.list = append(coch.p15.list, coch)
  }
  print(coch)
}

#Integrating Data
features.p8 = SelectIntegrationFeatures(object.list = coch.p8.list)
features.p15 = SelectIntegrationFeatures(object.list = coch.p15.list)

coch.p8.anchors <- FindIntegrationAnchors(object.list = coch.p8.list, anchor.features = features.p8)
coch.p8 = IntegrateData(anchorset = coch.p8.anchors)

coch.p15.anchors <- FindIntegrationAnchors(object.list = coch.p15.list, anchor.features = features.p15)
coch.p15 = IntegrateData(anchorset = coch.p15.anchors)

#Clustering
set.seed(4202023)
coch.p8 = ScaleData(coch.p8, verbose = FALSE)
coch.p8 = RunPCA(coch.p8)
coch.p8 <- FindNeighbors(coch.p8, dims = 1:50, reduction = "pca", verbose = FALSE)
coch.p8 <- FindClusters(coch.p8, resolution = .4, cluster.name = "unintegrated_clusters", verbose = FALSE)
coch.p8 <- RunUMAP(coch.p8, dims = 1:50, verbose = FALSE)

set.seed(4202023)
coch.p15 = ScaleData(coch.p15, verbose = FALSE)
coch.p15 = RunPCA(coch.p15)
coch.p15 <- FindNeighbors(coch.p15, dims = 1:50, reduction = "pca", verbose = FALSE)
coch.p15 <- FindClusters(coch.p15, resolution = .25, cluster.name = "unintegrated_clusters", verbose = FALSE)
coch.p15 <- RunUMAP(coch.p15, dims = 1:50, verbose = FALSE)

#SAVE DATA HERE
SaveH5Seurat(coch.p8, "data/processed/cochlear_p8" ,overwrite = TRUE)
SaveH5Seurat(coch.p15, "data/processed/cochlear_p15" ,overwrite = TRUE)
