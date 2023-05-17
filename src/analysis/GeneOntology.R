library(Seurat)
library(SeuratDisk)
library(icesTAF)

p8 =  LoadH5Seurat("data/processed/cochlear_p8_Annotated.h5Seurat")
p15 = LoadH5Seurat("data/processed/cochlear_p15_Annotated.h5Seurat")

mkdir('data/processed/Ontology_input/')

p8 <- RenameIdents(object = p8, 
                   "0" = "Deiters",
                   "1" = "Deiters",
                   "2" = "Deiters",
                   "3" = "Deiters",
                   "4" = "Interdental",
                   "5" = "LGER",
                   "6" = "Hair Cell",
                   "7" = "Undefined",
                   "8" = "Oc90+",
                   "9" = "Interdental",
                   "10" = "Oc90+",
                   "11" = "Undefined",
                   "12" = "OS",
                   "13" = "Undefined",
                   "14" = "Hair Cell",
                   "15" = "LGER",
                   "16" = "LGER",
                   "17" = "Glia",
                   "18" = "OS",
                   "19" = "Glia",
                   "20" = "Undefined",
                   "21" = "Undefined",
                   "22" = "Undefined",
                   "23" = "Undefined")
DE_p8_Deiters <- FindMarkers(p8, ident.1 = "rosa26WT_p8", group.by = 'groups', subset.ident = "Deiters")
write.csv(DE_p8_Deiters, 'data/processed/Ontology_input/DE_p8_Deiters.csv', row.names=TRUE)
DE_p8_Glia <- FindMarkers(p8, ident.1 = "rosa26WT_p8", group.by = 'groups', subset.ident = "Glia")
write.csv(DE_p8_Glia, 'data/processed/Ontology_input/DE_p8_Glia.csv', row.names=TRUE)
DE_p8_OC90 <- FindMarkers(p8, ident.1 = "rosa26WT_p8", group.by = 'groups', subset.ident = "Oc90+")
write.csv(DE_p8_OC90, 'data/processed/Ontology_input/DE_p8_roof.csv', row.names=TRUE)
DE_p8_Interdental <- FindMarkers(p8, ident.1 = "rosa26WT_p8", group.by = 'groups', subset.ident = "Interdental")
write.csv(DE_p8_Interdental, 'data/processed/Ontology_input/DE_p8_Interdental.csv', row.names=TRUE)
DE_p8_LGER <- FindMarkers(p8, ident.1 = "rosa26WT_p8", group.by = 'groups', subset.ident = "LGER")
write.csv(DE_p8_LGER, 'data/processed/Ontology_input/DE_p8_LGER.csv', row.names=TRUE)


p15 <- RenameIdents(object = p15, 
                    "0" = "Interdental",
                    "1" = "Undefined",
                    "2" = "Deiters",
                    "3" = "Undefined",
                    "4" = "Deiters",
                    "5" = "LGER",
                    "6" = "Glia",
                    "7" = "Deiters",
                    "8" = "LGER",
                    "9" = "Undefined",
                    "10" = "Undefined",
                    "11" = "Hair cells",
                    "12" = "Undefined",
                    "13" = "Undefined",
                    "14" = "Oc90+",
                    "15" = "Undefined",
                    "16" = "Undefined",
                    "17" = "LGER",
                    "18" = "Undefined",
                    "19" = "Undefined",
                    "20" = "Undefined",
                    "21" = "Undefined")
DE_p15_Deiters <- FindMarkers(p15, ident.1 = "rosa26WT_p15", group.by = 'groups', subset.ident = "Deiters")
write.csv(DE_p15_Deiters, 'data/processed/Ontology_input/DE_p15_Deiters.csv', row.names=TRUE)
DE_p15_Glia <- FindMarkers(p15, ident.1 = "rosa26WT_p15", group.by = 'groups', subset.ident = "Glia")
write.csv(DE_p15_Glia, 'data/processed/Ontology_input/DE_p15_Glia.csv', row.names=TRUE)
DE_p15_OC90 <- FindMarkers(p15, ident.1 = "rosa26WT_p15", group.by = 'groups', subset.ident = "Oc90+")
write.csv(DE_p15_OC90, 'data/processed/Ontology_input/DE_p15_roof.csv', row.names=TRUE)
DE_p15_Interdental <- FindMarkers(p15, ident.1 = "rosa26WT_p15", group.by = 'groups', subset.ident = "Interdental")
write.csv(DE_p15_Interdental, 'data/processed/Ontology_input/DE_p15_Interdental.csv', row.names=TRUE)
DE_p15_LGER <- FindMarkers(p15, ident.1 = "rosa26WT_p15", group.by = 'groups', subset.ident = "LGER")
write.csv(DE_p15_LGER, 'data/processed/Ontology_input/DE_p15_LGER.csv', row.names=TRUE)