---
title: "CD4 T in 6xAS"
author: "SWG"
date: '2022 10 23'
output: R_script
---

Idents(SO.CD4.harmony) <- 'integrated_snn_res.0.8'
DimPlot(SO.CD4.harmony, label = T, reduction = 'umap')
DefaultAssay(SO.CD4.harmony) <- 'RNA'
Idents(SO.CD4.harmony) <- "integrated_snn_res.0.6"
SO.CD4 <- subset(SO.T_NK.harmony, idents = c('0', '1', '4', '12', '14'))
SO.CD4 <- subset(SO.CD4.harmony, idents = c('8'), invert = T)

DefaultAssay(SO.CD4) <- 'RNA'
SO.CD4.list <- SplitObject(object = SO.CD4, split.by = "ident.study")
SO.CD4.list <- lapply(X = SO.CD4.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features.test <- SelectIntegrationFeatures(object.list = SO.CD4.list)
test.anchors <- FindIntegrationAnchors(object.list = SO.CD4.list, anchor.features = features.test) #defaultassay 사용됨
SO.CD4.combined <- IntegrateData(anchorset = test.anchors)

DefaultAssay(SO.CD4.combined) <- 'integrated'
SO.CD4.combined <- SO.CD4.combined %>%
  ScaleData %>%
  RunPCA  %>%
  FindNeighbors %>%
  FindClusters %>%
  RunUMAP(dims= 1:7)

DimPlot(SO.CD4.combined, split.by = 'ident.study', reduction = 'umap') + NoLegend()

library(harmony)
SO.CD4.harmony <- SO.CD4.combined %>%
  RunHarmony("ident.name", assay.use="integrated")

# UMAP and clustering with harmonized PCs
SO.CD4.harmony <- RunUMAP(SO.CD4.harmony, reduction='harmony', dims = 1:10) # 1:30 → 1:35
SO.CD4.harmony <- FindNeighbors(SO.CD4.harmony, reduction='harmony')
SO.CD4.harmony <- FindClusters(SO.CD4.harmony, resolution = 0.6) # 0.3 → 0.1
#save
saveRDS(seurat_obj, file=paste0(out_data_dir, 'seurat_object_sct_harmony.rds'))
p5 <- DimPlot(SO.CD4.combined, group.by = 'ident.name', reduction = 'umap') + ggtitle('before harmony')
p6 <- DimPlot(SO.CD4.harmony, group.by = 'ident.name', reduction = 'umap') + ggtitle('after harmony')
CombinePlots(list(p5, p6))

SO.CD4.harmony <- FindNeighbors(SO.CD4.harmony, reduction='harmony', dims = 1:10)
SO.CD4.harmony <- FindClusters(SO.CD4.harmony, resolution = 0.6) # 선택
SO.CD4.harmony <- RunUMAP(SO.CD4.harmony, reduction='harmony', dims = 1:10) # 1:30 → 1:35


DimPlot(SO.CD4.combined, group.by = 'ident.name', split.by = 'ident.name',reduction = 'umap', ncol = 3) + ggtitle('after study')
DimPlot(SO.CD4.harmony, group.by = 'ident.name', split.by = 'ident.name',reduction = 'umap', ncol = 3) + ggtitle('after harmony')
DimPlot(SO.CD4.combined, group.by = 'origin', split.by = 'origin', reduction = 'umap')
DimPlot(SO.CD4.combined, group.by = 'ident.sample', split.by = 'ident.sample', reduction = 'umap', ncol = 4)
DimPlot(SO.CD4.harmony, group.by = 'origin', split.by = 'origin', reduction = 'umap')
DimPlot(SO.CD4.harmony, group.by = 'ident.sample', split.by = 'ident.sample', reduction = 'umap', ncol = 4) + NoLegend()
DimPlot(SO.CD4.harmony, group.by = 'ident.name', split.by = 'ident.name',reduction = 'umap', ncol = 3) + NoLegend()


DefaultAssay(SO.CD4.harmony) <- "FB"
Idents(SO.spare1) <- 'ident.study'
test1 <- subset(SO.spare1, idents = 'old')
DefaultAssay(test1) <- "FB"
p1 <- FeaturePlot(test1, "CD45RO-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD45RO protein")
p2 <- FeaturePlot(test1, "CD27-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD27 protein")
p3 <-FeaturePlot(test1, "CD134-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD134 protein")
p4 <- FeaturePlot(test1, "CD270-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD270(TNFRSF14) protein")
p5 <- FeaturePlot(test1, "CD38-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD38 protein")
p6 <- FeaturePlot(test1, "CD39-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD39 protein")
CombinePlots(list(p1, p2, p3, p4, p5, p6), ncol = 3)
DefaultAssay(SO.CD4.harmony) <- 'integrated'
p3
p2|p3|p4

DefaultAssay(SO.CD4.harmony) <- 'RNA'
OxGi <- c("TNFRSF4", "TNFRSF18", "CD252", "GITRL", "GPR65")
FeaturePlot(SO.CD4.harmony, features = "IL2RA")
FeaturePlot(SO.CD4.harmony, features = c("CD8A", "CD4"))

DotPlot(SO.T_NK.harmony, features = goi4, cols = c("lightgrey", "red"), dot.scale = 8)+
  RotatedAxis()
DimPlot(SO.T_NK.harmony, reduction = 'umap', label = T)
DefaultAssay(SO.CD4.harmony) <- 'RNA'
DotPlot(test.CD4.harmony, features = goi4, cols = c("lightgrey", "red"), dot.scale = 8)+
  RotatedAxis()
fig3 <- c("IL17A", "IFNG", "CCR6", "CXCR3", "TBX21", "KLRB1", "TNFRSF4", "TNFRSF18", "GPR65")
DotPlot(SO.CD4.harmony, features = fig3,cols = c("lightgrey", "red"), dot.scale = 8)+
  RotatedAxis()
DimPlot(SO.CD4.harmony, label = T)
table(Idents(SO.CD4.harmony))
GetAssay(SO.CD4.harmony)

DotPlot(SO.CD4.harmony, features = "FOXP3", cols = c("lightgrey", "red"), dot.scale = 8)+
  RotatedAxis()
view(p1$data)

FeaturePlot(test.CD4.harmony, features = fig3)
goi4.1 <- c("GPR65", "TNFRSF4", "TNFRSF18", "CCR6", "FOXP3", "TBX21", "CXCR3", "CCR7", "IFNG")

# split further
DefaultAssay(SO.CD4.harmony) <- 'integrated'
test <- FindClusters(SO.CD4.harmony, resolution = c((7:12)/10))
p1 <- DimPlot(test, reduction='umap', group.by = "integrated_snn_res.0.7", label=T) + NoLegend()
p2 <- DimPlot(test, reduction='umap', group.by = "integrated_snn_res.0.8", label=T)+ NoLegend() 
p3 <- DimPlot(test, reduction='umap', group.by = "integrated_snn_res.0.9", label=T)+ NoLegend()
p4 <- DimPlot(test, reduction='umap', group.by = "integrated_snn_res.1", label=T)+ NoLegend()
p5 <- DimPlot(test, reduction='umap', group.by = "integrated_snn_res.1.1", label=T)+ NoLegend()
p6 <- DimPlot(test, reduction='umap', group.by = "integrated_snn_res.1.2", label=T)+ NoLegend()

CombinePlots(plots = list(p1, p2, p3, p4, p5, p6), ncol = 3)

goi4

```{r fig.height = 6, fig.width = 7}
DimPlot(SO.CD4.harmony, group.by = "hash.manual", label=T) + NoLegend()
```
DefaultAssay(SO.T_NK.harmony) <- 'RNA'
Idents(SO.CD4.harmony) <- "integrated_snn_res.0.8"
FeaturePlot(SO.T_NK.harmony, features = c("FAS", "SELL", "IL2RB", "CD27"))
FeaturePlot(test.CD4.harmony, features = goi3, ncol = 3)
FeaturePlot(SO.CD4.harmony, features = OxGi, ncol = 3)


# Cluster proportion
pie3D(x=table(Idents(SO.CD4.harmony)), 
      labels=CD4.piepct, init.angle = 90, clockwise = T, explode = 0.1,
      main = "Cluster proportion", col = rainbow(length(table(Idents(SO.CD4.harmony)))))
CD4.piepct<- round(100 * table(Idents(SO.CD4.harmony)) / sum(table(Idents(SO.CD4.harmony))), 1)
legend("upright",SO.CD4.harmony@active.ident,
       cex = 0.5, fill = rainbow(length(table(Idents(SO.CD4.harmony)))))

mytable <- table(Idents(SO.CD4.harmony))
lbls <- paste(names(mytable), "\n", CD4.piepct, sep = "")
pie3D(mytable, labels = lbls, explode = 0.1,
      main = "Cluster proportion\n (pct)")


# Co-positive
SO.CD4.harmony$copositive_gpr65_ox40 <- "no"
copos = SO.CD4.harmony & SO.CD4.harmony@assays$RNA@counts["TNFRSF4",] > 0
SO.CD4.harmony$copositive_gpr65_ox40[copos] <- "yes"
DimPlot(SO.CD4.harmony, group.by = "copositive_gpr65_ox40")
SO.CD4.harmony$copositive_pseudo = "sample(copos)"
DimPlot(SO.CD4.harmony, group.by = "copositive_pseudo")

SO.CD4.harmony$cop_17_OX <- "x"
cop_17OX <- SO.CD4.harmony@assays$RNA@counts["RORC",] > 0  & SO.CD4.harmony@assays$RNA@counts["TNFRSF4",] > 0
SO.CD4.harmony$cop_17_OX[cop_17OX] <- "17_OX"
DimPlot(SO.CD4.harmony, group.by = "cop_17_OX")
DimPlot(SO.CD4.harmony, group.by = "cop_17OXGI")

cop_OXGI <- SO.CD4.harmony@assays$RNA@counts["TNFRSF4",] > 0  & SO.CD4.harmony@assays$RNA@counts["TNFRSF18",] > 0
SO.CD4.harmony$cop_OX_GI <- "etc."
SO.CD4.harmony$cop_OX_GI[cop_OXGI] <- "OX_GI"
DimPlot(SO.CD4.harmony, group.by = "cop_OX_GI", split.by = 'ident.sample', ncol =4)
DimPlot(SO.CD4.harmony, group.by = "cop_OX_GI", pt.size = 0.8)
DimPlot(SO.CD4.harmony, group.by = "cop_OX_GI", split.by = "origin",pt.size = 0.8)

cop_OXGPR65 <- SO.CD4.harmony@assays$RNA@counts["TNFRSF4",] > 0  & SO.CD4.harmony@assays$RNA@counts["GPR65",] > 0
SO.CD4.harmony$cop_OX_GPR65 <- "etc."
SO.CD4.harmony$cop_OX_GPR65[cop_OXGPR65] <- "OX_GPR65"

cop_GIGPR65 <- SO.CD4.harmony@assays$RNA@counts["TNFRSF18",] > 0  & SO.CD4.harmony@assays$RNA@counts["GPR65",] > 0
SO.CD4.harmony$cop_GI_GPR65 <- "etc."
SO.CD4.harmony$cop_GI_GPR65[cop_GIGPR65] <- "GI_GPR65"

cop_triple <-SO.CD4.harmony@assays$RNA@counts["TNFRSF18",] > 0  & SO.CD4.harmony@assays$RNA@counts["GPR65",] > 0 &
  SO.CD4.harmony@assays$RNA@counts["TNFRSF4",] > 0
SO.CD4.harmony$cop_OXGI_GPR65 <- "etc."
SO.CD4.harmony$cop_OXGI_GPR65[cop_triple] <- "cop_triple"

DimPlot(SO.CD4.harmony, group.by = "cop_OX_GPR65", pt.size = 0.8)
DimPlot(SO.CD4.harmony, group.by = "cop_GI_GPR65", pt.size = 0.8)
DimPlot(SO.CD4.harmony, group.by = "cop_OXGI_GPR65", pt.size = 0.8)

a <- table(SO.CD4.harmony$manual_cluster, SO.CD4.harmony$cop_OX_GI)
b <- table(SO.CD4.harmony$origin, SO.CD4.harmony$cop)

write.xlsx(a, 'oxgi.xlsx')

# CCR6, CXCR3
cop_C6CX3 <- SO.CD4.harmony@assays$RNA@counts["CCR6",] > 0  & SO.CD4.harmony@assays$RNA@counts["CXCR3",] > 0
SO.CD4.harmony$cop_C6_CX3 <- "etc."
SO.CD4.harmony$cop_C6_CX3[cop_C6CX3] <- "CCR6_CXCR3"
DimPlot(SO.CD4.harmony, group.by = "cop_C6_CX3", pt.size = 0.4)
a <- table(SO.CD4.harmony$manual_cluster, SO.CD4.harmony$cop_C6_CX3)
write.xlsx(a, 'CCR6_CXCR3.xlsx')

cop_Tp_KLRB1 <- SO.CD4.harmony@assays$RNA@counts["CCR6",] > 0  & SO.CD4.harmony@assays$RNA@counts["CXCR3",] > 0 &
  SO.CD4.harmony@assays$RNA@counts["KLRB1",]
SO.CD4.harmony$cop_KLRB1 <- "etc."
SO.CD4.harmony$cop_KLRB1[cop_Tp_KLRB1] <- "Tp_KLRB1"
DimPlot(SO.CD4.harmony, group.by = "cop_KLRB1", pt.size = 0.4)
b <- table(SO.CD4.harmony$manual_cluster, SO.CD4.harmony$cop_KLRB1)
write.xlsx(b, 'CCR6_Triple_KLRB1.xlsx')

table(SO.CD4.harmony$origin, SO.CD4.harmony$cop_KLRB1)

# IL-17A too few
b <- table(SO.CD4.harmony$cop_17_OX[cop_17OX])
table(Idents(SO.CD4.harmony))
a <- table(SO.CD4.harmony@assays$RNA@counts["RORC",] > 0)
d <- table(SO.CD4.harmony$cop_17OXGI[cop_17OXGI])
c <- table(SO.CD4.harmony$cop_17GI[cop_17GI])

SO.CD4.harmony$cop_17GI <- "x"
cop_17GI <- SO.CD4.harmony@assays$RNA@counts["RORC",] > 0 &
  SO.CD4.harmony@assays$RNA@counts["TNFRSF18",] > 0
SO.CD4.harmony$cop_17GI[cop_17GI] <- "cop_17_GI"

df = data.frame(x = c( "RORC_OX", "RORC_GI","RORC_OXGI"), y= c(b, c, d))
df
view(df)

table(SO.CD4.harmony@assays$RNA@counts["TNFRSF4",] >0 & SO.CD4.harmony@assays$RNA@counts["TNFRSF18",] >0)
table(SO.CD4.harmony@assays$RNA@counts["TNFRSF4",] >0 & SO.CD4.harmony@assays$RNA@counts["TNFRSF18",] >0 &
        SO.CD4.harmony@assays$RNA@counts["FOXP3",] == 0 & SO.CD4.harmony@assays$RNA@counts["RORC",] >0)
SO.CD4.harmony$cop_all <- "etc."
cop_4all <- SO.CD4.harmony@assays$RNA@counts["TNFRSF4",] >0 & SO.CD4.harmony@assays$RNA@counts["TNFRSF18",] >0 &
  SO.CD4.harmony@assays$RNA@counts["FOXP3",] == 0 & SO.CD4.harmony@assays$RNA@counts["RORC",] >0
SO.CD4.harmony$cop_all[cop_4all] <- "target"
DimPlot(SO.CD4.harmony, group.by = "cop_all")

table(Idents(SO.CD4.harmony), SO.CD4.harmony$cop_all)
table(SO.CD4.harmony$origin, SO.CD4.harmony@active.ident)
table(SO.CD4.harmony$origin, SO.CD4.harmony$seurat_clusters)
table(SO.CD4.harmony$origin)
a <- table(SO.CD4.harmony$ident.sample, SO.CD4.harmony$seurat_clusters)
write.xlsx(a, 'cluster proportion')

# vlnplot
library(patchwork)
plots <- VlnPlot(SO.CD4.harmony, features = c("TNFRSF4", "TNFRSF18", "GPR65"), split.by = "origin", group.by = "manual_cluster",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# annotation
FeaturePlot(SO.CD4.harmony, features = c("PRF1", "CCL4", "CCL5", "GZMA", "GZMH", "GNLY", "NKG7", "CST7"))
FeaturePlot(SO.CD4.harmony, features = "PTPRC", min.cutoff = 'q50')

# export by exel
install.packages("xlsx")
library("xlsx")
library(openxlsx)
df <- table(SO.CD4.harmony@active.ident, SO.CD4.harmony$ident.name)
df
write.xlsx(df, 'pie2')


### Topgene 
# 1) total
test$manual_cluster <- c(
  "0"="T cell",
  "1"="SF Monocyte",
  "2"="PB Monocyte",
  "3"="NK",
  "4"="DC",
  "5"="B cell",
  "6"="Uncategorized",
  "7"="pDC",
  "8"="activated T cell",
  "9"="Mk"
)[test$seurat_clusters]

Idents(SO.sum.harmony) <- "integrated_snn_res.0.1"
all_markers_seuratcluster <- FindAllMarkers(SO.sum.harmony)
all_markers_seuratcluster %>% write_csv("S.allmarkers_seuratcluster.csv")

Idents(SO.sum.harmony) <- "manual_cluster"
S.all_markers_manualcluster <- FindAllMarkers(SO.sum.harmony)
S.all_markers_manualcluster %>% write_csv("S.allmarkers_manualcluster.csv")

str(S.all_markers_manualcluster)
S.all_markers_manualcluster %>% group_by(cluster) %>% filter(, between(row_number(), 1,20)) %>% arrange(as.character(cluster)) %>%
  {colnames(.)=c("p-value", "average log2 fold-change", "% of expressing cells in the cluster", "% of expressing cells not in the cluster", "Adjusted p-value", "Cluster", "Gene Symbol");.} %>%
  write_csv("tot_marker_top_20genes_per_manual_cluster.csv")
```

# 2) T, NK
DimPlot(SO.T_NK.harmony, label = T)
SO.T_NK.harmony <- test1
SO.T_NK.harmony$manual_cluster <- c(
  "0"="CD4, naive",
  "1"="CD4, memory",
  "2"="CD4, memory",
  "3"="CD8, memory",
  "4"="NK",
  "5"="CD4, naive",
  "6"="CD8, memory",
  "7"="CD4, memory",
  "8"="CD4, naive",
  "9"="CD4, memory",
  "10"="CD8, naive",
  "11"="NK", 
  "12"="CD4, Treg",
  "13"="CD8, memory",
  "14"="CD4, memory",
  "15"="CD8, memory",
  "16"="NK",
  "17"="Gamma-Delta/MAIT",
  "18"="NK",
  "19"="CD8, memory",
  "20"="CD4, memory",
  "21"="Rapidly proliferating",
  "22"="NK",
  "23"="CD4, memory",
  "24"="Rapidly proliferating",
  "25"="Uncategorized",
  "26"="CD8, memory"
)[SO.T_NK.harmony$seurat_clusters]

Idents(SO.T_NK.harmony) <- 'manual_cluster'
DimPlot(SO.T_NK.harmony, reduction = 'umap')

Idents(SO.T_NK.harmony) <- "integrated_snn_res.1.5"
T.NK_markers_seuratcluster <- findallmarkers_mc(SO.T_NK.harmony)
a%>% write_csv("T.NK_markers_seuratcluster.csv")
# T.NK_markers_seuratcluster%>% write.xlsx("T.NK_markers_seuratcluster.csv")

library(data.table)
a <- rbindlist(T.NK_markers_seuratcluster)
head(rbindlist(T.NK_markers_seuratcluster))
a%>% write_csv("T.NK_markers_seuratcluster.csv")


Idents(SO.T_NK.harmony) <- "manual_cluster"
T.NK_markers_manualcluster <- FindAllMarkers(SO.T_NK.harmony)
T.NK_markers_manualcluster %>% write_csv("T.NK_markers_manualcluster.csv")

DotPlot(SO.T_NK.harmony, features = goi.1) +
  RotatedAxis()
goi.1 <- c("CD3D", "CD3E", "CD4", "CD8A", "CCR7", "FOXP3", "IL2RA",
           "KLRB1", "TRBV6-4", "FCGR3A", "PPBP", "HBB")

T.NK_markers_manualcluster %>% group_by(cluster) %>% filter(, between(row_number(), 1,20)) %>% arrange(as.character(cluster)) %>% {colnames(.)=c("p-value", "average log2 fold-change", "% of expressing cells in the cluster", "% of expressing cells not in the cluster", "Adjusted p-value", "Cluster", "Gene Symbol");.} %>% write_csv("T.NK_marker_top_20genes_per_manual_cluster.csv")
T.NK_manual_list <- rbindlist(T.NK_markers_manualcluster)
head(T.NK_manual_list)
T.NK_manual_list%>% write_csv("T.NK_markers_manualcluster.csv")
Idents(SO.T_NK.harmony) <- 'manual_cluster'
DotPlot(SO.T_NK.harmony, features = goi) + RotatedAxis()

# 3) CD4 T 
SO.CD4.harmony <- test1.3.harmony
Idents(SO.CD4.harmony) <- 'seurat_clusters'
SO.CD4.harmony.all <- SO.CD4.harmony
SO.CD4.harmony <- subset(SO.CD4.harmony, idents = 6, invert = T)
DimPlot(SO.CD4.harmony, reduction= 'umap', label=T)
SO.CD4.harmony$manual_cluster <- c(
  "0"="Naive",
  "1"="Th1",
  "2"="Naive",
  "3"="Naive",
  "4"="Th17",
  "5"="CM",
  "6"="E/EM",
  "7"="CM",
  "8"="Naive",
  "9"="E/ME",
  "10"="EMRA(SF)",
  "11"="Th1", 
  "12"="pTh17",
  "13"="Naive",
  "14"="Uncategorized"
)[SO.CD4.harmony$seurat_clusters]

Idents(SO.CD4.harmony) <- 'manual_cluster'
DimPlot(SO.CD4.harmony, reduction = 'umap', label = T) + NoLegend()
SO.CD4.harmony <- FindClusters(SO.CD4.harmony, resolution = 1.1)
Idents(SO.CD4.harmony) <- "integrated_snn_res.1.1"
CD4_markers_seuratcluster.del <- FindAllMarkers(SO.CD4.harmony, logfc.threshold = 0.2)
?FindAllMarkers
CD4_markers_seuratcluster%>% write_csv("CD4_markers_seuratcluster.csv")
# T.NK_markers_seuratcluster%>% write.xlsx("T.NK_markers_seuratcluster.csv")

CD4_seurat_list <- dplyr::bind_rows(CD4_markers_seuratcluster)
head(CD4_seurat_list)
CD4_seurat_list%>% write_csv("CD4_markers_seuratcluster.csv")

a <- table(SO.CD4.harmony$ident.name, SO.CD4.harmony$seurat_clusters)
write.xlsx(a, 'CD4_name.xlsx')
table(SO.CD4.harmony$origin, SO.CD4.harmony$seurat_clusters)
table(SO.CD4.harmony$origin)

Idents(SO.CD4.harmony) <- "manual_cluster"
levels(as.factor(SO.CD4.harmony$manual_cluster))
CD4_markers_manualcluster <- FindAllMarkers(SO.CD4.harmony, logfc.threshold = 0.2)
CD4_markers_manualcluster %>% write_csv("CD4_markers_manualcluster.csv")


CD4_seurat_list <- dplyr::bind_rows(CD4_markers_seuratcluster)
head(CD4_seurat_list)
CD4_seurat_list%>% write_csv("CD4_markers_seuratcluster.csv")


CD4_markers_manualcluster %>% group_by(cluster) %>% filter(, between(row_number(), 1,20)) %>% arrange(as.character(cluster)) %>% {colnames(.)=c("p-value", "average log2 fold-change", "% of expressing cells in the cluster", "% of expressing cells not in the cluster", "Adjusted p-value", "Cluster", "Gene Symbol");.} %>% write_csv("CD4_marker_top_20genes_per_manual_cluster.csv")

CD4_manual_list <- rbindlist(CD4_markers_manualcluster)
head(CD4_manual_list)
CD4_manual_list%>% write_csv("CD4_markers_manualcluster.csv")

DefaultAssay(SO.CD4.harmony) <- 'RNA'
DotPlot(SO.CD4.harmony, features = c("PDE4D", "GPR65"), cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis()
SO.CD4.harmony <- FindClusters(test1, resolution = 1.2)
FeaturePlot(SO.CD4.harmony, features = c("PASK", "CCR7", "SELL", "KLF2"))
FeaturePlot(SO.CD4.harmony, features = goi4.1)

# proliferating?
FeaturePlot(SO.T.NK, features = c("TCF7", "CDKN1A", "MKI67", "CACYBP"))
FeaturePlot(SO.CD4.harmony, features = c("CDKN1A", "CDKN2D", "CDK2AP1", "CDK5RAP2")) # 10번은 proliferating
FeaturePlot(SO.T.NK, features = "HBB")

DimPlot(SO.T.NK, group.by = "origin", split.by = 'origin',reduction = 'umap')




