---
title: "Clustering and annotation"
author: "kjyi,jslee"
date: '2020 12 24'
output: R_script
---
  
# AS knee single cell analysis
  
```{r}
library(Seurat, lib.loc = "~kjyi/R/x86_64-redhat-linux-gnu-library/3.6/")
library(tidyverse, lib.loc = "~kjyi/R/x86_64-redhat-linux-gnu-library/3.6/")
# setwd("~kjyi/Projects/as")
system("tree -L 2 .", intern=T) %>% print()
```

TCR-IMGT PB, SF1, SF2 (library replicates)

```{r}
SF1_RNA_totalseq_counts <- Seurat::Read10X("cellranger/SF1/outs/filtered_feature_bc_matrix/")
SF2_RNA_totalseq_counts <- Seurat::Read10X("cellranger/SF2/outs/filtered_feature_bc_matrix/")
PB_RNA_totalseq_counts <- Seurat::Read10X("cellranger/PB/outs/filtered_feature_bc_matrix/")

colnames(SF1_RNA_totalseq_counts$`Gene Expression`)  <- paste0("SF1_",colnames(SF1_RNA_totalseq_counts$`Gene Expression`))
colnames(SF2_RNA_totalseq_counts$`Gene Expression`)  <- paste0("SF2_",colnames(SF2_RNA_totalseq_counts$`Gene Expression`))
colnames(PB_RNA_totalseq_counts$`Gene Expression`)   <- paste0("PB_", colnames( PB_RNA_totalseq_counts$`Gene Expression`))
colnames(SF1_RNA_totalseq_counts$`Antibody Capture`) <- paste0("SF1_",colnames(SF1_RNA_totalseq_counts$`Antibody Capture`))
colnames(SF2_RNA_totalseq_counts$`Antibody Capture`) <- paste0("SF2_",colnames(SF2_RNA_totalseq_counts$`Antibody Capture`))
colnames(PB_RNA_totalseq_counts$`Antibody Capture`)  <- paste0("PB_", colnames( PB_RNA_totalseq_counts$`Antibody Capture`))

AS <- CreateSeuratObject(cbind(
  SF1_RNA_totalseq_counts$`Gene Expression`,
  SF2_RNA_totalseq_counts$`Gene Expression`,
  PB_RNA_totalseq_counts$`Gene Expression`))

featurebarcodes = rownames(SF1_RNA_totalseq_counts$`Antibody Capture`)
featurebarcodes = featurebarcodes[!grepl("Hash",featurebarcodes)]
hashbarcodes = rownames(SF1_RNA_totalseq_counts$`Antibody Capture`)
hashbarcodes = hashbarcodes[grepl("Hash",hashbarcodes)]

AS[["FB"]] <- CreateAssayObject(counts = cbind(
  SF1_RNA_totalseq_counts$`Antibody Capture`,
  SF2_RNA_totalseq_counts$`Antibody Capture`,
  PB_RNA_totalseq_counts$`Antibody Capture`)[featurebarcodes,])

AS[["HTO"]] <- CreateAssayObject(counts = cbind(
  SF1_RNA_totalseq_counts$`Antibody Capture`,
  SF2_RNA_totalseq_counts$`Antibody Capture`,
  PB_RNA_totalseq_counts$`Antibody Capture`)[hashbarcodes,])

```

```{r}
# names of feature barcodes
AS[["FB"]]
```

```{r}
AS[["HTO"]]
```

## Demultiplexing with Hash tag oligo antibody info (HTO)

<https://satijalab.org/seurat/v3.2/hashing_vignette.html> 이것 참고해서 함. 근데 normalize하면 이상하게 됨..

```{r}
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
AS <- NormalizeData(AS, assay = "HTO", normalization.method = "CLR")
AS <- HTODemux(AS, assay = "HTO", positive.quantile = 0.99)
```

```{r}
# Global classification results
table(AS$HTO_classification.global)
```

```{r}
FeatureScatter(AS, feature1 = "hto_Hash1-TotalSeqC", feature2 = "hto_Hash2-TotalSeqC")
FeatureScatter(AS, feature1 = "hto_Hash1-TotalSeqC", feature2 = "hto_Hash3-TotalSeqC")
FeatureScatter(AS, feature1 = "hto_Hash1-TotalSeqC", feature2 = "hto_Hash4-TotalSeqC")
FeatureScatter(AS, feature1 = "hto_Hash1-TotalSeqC", feature2 = "hto_Hash5-TotalSeqC")
FeatureScatter(AS, feature1 = "hto_Hash2-TotalSeqC", feature2 = "hto_Hash3-TotalSeqC")
FeatureScatter(AS, feature1 = "hto_Hash2-TotalSeqC", feature2 = "hto_Hash4-TotalSeqC")
FeatureScatter(AS, feature1 = "hto_Hash2-TotalSeqC", feature2 = "hto_Hash5-TotalSeqC")
FeatureScatter(AS, feature1 = "hto_Hash3-TotalSeqC", feature2 = "hto_Hash4-TotalSeqC")
FeatureScatter(AS, feature1 = "hto_Hash3-TotalSeqC", feature2 = "hto_Hash5-TotalSeqC")
FeatureScatter(AS, feature1 = "hto_Hash4-TotalSeqC", feature2 = "hto_Hash5-TotalSeqC")
```

```{r}

set.seed(42)

res.km_1 <- kmeans(x = t(AS@assays$HTO@data[c("Hash1-TotalSeqC","Hash3-TotalSeqC"),]),
                   centers = 4,)
plot(t(AS@assays$HTO@data[c("Hash1-TotalSeqC","Hash3-TotalSeqC"),]),
     col=res.km_1$cluster)


```

kmeans가 잘 안돼서 매뉴얼하게 시도 (보이는 대로 끊음)

```{r}

mycol1 = rep("missing",ncol(AS@assays$HTO@data))

mycol1[AS@assays$HTO@data[c("Hash1-TotalSeqC"),] > 1 & 
         AS@assays$HTO@data[c("Hash3-TotalSeqC"),] < 1.5 ] = "Hash1"

mycol1[AS@assays$HTO@data[c("Hash1-TotalSeqC"),] < 1.5 & 
         AS@assays$HTO@data[c("Hash3-TotalSeqC"),] > 1.5 ] = "Hash3"

mycol1[AS@assays$HTO@data[c("Hash1-TotalSeqC"),] > 1.5 & 
         AS@assays$HTO@data[c("Hash3-TotalSeqC"),] > 1.5 ] = "Doublet_1_3"
table(mycol1)
```

```{r}
plot(t(AS@assays$HTO@data[c("Hash1-TotalSeqC","Hash3-TotalSeqC"),]),
     col=factor(mycol1))
```

```{r}
plot(t(AS@assays$HTO@data[c("Hash2-TotalSeqC","Hash4-TotalSeqC"),]),
     col=factor(mycol1)); abline(v = 3, h = 3)
```

```{r}
mycol1[AS@assays$HTO@data[c("Hash2-TotalSeqC"),] > 3 & 
         AS@assays$HTO@data[c("Hash4-TotalSeqC"),] < 3 ] = "Hash2"

mycol1[AS@assays$HTO@data[c("Hash2-TotalSeqC"),] < 3 & 
         AS@assays$HTO@data[c("Hash4-TotalSeqC"),] > 3 ] = "Hash4"

mycol1[AS@assays$HTO@data[c("Hash2-TotalSeqC"),] > 3 & 
         AS@assays$HTO@data[c("Hash4-TotalSeqC"),] > 3 ] = "Doublet_2_4"
table(mycol1)
```

```{r}
par(mfrow=c(1,2))
plot(t(AS@assays$HTO@data[c("Hash1-TotalSeqC","Hash3-TotalSeqC"),]),
     col=factor(mycol1))
plot(t(AS@assays$HTO@data[c("Hash2-TotalSeqC","Hash4-TotalSeqC"),]),
     col=factor(mycol1)); abline(v = 3, h = 3)
```

```{r}
AS$batch <- colnames(AS) %>% str_replace("_.*","")
table(AS$batch, mycol1)
```

```{r}
AS$hash.manual <- mycol1

AS <- AS[,AS$hash.manual %in% c("Hash1","Hash2","Hash3","Hash4")]
write_rds(AS, "data/AS.hash_double_filtered.Rds", compress = "gz")
```

```{r}
length(mycol1)
AS
```

# hash tag doublet filter : 35938 → 33573

```{r}
#저장
write_rds(AS, "data/AS.hash_double_filtered.Rds", compress = "gz")
```

# Basic QC: mitochondria, RNA, UMI count

```{r fig.height=4, fig.width=10}
AS[["percent.mt"]] <- PercentageFeatureSet(AS, pattern = "^MT-")
VlnPlot(AS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,pt.size = 0, group.by = "hash.manual")
VlnPlot(AS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,pt.size = 0, group.by = "batch")
VlnPlot(AS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,pt.size = 0)
```

```{r fig.height = 12, fig.width = 5}
VlnPlot(AS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 1,pt.size = 0)
```

```{r}
plot(density(AS$percent.mt)); abline(v = 10, h = 0 )
```

```{r}
tmp = na.omit(cbind(AS$percent.mt, AS$nCount_RNA))
tmp = tmp[tmp[,1] != Inf & tmp[,2] != Inf,]
smoothScatter(tmp);rm(tmp)
```

```{r}
table(AS$percent.mt < 10)
```

```{r}
AS = AS[,AS$percent.mt < 10]
# save
write_rds(AS, "data/AS.hash_double_filtered.mito_filtered.Rds", compress = "gz")
```

## Basic Seurat processing: clustering and UMAP

```{r}
AS <- AS %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData %>%
  RunPCA %>%
  FindNeighbors %>%
  FindClusters %>%
  RunUMAP(dims= 1:20)
write_rds(AS, "data/AS.hash_double_filtered.mito_filtered.Rds", compress = "gz")
```

## explore and annotation of cluster

```{r fig.height = 6, fig.width = 7}
DimPlot(AS, label=T) + NoLegend()
```

```{r fig.height = 6, fig.width = 7}
DimPlot(AS, group.by = "batch", label=T) + NoLegend()
```

```{r fig.height = 6, fig.width = 7}
DimPlot(AS, group.by = "hash.manual", label=T) + NoLegend()
```

```{r fig.height=14, fig.width=10}

FeaturePlot(AS,features = c("hto_Hash1-TotalSeqC",
                            "hto_Hash2-TotalSeqC",
                            "hto_Hash3-TotalSeqC",
                            "hto_Hash4-TotalSeqC",
                            "hto_Hash5-TotalSeqC"))

```

```{r, fig.height=20, fig.width=10}
goi = c("CD3D","CD4","KRT7",
        "PTPRC","CD1A","DNTT",
        "CD4","CD3E","GZMB",
        "CD68","IFNG", "CD79A")
FeaturePlot(AS,features = goi,ncol=2)
```

# Normalize antibody feature plot (raw count / log-normalized count)

```{r}
AS <- NormalizeData(AS, assay = "FB")
write_rds(AS, "data/AS.hash_double_filtered.mito_filtered.Rds", compress = "gz")
```

```{r, fig.height=30, fig.width=10}
all_totalseq_markers <- AS[["FB"]]@counts@Dimnames[[1]] %>% paste0("fb_",.)
FeaturePlot(AS,features = all_totalseq_markers, ncol=3)
```

```{r, fig.height=30, fig.width=10}
goi <- c("CD3E","CD4","CD8A","CD14","NCAM1", "IFNG", "TNF", "HLA-DQA2", "IL1B", "IL6", "FCGR3A")
FeaturePlot(AS,features = goi, ncol=2,keep.scale = "feature") # "all"
```