# AS 6 people(SF, PB pair) anaysis 4 + 2(lkj)

---
  title: "SF specific memory CD4 T in 6xAS"
author: "SWG"
date: '2022 08 31'
output: R_script
---
  
  ```

load("~/Assf/script/220922.RData")

# Pre-setting
```{r Loading packages, echo= FALSE}
setwd("/home/users/songwoogil/Assf")

library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(monocle3)
library(magrittr)
library(tidyr)
library(ggpubr)
library(tibble)
library(SeuratWrappers)
# library(SeuratWrappers, lib.loc = "/home/users/seongjinchoi/R/x86_64-redhat-linux-gnu-library/3.6/")
# remotes::install_github('satijalab/seurat-wrappers')
library(sctransform)
library(ggpubr)
library(pheatmap)
library(stringr)
library(viridisLite)
library(viridis) # required viridisLite

library(tidyverse)


## New pair
# Read raw count & Making Seurat object
```{r, eval = FALSE}
## Read raw count
setwd("/home/users/songwoogil/Assf")
new.sample <- list('counts/PBMC_F46_SFMC_F56_M57_GEX',
                   'counts/PBMC_M41_SFMC_F53_M52_GEX',
                   'counts/PBMC_M52_SFMC_F41_F46_GEX',
                   'counts/PBMC_M57_SFMC_F83_M41_GEX')

new.sample.demux <- list('demuxlet/PBMC_F46_SFMC_F56_M57_GEX.best',
                         'demuxlet/PBMC_M41_SFMC_F53_M52_GEX.best',
                         'demuxlet/PBMC_M52_SFMC_F41_F46_GEX.best',
                         'demuxlet/PBMC_M57_SFMC_F83_M41_GEX.best')

names(new.sample) <- c("PBMC_F46_SFMC_F56_M57", "PBMC_M41_SFMC_F53_M52", "PBMC_M52_SFMC_F41_F46", 
                       "PBMC_M57_SFMC_F83_M41")
new.sample

for(i in 1:length(new.sample)){
  tmp.matrix <- as.character(new.sample[[i]])
  tmp.matrix <- Read10X(tmp.matrix) ##
  
  tmp <- CreateSeuratObject(counts = tmp.matrix, project = "ASSF")
  
  tmp.demux <- as.character(new.sample.demux[[i]])
  tmp.demux <- read.csv(file = tmp.demux, sep = "\t", header = TRUE) ##
  tmp.demux <- tmp.demux %>% mutate(., singlet = sapply(str_split(tmp.demux$BEST, "-"), "[", 1)) 
  tmp.demux <- tmp.demux %>% filter(., BARCODE %in% colnames(tmp))
  # %>% filter(., singlet == "SNG")
  tmp <- tmp[, tmp.demux$BARCODE]
  
  tmp$orig.gem <- names(new.sample[i])
  tmp$singlet.demux <- tmp.demux$singlet
  tmp$ident.demux <- tmp.demux$SNG.1ST
  
  if(i == 1){
    SO.new <- tmp
  } else{
    SO.new <- merge(SO.new, y = tmp,
                    merge.data = TRUE,
                    merge.dr = NULL)
  }
  
  rm(tmp.matrix)
  rm(tmp.demux)
  rm(tmp)
}

Idents(SO.new) <- "ident.demux"


# test.RAin$ident.study <- ifelse(test = test.RAin$orig.ident %in% "ASSF", yes = "new", no = "old")
table(test.RAin$ident.study)
table(test.RAin$ident.sample)

SO.sum$ident.sample <- NA
SO.sum$ident.name <- NA
SO.sum$age <- NA
SO.sum$gender <- NA
SO.sum$ident.dz <- NA
SO.sum$origin <- NA
SO.sum$ident.study <- NA
SO.sum$CRP <- NA


# Combine data through ident.sample
SO.sum$ident.sample[SO.sum$ident.demux == "PBMC_AS517"] <- "PBMC_AS517"
SO.sum$ident.sample[SO.sum$ident.demux == "SFMC_SF433"] <- "SFMC_SF433"
SO.sum$ident.sample[SO.sum$ident.demux == "PBMC_AS426"] <- "PBMC_AS426"
SO.sum$ident.sample[SO.sum$ident.demux == "SFMC_SF308"] <- "SFMC_SF308"
SO.sum$ident.sample[SO.sum$ident.demux == "PBMC_AS210"] <- "PBMC_AS210"
SO.sum$ident.sample[SO.sum$ident.demux == "SFMC_SF355"] <- "SFMC_SF355"
SO.sum$ident.sample[SO.sum$ident.demux == "PBMC_AS509"] <- "PBMC_AS509"
SO.sum$ident.sample[SO.sum$ident.demux == "SFMC_SF409"] <- "SFMC_SF409"
SO.sum$ident.sample[SO.sum$ident.demux == "SFMC_1946"] <- "SFMC_RA1946"
SO.sum$ident.sample[SO.sum$ident.demux == "SFMC_3098"] <- "SFMC_RA3098"
SO.sum$ident.sample[SO.sum$ident.demux == "SFMC_5210"] <- "SFMC_RA5210"
SO.sum$ident.sample[SO.sum$ident.demux == "SFMC_3625"] <- "SFMC_RA3625"

SO.sum$ident.sample[SO.sum$orig.ident == "SF1"] <- "SF1"
SO.sum$ident.sample[SO.sum$orig.ident == "SF2"] <- "SF2"
SO.sum$ident.sample[SO.sum$orig.ident == "PB"] <- "PB"

# ident.name
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "PBMC_AS517"] <- "AS_1"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "SFMC_SF433"] <- "AS_1"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "PBMC_AS426"] <- "AS_2"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "SFMC_SF308"] <- "AS_2"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "PBMC_AS210"] <- "AS_3"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "SFMC_SF355"] <- "AS_3"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "PBMC_AS509"] <- "AS_4"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "SFMC_SF409"] <- "AS_4"

SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "SFMC_RA1946"] <- "JJS"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "SFMC_RA3098"] <- "JMH"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "SFMC_RA5210"] <- "CYG"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "SFMC_RA3625"] <- "COS"

SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "SF1"] <- "AS_o"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "SF2"] <- "AS_o"
SO.sum.combined$ident.name[SO.sum.combined$ident.sample == "PB"] <- "AS_o"

# Age
SO.sum$age[SO.sum$ident.sample == "PBMC_AS517"] <- "46"
SO.sum$age[SO.sum$ident.sample == "SFMC_SF433"] <- "46"
SO.sum$age[SO.sum$ident.sample == "PBMC_AS426"] <- "52"
SO.sum$age[SO.sum$ident.sample == "SFMC_SF308"] <- "52"
SO.sum$age[SO.sum$ident.sample == "PBMC_AS210"] <- "41"
SO.sum$age[SO.sum$ident.sample == "SFMC_SF355"] <- "41"
SO.sum$age[SO.sum$ident.sample == "PBMC_AS509"] <- "57"
SO.sum$age[SO.sum$ident.sample == "SFMC_SF409"] <- "57"

SO.sum$age[SO.sum$ident.sample == "SFMC_RA1946"] <- "56"
SO.sum$age[SO.sum$ident.sample == "SFMC_RA3098"] <- "41"
SO.sum$age[SO.sum$ident.sample == "SFMC_RA5210"] <- "53"
SO.sum$age[SO.sum$ident.sample == "SFMC_RA3625"] <- "83"


# sex
SO.sum$gender <- as.character(SO.sum$ident.sample)
SO.sum$gender[SO.sum$ident.sample == "PBMC_AS517"] <- "F"
SO.sum$gender[SO.sum$ident.sample == "SFMC_SF433"] <- "F"
SO.sum$gender[SO.sum$ident.sample == "PBMC_AS426"] <- "M"
SO.sum$gender[SO.sum$ident.sample == "SFMC_SF308"] <- "M"
SO.sum$gender[SO.sum$ident.sample == "PBMC_AS210"] <- "M"
SO.sum$gender[SO.sum$ident.sample == "SFMC_SF355"] <- "M"
SO.sum$gender[SO.sum$ident.sample == "PBMC_AS509"] <- "M"
SO.sum$gender[SO.sum$ident.sample == "SFMC_SF409"] <- "M"

SO.sum$gender[SO.sum$ident.sample == "SFMC_RA1946"] <- "F"
SO.sum$gender[SO.sum$ident.sample == "SFMC_RA3098"] <- "F"
SO.sum$gender[SO.sum$ident.sample == "SFMC_RA5210"] <- "F"
SO.sum$gender[SO.sum$ident.sample == "SFMC_RA3625"] <- "F"

# Dz
SO.sum$ident.dz[SO.sum$ident.sample == "PBMC_AS517"] <- "AS"
SO.sum$ident.dz[SO.sum$ident.sample == "SFMC_SF433"] <- "AS"
SO.sum$ident.dz[SO.sum$ident.sample == "PBMC_AS426"] <- "AS"
SO.sum$ident.dz[SO.sum$ident.sample == "SFMC_SF308"] <- "AS"
SO.sum$ident.dz[SO.sum$ident.sample == "PBMC_AS210"] <- "AS"
SO.sum$ident.dz[SO.sum$ident.sample == "SFMC_SF355"] <- "AS"
SO.sum$ident.dz[SO.sum$ident.sample == "PBMC_AS509"] <- "AS"
SO.sum$ident.dz[SO.sum$ident.sample == "SFMC_SF409"] <- "AS"

SO.sum$ident.dz[SO.sum$ident.sample == "SFMC_RA1946"] <- "RA"
SO.sum$ident.dz[SO.sum$ident.sample == "SFMC_RA3098"] <- "RA" # 얘만 sero(+)라고?
SO.sum$ident.dz[SO.sum$ident.sample == "SFMC_RA5210"] <- "RA"
SO.sum$ident.dz[SO.sum$ident.sample == "SFMC_RA3625"] <- "RA"

SO.sum$ident.dz[SO.sum$ident.sample == "SF1"] <- "AS"
SO.sum$ident.dz[SO.sum$ident.sample == "SF2"] <- "AS"
SO.sum$ident.dz[SO.sum$ident.sample == "PB"] <- "AS"

# origin
SO.sum$origin[SO.sum$ident.sample == "PBMC_AS517"] <- "PB"
SO.sum$origin[SO.sum$ident.sample == "SFMC_SF433"] <- "SF"
SO.sum$origin[SO.sum$ident.sample == "PBMC_AS426"] <- "PB"
SO.sum$origin[SO.sum$ident.sample == "SFMC_SF308"] <- "SF"
SO.sum$origin[SO.sum$ident.sample == "PBMC_AS210"] <- "PB"
SO.sum$origin[SO.sum$ident.sample == "SFMC_SF355"] <- "SF"
SO.sum$origin[SO.sum$ident.sample == "PBMC_AS509"] <- "PB"
SO.sum$origin[SO.sum$ident.sample == "SFMC_SF409"] <- "SF"

SO.sum$origin[SO.sum$ident.sample == "SFMC_RA1946"] <- "SF"
SO.sum$origin[SO.sum$ident.sample == "SFMC_RA3098"] <- "SF"
SO.sum$origin[SO.sum$ident.sample == "SFMC_RA5210"] <- "SF"
SO.sum$origin[SO.sum$ident.sample == "SFMC_RA3625"] <- "SF"

SO.sum$origin[SO.sum$ident.sample == "SF1"] <- "SF"
SO.sum$origin[SO.sum$ident.sample == "SF2"] <- "SF"
SO.sum$origin[SO.sum$ident.sample == "PB"] <- "PB"

# ident.study(old vs new)
SO.sum$ident.study[SO.sum$ident.sample == "PBMC_AS517"] <- "new"
SO.sum$ident.study[SO.sum$ident.sample == "SFMC_SF433"] <- "new"
SO.sum$ident.study[SO.sum$ident.sample == "PBMC_AS426"] <- "new"
SO.sum$ident.study[SO.sum$ident.sample == "SFMC_SF308"] <- "new"
SO.sum$ident.study[SO.sum$ident.sample == "PBMC_AS210"] <- "new"
SO.sum$ident.study[SO.sum$ident.sample == "SFMC_SF355"] <- "new"
SO.sum$ident.study[SO.sum$ident.sample == "PBMC_AS509"] <- "new"
SO.sum$ident.study[SO.sum$ident.sample == "SFMC_SF409"] <- "new"

SO.sum$ident.study[SO.sum$ident.sample == "SFMC_RA1946"] <- "new"
SO.sum$ident.study[SO.sum$ident.sample == "SFMC_RA3098"] <- "new"
SO.sum$ident.study[SO.sum$ident.sample == "SFMC_RA5210"] <- "new"
SO.sum$ident.study[SO.sum$ident.sample == "SFMC_RA3625"] <- "new"

SO.sum$ident.study[SO.sum$ident.sample == "SF1"] <- "old"
SO.sum$ident.study[SO.sum$ident.sample == "SF2"] <- "old"
SO.sum$ident.study[SO.sum$ident.sample == "PB"] <- "old"


# CRP
SO.sum$CRP[SO.sum$ident.sample == "PBMC_AS517"] <- "0.65"
SO.sum$CRP[SO.sum$ident.sample == "SFMC_SF433"] <- "0.65"
SO.sum$CRP[SO.sum$ident.sample == "PBMC_AS426"] <- "7.37"
SO.sum$CRP[SO.sum$ident.sample == "SFMC_SF308"] <- "7.37"
SO.sum$CRP[SO.sum$ident.sample == "PBMC_AS210"] <- "4.11"
SO.sum$CRP[SO.sum$ident.sample == "SFMC_SF355"] <- "4.11"
SO.sum$CRP[SO.sum$ident.sample == "PBMC_AS509"] <- "4.25"
SO.sum$CRP[SO.sum$ident.sample == "SFMC_SF409"] <- "4.25"

Idents(SO.new) <- "ident.demux"
SO.new <- subset(SO.new, idents = c("SFMC_1946", "SFMC_3098", "SFMC_3625", "SFMC_5210"), invert = TRUE)
levels(SO.new)


# QC metrics and filter cells
SO.new[["percent.mt"]] <- PercentageFeatureSet(SO.new, pattern = "^MT-")
VlnPlot(SO.new, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "orig.ident")
VlnPlot(SO.new, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0, y.max = 10000,group.by = "singlet.demux") # 왜 확인하누?
VlnPlot(SO.new, features = "percent.mt", ncol = 3, pt.size = 0,group.by = "singlet.demux") # 왜 확인하누?
VlnPlot(SO.new, features = "percent.mt", ncol = 3, pt.size = 0,group.by = "ident.demux")

Idents(SO.new) <- "singlet.demux"
SO.new <- subset(SO.new, idents = c("SNG", "AMB"))
Idents(SO.new) <- "ident.demux"
table(SO.new$singlet.demux)
table(Idents(SO.new), SO.new$singlet.demux)

# QC after Merge
SO.sum <- merge(AS, y = SO.new, add.cell.ids = c("old", "new"), project = "oldnew")

VlnPlot(SO.sum, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "ident.sample")

Idents(SO.new) <- "ident.demux" 
levels(SO.new)
FeatureScatter(SO.sum, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1) +
  ylim(0,15) +
  xlim(0,40000) +
  geom_hline(yintercept = c(0.5,10), linetype = "dashed") 

FeatureScatter(SO.new, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1) + 
  ylim(0,5000) +
  xlim(0,30000) +
  geom_hline(yintercept = c(500, 3300), linetype = "dashed")+
  geom_vline(xintercept = c(1000, 12000), linetype = "dashed")

FeatureScatter(SO.sum, feature1 = "percent.mt", feature2 = "nFeature_RNA", pt.size = 0.1) + 
  ylim(0,10000) +
  xlim(0,15) +
  geom_hline(yintercept = c(500, 5000), linetype = "dashed") +
  geom_vline(xintercept = c(0.5, 10), linetype = "dashed") 

SO.sum.filtered <- subset(SO.sum, nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt <10)
SO.new.filtered <- subset(SO.new, subset = nFeature_RNA > 500 & nFeature_RNA < 2800)

SO.sum <- SO.sum.filtered
rm(SO.sum.filtered)
```


# 
SO.sum$tot.gem <- NA
SO.sum$tot.gem[SO.sum$ident.study == "old"] <- "old"
SO.sum$tot.gem[SO.sum$orig.gem == "PBMC_F46_SFMC_F56_M57"] <- "new46"
SO.sum$tot.gem[SO.sum$orig.gem == "PBMC_M41_SFMC_F53_M52"] <- "new41"
SO.sum$tot.gem[SO.sum$orig.gem == "PBMC_M52_SFMC_F41_F46"] <- "new52"
SO.sum$tot.gem[SO.sum$orig.gem == "PBMC_M57_SFMC_F83_M41"] <- "new57"


# Integration
SO.sum.list <- SplitObject(object = SO.sum, split.by = "ident.study")
SO.sum.list <- lapply(X = SO.sum.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features.AS <- SelectIntegrationFeatures(object.list = SO.sum.list)
SO.sum.anchors <- FindIntegrationAnchors(object.list = SO.sum.list, anchor.features = features.AS) #defaultassay 사용됨
SO.sum.combined <- IntegrateData(anchorset = SO.sum.anchors)


SO.sum <- SO.sum %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData %>%
  RunPCA %>%
  FindNeighbors %>%
  FindClusters %>%
  RunUMAP(dims= 1:20)
DefaultAssay(SO.sum.combined) <- 'integrated'
SO.sum.combined <- SO.sum.combined %>%
  ScaleData %>%
  RunPCA %>%
  FindNeighbors %>%
  FindClusters %>%
  RunUMAP(dims= 1:20)
p3 <- DimPlot(SO.sum, group.by = 'ident.study', reduction = 'umap') + ggtitle('before study')
p4 <- DimPlot(SO.sum.combined, group.by = 'ident.study', reduction = 'umap') + ggtitle('after study')
CombinePlots(list(p3, p4))


# batch comparison
Idents(SO.sum.combined) <- 'ident.study'
test <- subset(SO.sum.combined, idents = 'new')
DimPlot(test, group.by = 'orig.gem', reduction = 'umap')
DimPlot(test, group.by = 'orig.gem', split.by = 'orig.gem',reduction = 'umap') 
DimPlot(test, group.by = 'ident.name', reduction = 'umap')
DimPlot(test, group.by = 'ident.name', split.by = 'ident.name',reduction = 'umap') 


# Harmony
library(harmony)
SO.sum.harmony <- SO.sum.combined %>%
  RunHarmony("ident.name", assay.use="integrated")

# UMAP and clustering with harmonized PCs
SO.sum.harmony <- RunUMAP(SO.sum.harmony, reduction='harmony', dims = 1:35) # 1:30 → 1:35
SO.sum.harmony <- FindNeighbors(SO.sum.harmony, reduction='harmony')
SO.sum.harmony <- FindClusters(SO.sum.harmony, resolution = 0.1) # 0.3 → 0.1
#save
saveRDS(seurat_obj, file=paste0(out_data_dir, 'seurat_object_sct_harmony.rds'))
p5 <- DimPlot(SO.sum.combined, group.by = 'ident.name', reduction = 'umap') + ggtitle('before harmony')
p6 <- DimPlot(SO.sum.harmony, group.by = 'ident.name', reduction = 'umap') + ggtitle('after harmony')
CombinePlots(list(p5, p6))

DimPlot(SO.sum.harmony, reduction = 'umap', label = TRUE)
DimPlot(SO.sum.combined, group.by = 'ident.study', split.by = 'ident.study',reduction = 'umap') + ggtitle('after study')
DimPlot(SO.sum.harmony, group.by = 'ident.study', split.by = 'ident.study',reduction = 'umap') + ggtitle('after harmony')
DimPlot(SO.sum.combined, group.by = 'origin', split.by = 'origin', reduction = 'umap')
DimPlot(SO.sum.combined, group.by = 'ident.sample', split.by = 'ident.sample', reduction = 'umap', ncol = 4)
DimPlot(SO.sum.harmony, group.by = 'origin', split.by = 'origin', reduction = 'umap')
DimPlot(SO.sum.harmony, group.by = 'ident.sample', split.by = 'ident.sample', reduction = 'umap', ncol = 4)

test1 <- RunUMAP(SO.sum.harmony, reduction='harmony', dims = 1:26)
test2 <- RunUMAP(SO.sum.harmony, reduction='harmony', dims = 1:27)
test3 <- RunUMAP(SO.sum.harmony, reduction='harmony', dims = 1:30)
test4 <- RunUMAP(SO.sum.harmony, reduction='harmony', dims = 1:35)

p1 <- DimPlot(test1, reduction='umap', label=T) + NoLegend()
p2 <- DimPlot(test2, reduction='umap', label=T)+ NoLegend() 
p3 <- DimPlot(test3, reduction='umap', label=T)+ NoLegend()
p4 <- DimPlot(test4, reduction='umap', label=T)+ NoLegend()
CombinePlots(plots = list(p1, p2, p3, p4), ncol = 2)
p3

DimPlot(SO.sum.harmony, reduction = 'umap')

#
DefaultAssay(SO.sum.harmony) <- "integrated" # scaledata, RunPCA, Findclusters는 모두 integrated
test <- ScaleData(SO.sum.harmony, verbose = FALSE) %>% RunPCA() %>% 
  FindNeighbors %>%
  FindClusters(, resolution = c((1:6)/10)) %>%
  RunUMAP(dims = 1:20)

test.res <- test4 %>%
  FindNeighbors %>%
  FindClusters(resolution = c((1:6)/10))
p1 <- DimPlot(test.res, reduction='umap', group.by = "integrated_snn_res.0.1", label=T) + NoLegend()
p2 <- DimPlot(test.res, reduction='umap', group.by = "integrated_snn_res.0.2", label=T)+ NoLegend() 
p3 <- DimPlot(test.res, reduction='umap', group.by = "integrated_snn_res.0.3", label=T)+ NoLegend()
p4 <- DimPlot(test.res, reduction='umap', group.by = "integrated_snn_res.0.4", label=T)+ NoLegend()
p5 <- DimPlot(test.res, reduction='umap', group.by = "integrated_snn_res.0.5", label=T)+ NoLegend()
p6 <- DimPlot(test.res, reduction='umap', group.by = "integrated_snn_res.0.6", label=T)+ NoLegend()

CombinePlots(plots = list(p1, p2, p3, p4, p5, p6), ncol = 3)


# Confirm
cluster_var <- 'seurat_clusters'
clusters <- unique(SO.sum.harmony@meta.data[[cluster_var]])
clusters <- clusters[order(clusters)]
df <- data.frame()
for(i in 1:length(clusters)){
  cur_df <- as.data.frame(SO.sum.harmony@meta.data %>% subset(seurat_clusters == clusters[i]) %>% .$ident.name %>% table() /
                            table(SO.sum.harmony$ident.name))
  
  cur_df$Freq <- cur_df$Freq * 1/(sum(cur_df$Freq))
  
  cur_df$cluster <- clusters[i]
  df <- rbind(df, cur_df)
  print(i)
}

pdf(paste0(fig_dir, 'barplot_batches_by_cluster.pdf'), height=4, width=8)
png(paste0(fig_dir, 'pngs/barplot_batches_by_cluster.png'), height=4, width=8, res=250, units='in')
p <- ggplot(df, aes(y=Freq, x=cluster, fill=.)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  ylab('normalized proportion') +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank()
  )
print(p)

saveRDS(seurat_obj, file=paste0(out_data_dir, 'seurat_object_sct_harmony.rds'))


```
### Cluster annotation
rownames(AS[["FB"]])
DefaultAssay(AS) <- "FB"
# Note that the following command is an alternative but returns the same result
SO.sum.harmony<- NormalizeData(SO.sum, normalization.method = "CLR", margin = 2, assay = "FB")
test <- SO.sum.harmony
DefaultAssay(test) <- "FB"
rownames(test[["FB"]])
test <- NormalizeData(test.SO, normalization.method = "CLR", margin = 2, assay = "FB")
DefaultAssay(SO.sum.harmony) <- "FB"

FeaturePlot(.harmony, "CD14-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD14 protein")
DefaultAssay(SO.T.NK) <- "FB"
Idents(test.CD4.harmony) <- 'ident.study'
test <- subset(test.CD4.harmony, idents = 'old')
DefaultAssay(test) <- "FB"
p1 <- FeaturePlot(test, "CD14-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD14 protein")
p2 <- FeaturePlot(test, "CD19-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD19 protein")
p3 <-FeaturePlot(test, "CD4-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD4 protein")
p4 <- FeaturePlot(test, "CD8-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD8 protein")
p5 <- FeaturePlot(test, "CD56-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD56 protein")
p6 <- FeaturePlot(test, "TCRalphabeta-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("TCRab protein")
CombinePlots(list(p1, p2, p3, p4, p5, p6), ncol = 3)
p1|p2|p3

# FeaturePlot after subset
Idents(test4) <- "ident.study"
test.SO <- subset(test4, idents = "old")
Idents(test.SO) <- "FB"
FeaturePlot(test.SO, "CD4-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD4 protein")


# Now, we will visualize CD14 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(SO.sum.harmony) <- "RNA"
p1 <- FeaturePlot(SO.sum.harmony, "CD19") + ggtitle("CD19 RNA")
DefaultAssay(SO.sum.harmony) <- "FB"
p2 <- FeaturePlot(SO.sum.harmony, "CD19-TotalSeqC", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
p1|p2

# Now, we can include the key in the feature name, which overrides the default assay
p1 <- FeaturePlot(SO.sum.harmony, "adt_CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
p2 <- FeaturePlot(SO.sum.harmony, "rna_CD19") + ggtitle("CD19 RNA")
p1 | p2


# Annotation 
DefaultAssay(SO.sum.harmony) <- "RNA"
MtP.1 <- c("CD14", "HLA-DRB1", "CD3E", "NCAM1", "CD19", "ITGAX", "TIFAB", "HBB", "PPBP")
FeaturePlot(SO.sum.harmony, features = "percent.mt", min.cutoff = "q9")
FeaturePlot(SO.sum.harmony, features = MtP.mitotic, min.cutoff = "q9")

NvM.data <- DotPlot(MM.u, features = MtP.NvM, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()
mitotic.data <- DotPlot(SO.sum.harmony, features = MtP.mitotic, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis() 
DotPlot(SO.sum.harmony, features = c("CD14", "HLA-DRB1", "CD3E", "NCAM1", "CD19", "ITGAX", "TIFAB", "HBB", "PPBP"),
        cols = c("lightgrey", "red"),split.by = 'ident.study', dot.scale = 8)+
  RotatedAxis()
DotPlot(SO.sum.harmony, features = tutorial.mrk,
        cols = c("lightgrey", "red"), dot.scale = 8)+
  RotatedAxis()


MtP.mitotic <- c("TCF7", "CDKN1A")
FeaturePlot(MM.u, features = "percent.mt")

goi = c("CD3D","CD4","KRT7",
        "PTPRC","CD1A","DNTT",
        "CD4","CD3E","GZMB",
        "CD68","IFNG", "CD79A")
FeaturePlot(AS,features = goi,ncol=3)
FeaturePlot(SO.sum.harmony, features = goi, ncol = 3)

Idents(SO.sum.harmony) <- "integrated_snn_res.0.1"
new.cluster.ids <- c("T", "SF Monocyte", "PB Mononcyte", "NK", "DC", "B", "Uncategorized", "pDC", "activated T", "Mk")
names(new.cluster.ids) <- levels(SO.sum.harmony)
SO.sum.harmony <- RenameIdents(SO.sum.harmony, new.cluster.ids)
DimPlot(SO.sum.harmony, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(SO.sum.harmony, reduction = "umap", pt.size = 0.5)


DotPlot(immune.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
  RotatedAxis()

#
findallmarkers_mc <- function(seurat_obj, ...){
  n_clust <- unique(Idents(seurat_obj))
  mcFindMarkers <- function(i){
    ident1 <- i
    ident2 <- n_clust[n_clust != i]
    table <- FindMarkers(seurat_obj,
                         ident.1 = ident1, ident.2 = ident2, ...)
    table$Gene.name.uniq <- rownames(table)
    table$cluster <- rep(i, nrow(table))
    return(table)
  }
  marker_results <- list()[n_clust]
  ptm <- proc.time()
  marker_results <- parallel::mclapply(n_clust, mcFindMarkers, mc.cores = 16)
  time_diff <- proc.time() - ptm
  time_diff
  marker_results
}

findallmarkers_mc(MM.SCT.combined)
MM.markers <- findallmarkers_mc(MM.SCT.combined)


# T/NK
DefaultAssay(SO.sum.harmony) <- 'RNA'
SO.T.NK <- subset(SO.sum.harmony, idents = c("0", "3", "9"))
Idents(SO.sum.harmony)
test <- subset(SO.sum.harmony, idents = c("T", "NK", "activated T"))
test <- SO.sum.harmony[,SO.sum.harmony$integrated_snn_res.0.1 %in% c(0, 3, 9)]


DefaultAssay(SO.T.NK) <- 'integrated'
SO.T.NK <- SO.T.NK %>%
  FindVariableFeatures %>%
  ScaleData %>%
  RunPCA  %>%
  FindNeighbors(reduction = 'harmony') %>%
  FindClusters %>%
  RunUMAP(reduction = 'harmony',dims= 1:7)

Idents(SO.T.NK) <- "seurat_clusters"
DimPlot(SO.T.NK, reduction = 'umap', label= T)

DefaultAssay(test) <- 'RNA'
test.list <- SplitObject(object = test, split.by = "ident.study")
test.list <- lapply(X = test.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features.test <- SelectIntegrationFeatures(object.list = test.list)
test.anchors <- FindIntegrationAnchors(object.list = test.list, anchor.features = features.test) #defaultassay 사용됨
test.combined <- IntegrateData(anchorset = test.anchors)


DefaultAssay(test.combined) <- 'integrated'
test.combined <- test.combined %>%
  ScaleData %>%
  RunPCA  %>%
  FindNeighbors(reduction = 'pca') %>%
  FindClusters %>%
  RunUMAP(reduction = 'pca',dims= 1:7)

DimPlot(test.combined, reduction = 'umap', label = T)

library(harmony)
test.harmony <- test.combined %>%
  RunHarmony("ident.name", assay.use="integrated")


test3 <- FindClusters(test3, resolution = 0.6)
DimPlot(test3, label = T, reduction = 'umap')

SO.T_NK.harmony <- test.harmony

Idents(SO.T.NK) <- 'ident.study'
test1 <- subset(SO.T.NK, idents = 'old')
DefaultAssay(test1) <- "FB"
p1 <- FeaturePlot(test1, "TCRalphabeta-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("TCRab protein")
p2 <- FeaturePlot(test1, "CD4-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD4 protein")
p3 <-FeaturePlot(test1, "CD8-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD8 protein")
p4 <- FeaturePlot(test1, "CD56-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD56 protein")
p5 <- FeaturePlot(test1, "CCR7-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CCR7 protein")
p6 <- FeaturePlot(test1, "CD45RO-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD45RO protein")
p7 <- FeaturePlot(test1, "CD25-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD25 protein")
p8 <- FeaturePlot(test1, "CD69-TotalSeqC", cols = c("lightgrey", "darkgreen")) +
  ggtitle("CD69 protein")
CombinePlots(list(p1, p2, p3, p4, p5, p6, p7, p8), ncol = 4)

DefaultAssay(SO.T.NK) <- 'RNA'
goi2 = c("CD3D","CD3E","CD4","CD8A","NCAM1","FCGR3A","CCR7","CD25","FOXP3")
goi4 = c("GPR65","TNFRSF4","TNFRSF18","CCR6","FOXP3","TBX21","CXCR3","CCR7")
FeaturePlot(SO.T.NK,features = goi4,ncol=4)
FeaturePlot(test3,features = goi2,ncol=4)
DefaultAssay(test) <- 'RNA'
FeaturePlot(test,features = goi2,ncol=4)
DimPlot(test, group.by = 'origin', split.by = 'origin', reduction = 'umap')

print(SO.T.NK[["pca"]], dims = 1:5, nfeatures = 5)

Idents(SO.T.NK) <- 'integrated_snn_res.0.4'

goi2 = c("CD3D","CD3E","CD4","CD8A","NCAM1","FCGR3A","CCR7","CD25","FOXP3")


# How many cells are in each cluster
table(Idents(MM.SCT.combined))
prop.table(table(Idents(SO.sum.harmony)))

install.packages("plotrix")
library(plotrix)
pie3D(x=table(Idents(SO.sum.harmony)), 
      labels=new.cluster.ids, init.angle = 90, clockwise = T, explode = 0.1,
      main = "Cluster proportion", col = rainbow(length(table(Idents(SO.sum.harmony)))))

piepercent<- round(100 * table(Idents(SO.sum.harmony)) / sum(table(Idents(SO.sum.harmony))), 1)
pie(x=table(Idents(SO.sum.harmony)), labels = piepercent,
    main = "Cluster proportion", col = rainbow(length(table(Idents(SO.sum.harmony)))), clockwise = T)
legend("right",new.cluster.ids,
       cex = 0.5, fill = rainbow(length(table(Idents(SO.sum.harmony)))))
view(piepercent)

Idents(test)
table(test$singlet.demux)
table(Idents(test), test$singlet.demux)



## CD4 T cells
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



