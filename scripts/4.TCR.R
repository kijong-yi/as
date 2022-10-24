---
  title: "TCR clonotype"
author: "SWG"
date: "2022 10 14"
output: "R_script"
---
  
#
setwd("/home/users/songwoogil/Assf")
colnames(SO.CD4.TCR)

# Loading libraries
getwd()
devtools::install_github("ncborcherding/scRepertoire")
library(scRepertoire)
library(Seurat)

S1 <- read.csv("counts/PBMC_F46_SFMC_F56_M57_TCR/outs/filtered_contig_annotations.csv")
S2 <- read.csv("counts/PBMC_M41_SFMC_F53_M52_TCR/outs/filtered_contig_annotations.csv")
S3 <- read.csv("counts/PBMC_M52_SFMC_F41_F46_TCR/outs/filtered_contig_annotations.csv")
S4 <- read.csv("counts/PBMC_M57_SFMC_F83_M41_TCR/outs/filtered_contig_annotations.csv")

tcr1 <- read.csv("lkj/TCR-IMGT/SFtechrep1-TCR/outs/filtered_contig_annotations.csv")
tcr2 <- read.csv("lkj/TCR-IMGT/SFtechrep2-TCR/outs/filtered_contig_annotations.csv")
tcr3 <- read.csv("lkj/TCR-IMGT/PB-TCR/outs/filtered_contig_annotations.csv")
clono1 <- read.csv("lkj/TCR-IMGT/SFtechrep1-TCR/outs/clonotypes.csv")
clono2 <- read.csv("lkj/TCR-IMGT/SFtechrep2-TCR/outs/clonotypes.csv")
clono3 <- read.csv("lkj/TCR-IMGT/PB-TCR/outs/clonotypes.csv")


# Create a function to trim unwanted "-1" and append sample information before barcodes
barcoder <- function(df, prefix,  trim="\\-1"){
  
  df$barcode <- gsub(trim, "", df$barcode)
  df$barcode <- paste0(prefix, df$barcode)
  
  df
}

# Override barcode data with the suitable replacement
S1 <- barcoder(S1, prefix = "46new_")
S2 <- barcoder(S2, prefix = "41new_")
S3 <- barcoder(S3, prefix = "52new_")
S4 <- barcoder(S4, prefix = "57new_")
## removeNA, removeMulti 불가?


contig_list <- list(S1, S2, S3, S4, S5, S6, S7)

data("contig_list")
head(contig_list[[1]])

SO.CD4.TCR <- SO.CD4.harmony
SO.CD4.harmony <- SO.spare
SO.T_NK.harmony <- SO.spare1
head(colnames(SO.T_NK.TCR@assays[["RNA"]]@counts))
# x : WhichCells(SO.CD4.TCR, idents = "new46")<- paste0("AS1", WhichCells(SO.CD4.TCR, idents = "new46"))
library(stringr) # good
rownames(SO.CD4.harmony@meta.data)
colnames(SO.CD4.harmony)
rownames(SO.CD4.harmony@meta.data)<- ifelse(grepl('-1_2', rownames(SO.CD4.harmony@meta.data)) == T, 
                                            paste0("57", rownames(SO.CD4.harmony@meta.data)), 
                                            rownames(SO.CD4.harmony@meta.data))
colnames(SO.CD4.harmony)<- ifelse(grepl('-1_2_1', colnames(SO.CD4.harmony@assays[["RNA"]]@counts)) == T, 
                                  paste0("41", colnames(SO.CD4.harmony@assays[["RNA"]]@counts)), 
                                  colnames(SO.CD4.harmony@assays[["RNA"]]@counts))
SO.T_NK.harmony@assays[["RNA"]]@counts
rownames(SO.T_NK.harmony@meta.data)<- ifelse(grepl('-1_2', rownames(SO.T_NK.harmony@meta.data)) == T, 
                                             paste0("57", rownames(SO.T_NK.harmony@meta.data)), 
                                             rownames(SO.T_NK.harmony@meta.data))
colnames(SO.T_NK.harmony@assays[["RNA"]]@counts)<- ifelse(grepl('-1_2_1', colnames(SO.T_NK.harmony@assays[["RNA"]]@counts)) == T, 
                                                          paste0("41", colnames(SO.T_NK.harmony@assays[["RNA"]]@counts)), 
                                                          colnames(SO.T_NK.harmony@assays[["RNA"]]@counts))
SO.T_NK.TCR <- SO.T_NK.harmony
SO.CD4.TCR <- SO.CD4.harmony


# 꼬다리 cut 
rownames(SO.T_NK.harmony@meta.data) <- str_replace_all(rownames(SO.T_NK.harmony@meta.data), "-1_1_1", "-1")
rownames(SO.T_NK.harmony@meta.data) <- str_replace_all(rownames(SO.T_NK.harmony@meta.data), "old_", "")
colnames(SO.T_NK.harmony@assays[["RNA"]]@counts) <- str_replace_all(colnames(SO.T_NK.harmony@assays[["RNA"]]@counts), "-1_2_1", "-1")
colnames(SO.T_NK.harmony@assays[["RNA"]]@counts) <- str_replace_all(colnames(SO.T_NK.harmony@assays[["RNA"]]@counts), "old_", "")

tail(rownames(SO.CD4.harmony@meta.data))
tail(rownames(SO.CD4.TCR@meta.data))

SO.CD4.TCR <- SO.CD4.harmony
SO.T_NK.TCR <- SO.T_NK.harmony

DimPlot(Lymphocytes.TCRtyped,reduction = "umap")
DimPlot(CD4_all.TCRtyped)

combined <- combineTCR(contig_list, 
                       sample = c("46new", "41new", "52new", "57new",
                                  "SF1", "SF2", "PB"),
                       removeNA = T,
                       removeMulti = T,
                       cells ="T-AB")

# Barcode processing
tail(rownames(testCD4@meta.data))
v10 <- colnames(SO.spare)
v11 <- colnames(SO.spare1) # T/NK

Idents(SO.spare) <- 'tot.gem'
test <- subset(SO.spare, idents = "new41")
tail(colnames(test))
v11<- ifelse(grepl('-1_1_1', v11), 
             paste0("46", v11), 
             v11)
v11<- ifelse(grepl('-1_2_1', v11), 
             paste0("41", v11), 
             v11)
v11 <- str_replace(v11, "-1_1_1", "-1")
v11 <- str_replace(v11, "-1_2_1", "-1")

v11<- ifelse(grepl('-1_1', v11), 
             paste0("52", v11), 
             v11)
v11<- ifelse(grepl('-1_2', v11), 
             paste0("57", v11), 
             v11)
v11 <- str_replace(v11, "-1_1", "-1")
v11 <- str_replace(v11, "-1_2", "-1")

v11 <- str_replace(v11, "old_", "")

# Combine scRNAseq
testCD4<- combineExpression(combined, test1.3.harmony, 
                            cloneCall="gene", 
                            proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
testT_NK <- combineExpression(combined, test1.3.harmony, 
                              cloneCall="gene", 
                              proportion = FALSE)

colnames(test$tot.gem %in% 'new41')
Idents(test1.3.harmony) <- 'tot.gem'
test <- subset(test1.3.harmony, idents = c("new41", "new46"))
tail(colnames(test))

# color and drawing
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

DimPlot(SO.CD4.TCR, group.by = "cloneType") + NoLegend() +
  scale_color_manual(values=colorblind_vector(2)) + 
  theme(plot.title = element_blank())
slot(SO.CD4.harmony, "meta.data")$cloneType <- factor(slot(SO.CD4.harmony, "meta.data")$cloneType, 
                                                      levels = c("Expanded (>2)","Expanded (2)", "Orphan", "No information"))


DimPlot(SO.CD4.TCR, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="red") + 
  theme(plot.title = element_blank())

##
SO.CD4.TCR <- SO.CD4.harmony
ncol(testCD4)
SO.CD4.harmony$cloneType <- ifelse(SO.CD4.TCR$Frequency!=1,ifelse(SO.CD4.TCR$Frequency>2,"Expanded (>2)","Expanded (2)"),"Orphan")
SO.CD4.harmony$cloneType[is.na(SO.CD4.TCR$cloneType)] <- "No information"
SO.CD4.harmony$cloneType <- factor(SO.CD4.TCR$cloneType, levels=c("Expanded (>2)","Expanded (2)","Orphan","No information"))

SO.T_NK.TCR <- testT_NK
SO.T_NK.TCR$cloneType <- ifelse(SO.T_NK.TCR$Frequency!=1,ifelse(SO.T_NK.TCR$Frequency>2,"Expanded (>2)","Expanded (2)"),"Orphan")
SO.T_NK.TCR$cloneType[is.na(SO.T_NK.TCR$cloneType)] <- "No information"
SO.T_NK.TCR$cloneType <- factor(SO.T_NK.TCR$cloneType, levels=c("Expanded (>2)","Expanded (2)","Orphan","No information"))

color_palette <- c("Expanded (>2)"="red","Expanded (2)"="orange", "Orphan"="grey60","No information"="lightgrey")
color_palette2 <- c("Expanded (>2)"="red","Expanded (2)"="orange", "Orphan"="#eeeeee20","No information"="#11111120")
size_palette <- c("Expanded (>2)"=1,"Expanded (2)"=0.8, "Orphan"=0.3,"No information"=0.3)

DimPlot(SO.CD4.TCR,reduction = "umap",group.by = "cloneType", split.by = 'ident.sample',
        cols=color_palette[SO.CD4.TCR$cloneType],
        pt.size=size_palette[SO.CD4.TCR$cloneType],
        order=c("Expanded (>2)","Expanded (2)")) 

DimPlot(SO.T_NK.TCR,reduction = "umap",group.by = "cloneType",
        cols=color_palette[SO.T_NK.TCR$cloneType],
        pt.size=size_palette[SO.T_NK.TCR$cloneType],
        order=c("Expanded (>2)","Expanded (2)")) + NoLegend()
DimPlot(SO.T_NK.TCR)
colnames(testT_NK@assays[["RNA"]]@counts) <- rownames(testT_NK@meta.data)

Idents(SO.CD4.TCR) <- 'origin'
A <- subset(SO.CD4.TCR, idents = 'SF')
B <- subset(SO.CD4.TCR, idents = 'PB')

a <- table(A$manual_cluster, A$cloneType)
b <- table(B$manual_cluster, B$cloneType)
write.xlsx(a, 'SF.xlsx')
write.xlsx(b, 'PB.xlsx')
library(openxlsx)

pdf("figures/expanded_dimplot.pdf",12,5)
Seurat::DimPlot(SO.CD4.TCR,cols = "red")
v1 = SO.CD4.TCR$cloneType %in% c("Orphan","No information")
c1 = SO.CD4.TCR$ident.study %in% c("old")
table(v1)
table(SO.CD4.TCR@meta.data[["cloneType"]])
plot(SO.CD4.TCR@reductions$umap@cell.embeddings[c1,],
     pch=20, col=color_palette2[SO.CD4.TCR$cloneType][c1])
plot(SO.CD4.TCR@reductions$umap@cell.embeddings[v1,],
     pch=20, col=color_palette2[SO.CD4.TCR$cloneType][v1])
plot(SO.CD4.TCR@reductions$umap@cell.embeddings[v1&&c1,],
     pch=20, col=color_palette2[SO.CD4.TCR$cloneType][v1&&c1])
table(v1&c1)
plot(SO.CD4.TCR@reductions$umap@cell.embeddings[!v1,],
     pch=20, col=color_palette2[SO.CD4.TCR$cloneType][!v1])
points(SO.CD4.TCR@reductions$umap@cell.embeddings[!v1,],
       pch=20, col=color_palette2[SO.CD4.TCR$cloneType][!v1])

v2 = test1.3.harmony$cloneType %in% c("Orphan","No information")
c2 = SO.T_NK.TCR$ident.study %in% c("old")
plot(SO.T_NK.TCR@reductions$umap@cell.embeddings[c2,],
     pch=20, col=color_palette2[SO.T_NK.TCR$cloneType][c2])
plot(SO.T_NK.TCR@reductions$umap@cell.embeddings[v2&&c2,],
     pch=20, col=color_palette2[SO.T_NK.TCR$cloneType][v2&&c2])
plot(test1.3.harmony@reductions$umap@cell.embeddings[v2,],
     pch=20, col=color_palette2[test1.3.harmony$cloneType][v2])
plot(SO.T_NK.TCR@reductions$umap@cell.embeddings[!v2&&c2,],
     pch=20, col=color_palette2[SO.T_NK.TCR$cloneType][!v2&&c2])
table(v2&c2)
points(test1.3.harmony@reductions$umap@cell.embeddings[!v2,],
       pch=20, col=color_palette2[test1.3.harmony$cloneType][!v2])
plot(SO.T_NK.TCR@reductions$umap@cell.embeddings[c2,],
     pch=20, col=color_palette2[SO.T_NK.TCR$cloneType][c2])

c2 = SO.T_NK.TCR$ident.study %in% c("old")


table(SO.CD4.TCR$origin, SO.CD4.TCR$cloneType)
table(SO.CD4.TCR$manual_cluster, SO.CD4.TCR$cloneType)
table(SO.CD4.TCR$seurat_clusters, SO.CD4.TCR$cloneType)
table(SO.CD4.TCR$cop_OX_GI, SO.CD4.TCR$cloneType)
table(SO.CD4.TCR$origin, SO.CD4.TCR$cop_F_OxGi)

c1 <- table(SO.T_NK.TCR$manual_cluster, SO.T_NK.TCR$cloneType)
c2 <-table(SO.CD4.TCR$manual_cluster, SO.CD4.TCR$cloneType)
c3 <- table(SO.CD4.TCR$cop_OX_GI, SO.CD4.TCR$cloneType)
c4 <-table(SO.CD4.TCR$seurat_clusters, SO.CD4.TCR$cloneType)

write.xlsx(c1, 'clone1.xlsx')
write.xlsx(c2, 'clone2.xlsx')
write.xlsx(c3, 'clone3.xlsx')
write.xlsx(c4, 'clone4.xlsx')

levels(SO.CD4.TCR$manual_cluster)
v3 = SO.T_NK.TCR$origin %in% c("SF")
c3 = SO.CD4.TCR$manual_cluster %in% c("SF specific pTh17")
c4 = SO.CD4.TCR$Frequency %in% c("2")
c5 = SO.CD4.TCR$cop_F_OxGi %in% c("2_OXGI")

v4 = SO.T_NK.TCR@assays$RNA@counts %in% ["CD4",] > 0
table(v3 & v4)

plot(SO.CD4.TCR@reductions$umap@cell.embeddings[v3&&c3,],
     pch=20, col=color_palette2[SO.CD4.TCR$cloneType][v3&&c3])
table(v3&c3&c4&c5)

write.xlsx(p3, "t3.xlsx")

cop_FqOxGi <- SO.CD4.TCR$Frequency == '2'  & SO.CD4.TCR$cop_OX_GI == 'OX_GI' 
SO.CD4.TCR$cop_F_OxGi <- "etc."
SO.CD4.TCR$cop_F_OxGi[cop_FqOxGi] <- "2_OXGI"
DimPlot(SO.CD4.TCR, group.by = "cop_F_OxGi",
        pt.size=size_palette1[SO.CD4.TCR$cop_F_OxGi], cols = c("blue", "red"))
DimPlot(SO.CD4.TCR, cols = c("red", "blue"))

cop_FqOxGiG <- SO.CD4.TCR$cop_F_OxGi == "2_OXGI"  & SO.CD4.TCR@assays$RNA@counts["GPR65",] > 0 
SO.CD4.TCR$cop_F_OxGiG <- "etc."
SO.CD4.TCR$cop_F_OxGiG[cop_FqOxGiG] <- "2_OXGIG"
DimPlot(SO.CD4.TCR, group.by = "cop_F_OxGiG",
        pt.size=size_palette2[SO.CD4.TCR$cop_F_OxGiG], cols = c("blue", "red"))

table(SO.CD4.TCR$origin, SO.CD4.TCR$cop_F_OxGiG)
size_palette1 <- c("2_OXGI"=5,"etc."=0.3)
size_palette2<- c("2_OXGIG"=5,"etc."=0.3)

Idents(SO.CD4.TCR) <- 'cop_F_OxGi'
WhichCells(SO.CD4.TCR, idents = '2_OXGI')

