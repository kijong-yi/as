Wilcox(03. clustering_annotation에서 배울 게 있음)

# FindallMarkers 할 때 logfc.threshold는 0.0으로도 하시고 0.1로도 하시던데 기준이?

```{r}
CD4_all_markers <- FindAllMarkers(object = CD4_all, logfc.threshold = 0.10)
```
```{r}
write.csv(CD4_all_markers, "CD4_all_markers3.csv")
```

```{r, fig.height=30, fig.width=10}
CD4_all_totalseq_markers <- CD4_all[["FB"]]@counts@Dimnames[[1]] %>% paste0("fb_",.)
FeaturePlot(CD4_all,features = CD4_all_totalseq_markers, ncol=3)
```

```{r fig.height = 6, fig.width = 7}
DimPlot(CD4_all, group.by = "hash.manual", label=T) + NoLegend()
```

```{r fig.height = 10, fig.width = 20}
VlnPlot(CD4_all, features = c("CCR6","RORC","TNFRSF4","TNFRSF18","GPR65","KLRB1"), ncol = 2,pt.size = 0)
```


```{r}
CD4_16_vs_023457_wilcox <- FindMarkers(CD4_all, ident.1 = c("1","6"), ident.2 = c("0","2","3","4","5","7"), test.use = "wilcox", logfc.threshold = 0.0)
```

```{r}
CD4_16_vs_023457_wilcox
```

```{r}
vol
```

```{r}
write.csv(x = CD4_16_vs_023457_wilcox,file = "CD4_16_vs_023457_wilcox_logfc0.csv")
```

```{r}
CD4_all_markers <- FindAllMarkers(object = CD4_all, logfc.threshold = 0.0)
```

```{r}
write.csv(CD4_all_markers, "CD4_all_markers_logfc0.csv")
```


```{r}
CD4_SF <- CD4_all[,CD4_all$seurat_clusters %in%c(1,6)]
```

```{r}
CD4_SF <- CD4_SF %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData %>%
  RunPCA %>%
  FindNeighbors %>%
  FindClusters %>%
  RunUMAP(dims= 1:7)
write_rds(CD4_SF, "data/CD4_SF.Rds", compress = "gz")
```

```{r}
ElbowPlot(CD4_SF)
```

```{r}
DimPlot(CD4_SF)
```

```{r}
DimPlot(CD4_SF,reduction="pca")
```

```{r}
CD4_SF_markers <- FindAllMarkers(object = CD4_SF, logfc.threshold = 0.10)
```

```{r}
write.csv(CD4_SF_markers, "CD4_SF_markers.csv")
```

```{r fig.height = 6, fig.width = 7}
DimPlot(CD4_SF, group.by = "hash.manual", label=T) + NoLegend()
```

```{r}
CD4_hash1 <- CD4_all[,CD4_all$seurat_clusters %in%c(1)]
```

```{r}
CD4_hash1 <- CD4_hash1 %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData %>%
  RunPCA %>%
  FindNeighbors %>%
  FindClusters %>%
  RunUMAP(dims= 1:7)
write_rds(CD4_hash1, "data/CD4_hash1.Rds", compress = "gz")
```

```{r}
ElbowPlot(CD4_all)
```


```{r}
DimPlot(CD4_hash1,reduction="pca")
```


```{r fig.height = 6, fig.width = 7}
DimPlot(CD4_hash1, group.by = "hash.manual", label=T) + NoLegend()
```

```{r}
CD4_hash1_markers <- FindAllMarkers(object = CD4_hash1, logfc.threshold = 0.10)
```

```{r}
write.csv(CD4_hash1_markers, "CD4_hash1_markers.csv")
```


```{r}

Memory_1_vs_2_wilcox <- FindMarkers(Memory, ident.1 = "0", ident.2 = c("1","2"), test.use = "wilcox")

```

```{r}
Memory_1_vs_2_wilcox
```

```{r}
plot(x=Memory_1_vs_2_wilcox$avg_logFC, 
     y=-log10(Memory_1_vs_2_wilcox$p_val), 
     col=adjustcolor("black",alpha=0.2))
```

```{r}
Memory_1_vs_2_wilcox <- FindMarkers(Memory, 
                                    ident.1 = "0", 
                                    ident.2 = c("1","2"), 
                                    test.use = "wilcox",
                                    min.pct = 0,
                                    logfc.threshold = 0,
                                    min.diff.pct = 0)
```

```{r}
plot(x=Memory_1_vs_2_wilcox$avg_logFC, 
     y=-log10(Memory_1_vs_2_wilcox$p_val), 
     col=adjustcolor("black",alpha=0.2))
```

```{r}
write.csv(x = Memory_1_vs_2_wilcox,file = "memory_1_vs_2.csv")
```

```{r}
Memory_markers <- FindAllMarkers(object = Memory, logfc.threshold = 0.10)
```


```{r}
write.csv(Memory_markers, "memory_markers3.csv")
```

### 위 중간 과정에서 특별 gene expression 애들만 찾기
```{r}
ElbowPlot(CD4_hash3)
```

```{r}
Memory <- Lymphocytes[,Lymphocytes$seurat_clusters == "1"]
```


```{r}
Memory <- Memory %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData %>%
  RunPCA %>%
  FindNeighbors %>%
  FindClusters %>%
  RunUMAP(dims= 1:5)
write_rds(Memory, "data/Memory.Rds", compress = "gz")
```

```{r}
ElbowPlot(Memory)
```
```{r}
DimPlot(Memory)
```


```{r}
FeaturePlot(Memory, 
            c("GPR65","TNFRSF4","TNFRSF18","CCR6"))
```



```{r}
Memory$copositive_gpr65_ox40 <- "no"


#Memory$copositive_gpr65_ox40 [Memory@assays$RNA@counts["GPR65",] > 0] <- "gpr"

#Memory$copositive_gpr65_ox40 [Memory@assays$RNA@counts["TNFRSF4",] > 0] <- "ox40"

copos = Memory@assays$RNA@counts["GPR65",] > 0  & Memory@assays$RNA@counts["TNFRSF4",] > 0

Memory$copositive_gpr65_ox40[copos] <- "yes"


DimPlot(Memory, group.by = "copositive_gpr65_ox40")

Memory$copositive_pseudo = "sample(copos)"

DimPlot(Memory, group.by = "copositive_pseudo")
```


```{r}

Memory_1_vs_2_wilcox <- FindMarkers(Memory, ident.1 = "0", ident.2 = c("1","2"), test.use = "wilcox")

```

```{r}
Memory_1_vs_2_wilcox
```

```{r}
plot(x=Memory_1_vs_2_wilcox$avg_logFC, 
     y=-log10(Memory_1_vs_2_wilcox$p_val), 
     col=adjustcolor("black",alpha=0.2))
```

## 어떻게 findallmarkers 함수 사용하는가?
```{r}

AS$manual_cluster <- c(
  "0"="C1",
  "1"="C2",
  "2"="C3",
  "3"="C2",
  "4"="C4",
  "5"="C1",
  "6"="C1",
  "7"="C4",
  "8"="C4",
  "9"="C5",
  "10"="C2",
  "11"="C1",
  "12"="C4",
  "13"="C3",
  "14"="C6",
  "15"="C7",
  "16"="C9",
  "17"="C8",
  "18"="C1",
  "19"="C1",
  "20"="C7",
  "21"="C9",
  "22"="C4",
  "23"="C9"
)[AS$seurat_clusters]

Idents(AS) <- "RNA_snn_res.0.8"
all_markers_seuratcluster <- FindAllMarkers(AS)
all_markers_seuratcluster %>% write_csv("allmarkers_seuratcluster.csv")


Idents(AS) <- "manual_cluster"
all_markers_manualcluster <- FindAllMarkers(AS)
all_markers_manualcluster %>% write_csv("allmarkers_manualcluster.csv")

```



```{r}
all_markers_manualcluster %>% group_by(cluster) %>% slice(1:20) %>% arrange(as.character(cluster)) %>% {colnames(.)=c("p-value", "average log2 fold-change", "% of expressing cells in the cluster", "% of expressing cells not in the cluster", "Adjusted p-value", "Cluster", "Gene Symbol");.} %>% write_csv("marker_top_20genes_per_manual_cluster.csv")
```

```{r}

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
  # nice way to flatten list into a single DF
  # markers <- dplyr::bind_rows(markers)
  marker_results
}


Lymphocytes
Lymphocytes$manual_cluster <- c(
  "0"="CD4, naive",
  "1"="CD4, memory",
  "2"="CD8, memory",
  "3"="CD8, naive",
  "4"="Gamma-Delta/MAIT",
  "5"="NK",
  "6"="CD4, Treg",
  "7"="NK",
  "8"="CD8, memory",
  "9"="CD8, memory",
  "10"="CD8, memory",
  "11"="Rapid proliferating",
  "12"="Uncategorized"
)[Lymphocytes$seurat_clusters]

Idents(Lymphocytes) <- "RNA_snn_res.0.8"
lympho_markers_seuratcluster <- findallmarkers_mc(Lymphocytes)
lympho_markers_seuratcluster %>% write_csv("lympho_markers_seuratcluster.csv")


Idents(Lymphocytes) <- "manual_cluster"
lympho_markers_manualcluster <- findallmarkers_mc(Lymphocytes)
lympho_markers_manualcluster %>% write_csv("lympho_markers_manualcluster.csv")


lympho_markers_manualcluster %>% group_by(cluster) %>% slice(1:20) %>% arrange(as.character(cluster)) %>% {colnames(.)=c("p-value", "average log2 fold-change", "% of expressing cells in the cluster", "% of expressing cells not in the cluster", "Adjusted p-value", "Cluster", "Gene Symbol");.} %>% write_csv("lympho_marker_top_20genes_per_manual_cluster.csv")


##
Key(AS[["FB"]])
DefaultAssay(AS) <- "FB"
FeaturePlot(AS, "CD4", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
featurebarcodes


# batch correction 잘 됐는지 막대그래프로 보는(swaruplab의 harmony에서 발췌)
# new metadata 칸을 만든 후..
cluster_var <- 'seurat_clusters'
clusters <- unique(seurat_obj@meta.data[[cluster_var]])
clusters <- clusters[order(clusters)]
df <- data.frame()
for(i in 1:length(clusters)){
  cur_df <- as.data.frame(seurat_obj@meta.data %>% subset(seurat_clusters == clusters[i]) %>% .$Library.Batch %>% table() /
                            table(seurat_obj$Library.Batch))
  
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
