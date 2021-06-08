## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111803&fbclid=IwAR3cxRWdm6XLcMiYnWwjhKsy_d6YWtUdXEH9X1bsjDS18rR4ET_mrxuEnt8
## GSE111803: 高通量測序分析肺癌患者和健康對照的外泌體中的miRNA的變化

#00# 安裝套件(若套件已經安裝過可略此步驟)

## 安裝bioconductor及其packages
install.packages("BiocManager")
BiocManager::install() ## 安裝bioconductor
BiocManager::install(c("GEOquery", "DESeq2"))

## 安裝非Bioconductor之R-packages
install.packages(c("pheatmap", "ggfortify", "openxlsx"))

#01# 讀取需要用到的packages

library(GEOquery)  ## 下載GEO上的檔案
library(DESeq2)    ## Count Based RNAseq analysis
library(ggfortify) ## 繪製PCA圖
library(pheatmap)  ## 繪製heat map
library(openxlsx)  ## 讀取、匯出excel檔案

#02# 資料下載、讀取、前處理

## 透過"GEOquery"下載GSE上的raw data
getGEOSuppFiles("GSE111803") 
## 讀取rawdata
GSE111803 <- read.delim(gzfile("GSE111803/GSE111803_Readcount_TPM.txt.gz"))

## 對資料進行前處理
rownames(GSE111803) <- GSE111803[, 1]
GSE111803 <- GSE111803[, -1][, 1:10] ## 只選取read count原始值，而不使用經過TPM校正過的值
colnames(GSE111803) <- sapply(colnames(GSE111803), function(a){strsplit(a, "\\.")[[1]][1]})

## 根據GEO官網上的資料，將病患分成兩類
## (正常:control, 肺腺癌患者:luad)
label <- data.frame(type = c(rep('luad', 5), rep('control', 5)), 
                    row.names = colnames(GSE111803))

#03# 分析、篩選

## 使用DESeq2進行Count Based RNA-seq analysis
GSE111803_ds <- DESeqDataSetFromMatrix(countData = GSE111803, ## 格式轉換
                                       colData = label,
                                       design = ~ type)
GSE111803_ds <- DESeq(GSE111803_ds) ## 標準化
GSE111803_rs <- results(GSE111803_ds, alpha = 0.05) ## 計算結果
plotMA(GSE111803_rs) ## 繪製MA plot

## 根據p value(<0.05)和 log2 fold change(>1)篩選miRNA
GSE111803_rs <- subset(data.frame(GSE111803_rs), pvalue < 0.05)      ## p value(<0.05)
GSE111803_rs <- GSE111803_rs[abs(GSE111803_rs$log2FoldChange) > 1, ] ## log2 fold change(>1)

## Normalized counts transformation
GSE111803_nm <- assay(normTransform(GSE111803_ds))[row.names(GSE111803_rs), ]

#04# 透過圖表檢視分析結果

## 透過heatmap將數據視覺化
GSE111803_nm_zs <- t(apply(GSE111803_nm, 1, function(x){ ## 對每個miRNA做z-score
  return((x - mean(x, na.rm = F)) / sd(x, na.rm = F))
}))
pheatmap(GSE111803_nm_zs, angle_col = 0, cluster_rows = T, show_rownames = T,
         cluster_cols = T, annotation_col = label)

## PCA: 主成分分析
pca <- data.frame(type = label, t(GSE111803_nm), stringsAsFactors = F)
autoplot(prcomp(pca[-1]), data = pca, colour = 'type', label = T, label.size = 3)

#05# 匯出分析結果+圖表

## 創建儲存結果的資料夾
dir.create("data_result", showWarnings = F)

## heatmap匯出成pdf檔
dev.off()
pdf("data_result/heatmap.pdf")
pheatmap(GSE111803_nm_zs, angle_col = 0, cluster_rows = T, show_rownames = T,
         cluster_cols = T, annotation_col = label)
dev.off()
pdf("data_result/pca.pdf")
autoplot(prcomp(pca[-1]), data = pca, colour = 'type', label = T, label.size = 3)
dev.off()

## 將計算結果匯出成excel檔
result <- data.frame("miRNA" = rownames(GSE111803_nm), 
                     GSE111803_nm, 
                     "DESeq2_result" = rep("", nrow(GSE111803_nm)), 
                     GSE111803_rs, 
                     stringsAsFactors = F)
write.xlsx(result, 
           "data_result/GSE111803_result(normalized).xlsx")
