library(tidyverse)
library(data.table)
library(DEseq)
library(ggplot2)

#####DESeq2 jxq data
{
  rm(list = ls())
  iso_anno = fread("GBC_filter_rules_RulesFilter_result_classification.txt")
  iso_anno = filter(iso_anno, filter_result == "Isoform")
  Indir = "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/04_stringtie/prepDE"

  jxq_data = read.csv(paste0(Indir, '/', 'jxq','/', "transcript_count_matrix.csv"))
  jxq_data = merge(jxq_data, iso_anno, by.x = "transcript_id", by.y = "isoform", all.x = T)
  
  jxq_data_group = group_by(jxq_data, associated_gene) %>% 
    summarise(jxq_pair_pt1 = sum(jxq.pair.pt1), jxq_pair_pt10 = sum(jxq.pair.pt10),
              jxq_pair_pt2 = sum(jxq.pair.pt2), jxq_pair_pt3 = sum(jxq.pair.pt3),
              jxq_pair_pt4 = sum(jxq.pair.pt4), jxq_pair_pt5 = sum(jxq.pair.pt5),
              jxq_pair_pt6 = sum(jxq.pair.pt6), jxq_pair_pt7 = sum(jxq.pair.pt7),
              jxq_pair_pt8 = sum(jxq.pair.pt8), jxq_pair_pt9 = sum(jxq.pair.pt9),
              jxq_pair_t1 = sum(jxq.pair.t1), jxq_pair_t10 = sum(jxq.pair.t10),
              jxq_pair_t2 = sum(jxq.pair.t2), jxq_pair_t3 = sum(jxq.pair.t3),
              jxq_pair_t4 = sum(jxq.pair.t4), jxq_pair_t5 = sum(jxq.pair.t5),
              jxq_pair_t6 = sum(jxq.pair.t6), jxq_pair_t7 = sum(jxq.pair.t7),
              jxq_pair_t8 = sum(jxq.pair.t8), jxq_pair_t9 = sum(jxq.pair.t9))
  
  all_data = jxq_data_group
  all_data = all_data[, c(str_which(colnames(all_data), "pair.pt"),
                          str_which(colnames(all_data), "pair.t"))]

  all_data = apply(all_data, 2, as.integer) %>% as.data.frame()
  rownames(all_data) = jxq_data_group$associated_gene
  all_data = all_data[-which(rowSums(all_data) == 0),]
  colnames(all_data) = str_replace_all(colnames(all_data) , '\\.', "_")
  
  group_list = c(rep("pt", ncol(all_data)/2),
                 rep("t",  ncol(all_data)/2))
  group_list <- factor(group_list,levels = c("pt","t"))
  condition = group_list
  subject <- factor(c(1:(ncol(all_data)/2), 1:(ncol(all_data)/2)))
  coldata <- data.frame(row.names = colnames(all_data), condition, subject)
  dds <- DESeqDataSetFromMatrix(countData = all_data,
                                colData = coldata,
                                design = ~ subject + condition)
  
  #dds$condition<- relevel(dds$condition, ref = "pt") 
  dds <- DESeq(dds)
  contrast = c("condition","t","pt")
  res = results(dds, contrast)
  table(res$padj<0.05)
  res$gene= rownames(res)
  write.table(res, paste0( "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/","jxq_gene_DESeq2Paire_DEG.txt"),
              quote = F, sep = "\t", row.names = F)
 
}


###DESeq2  LZN data
{
  rm(list = ls())
  iso_anno = fread("GBC_filter_rules_RulesFilter_result_classification.txt")
  iso_anno = filter(iso_anno, filter_result == "Isoform")
  Indir = "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/04_stringtie/prepDE"
  LZN_data = read.csv(paste0(Indir, '/', 'LZN','/', "transcript_count_matrix.csv"))
  LZN_data = merge(LZN_data, iso_anno, by.x = "transcript_id", by.y = "isoform", all.x = T)
  LZN_data_group = group_by(LZN_data, associated_gene) %>% 
    summarise(LZN_pair_pt1 = sum(LZN.pair.pt1), LZN_pair_pt10 = sum(LZN.pair.pt10),
              LZN_pair_pt11 = sum(LZN.pair.pt11), LZN_pair_pt12 = sum(LZN.pair.pt12),
              LZN_pair_pt13 = sum(LZN.pair.pt13), LZN_pair_pt14 = sum(LZN.pair.pt14),
              LZN_pair_pt15 = sum(LZN.pair.pt15), LZN_pair_pt16  = sum(LZN.pair.pt16 ),
              LZN_pair_pt17 = sum(LZN.pair.pt17), LZN_pair_pt18 = sum(LZN.pair.pt18),
              LZN_pair_pt2 = sum(LZN.pair.pt2), LZN_pair_pt3 = sum(LZN.pair.pt3),
              LZN_pair_pt4 = sum(LZN.pair.pt4), LZN_pair_pt5 = sum(LZN.pair.pt5),
              LZN_pair_pt6 = sum(LZN.pair.pt6), LZN_pair_pt7 = sum(LZN.pair.pt7), 
              LZN_pair_pt7 = sum(LZN.pair.pt7), LZN_pair_pt8 = sum(LZN.pair.pt8), 
              LZN_pair_pt9 = sum(LZN.pair.pt9),
              
              LZN_pair_t1 = sum(LZN.pair.t1), LZN_pair_t10 = sum(LZN.pair.t10),
              LZN_pair_t11 = sum(LZN.pair.t11), LZN_pair_t12 = sum(LZN.pair.t12),
              LZN_pair_t13 = sum(LZN.pair.t13), LZN_pair_t14 = sum(LZN.pair.t14),
              LZN_pair_t15 = sum(LZN.pair.t15), LZN_pair_t16  = sum(LZN.pair.t16 ),
              LZN_pair_t17 = sum(LZN.pair.t17), LZN_pair_t18 = sum(LZN.pair.t18),
              LZN_pair_t2 = sum(LZN.pair.t2), LZN_pair_t3 = sum(LZN.pair.t3),
              LZN_pair_t4 = sum(LZN.pair.t4), LZN_pair_t5 = sum(LZN.pair.t5),
              LZN_pair_t6 = sum(LZN.pair.t6), LZN_pair_t7 = sum(LZN.pair.t7),
              LZN_pair_t7 = sum(LZN.pair.t7), LZN_pair_t8 = sum(LZN.pair.t8),
              LZN_pair_t9 = sum(LZN.pair.t9))
  
  all_data = LZN_data_group
  all_data = all_data[, c(str_which(colnames(all_data), "pair.pt"),
                          str_which(colnames(all_data), "pair.t"))]
  all_data = apply(all_data, 2, as.integer) %>% as.data.frame()
  rownames(all_data) = LZN_data_group$associated_gene
  all_data =  all_data[-which(rowSums(all_data) == 0),]
  colnames(all_data) = str_replace_all(colnames(all_data) , '\\.', "_")
  
  group_list = c(rep("pt", ncol(all_data)/2),
                 rep("t",  ncol(all_data)/2))
  group_list <- factor(group_list,levels = c("pt","t"))
  
  condition = group_list
  subject <- factor(c(1:(ncol(all_data)/2), 1:(ncol(all_data)/2)))
  coldata <- data.frame(row.names = colnames(all_data), condition, subject)
  dds <- DESeqDataSetFromMatrix(countData = all_data,
                                colData = coldata,
                                design = ~ subject + condition)
  #dds$condition<- relevel(dds$condition, ref = "pt") # 指定哪一组作为对照组
  dds <- DESeq(dds)
  contrast = c("condition","t","pt")
  res = results(dds, contrast)
  table(res$padj<0.05)
  res$gene= rownames(res)
  
  write.table(res, paste0( "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/","LZN_gene_DESeq2Paire_DEG.txt"),
              quote = F, sep = "\t", row.names = F)
}

  
##### wj data
{
  rm(list = ls())
  iso_anno = fread("/picb/sysgenomics/gaoli/cell_data/cell_pb_trans_merge/9_squanti_rmPT/GBC_filter_rules_RulesFilter_result_classification.txt")
  iso_anno = filter(iso_anno, filter_result == "Isoform")
  Indir = "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/04_stringtie/prepDE"
  
  wj_data = read.csv(paste0(Indir, '/', 'wj', '/',  "transcript_count_matrix.csv")) 
  wj_data = merge(wj_data, iso_anno, by.x = "transcript_id", by.y = "isoform", all.x = T)
  wj_data_group = group_by(wj_data, associated_gene) %>% 
    summarise(., wj_pair_pt1 = sum(wj.pair.pt1), wj_pair_pt10 = sum(wj.pair.pt10),
              wj_pair_pt2 = sum(wj.pair.pt2), wj_pair_pt3 = sum(wj.pair.pt3),
              wj_pair_pt4 = sum(wj.pair.pt4), wj_pair_pt5 = sum(wj.pair.pt5),
              wj_pair_pt6  = sum(wj.pair.pt6 ), wj_pair_pt7 = sum(wj.pair.pt7),
              wj_pair_pt8 = sum(wj.pair.pt8), wj_pair_pt9 = sum(wj.pair.pt9), 
              
              wj_pair_t1 = sum(wj.pair.t1), wj_pair_t10 = sum(wj.pair.t10),
              wj_pair_t2 = sum(wj.pair.t2), wj_pair_t3 = sum(wj.pair.t3),
              wj_pair_t4 = sum(wj.pair.t4), wj_pair_t5 = sum(wj.pair.t5),
              wj_pair_t6  = sum(wj.pair.t6), wj_pair_t7 = sum(wj.pair.t7),
              wj_pair_t8 = sum(wj.pair.t8), wj_pair_t9 = sum(wj.pair.t9))
              
  all_data = wj_data_group
  all_data = all_data[, c(str_which(colnames(all_data), "pair.pt"),
                          str_which(colnames(all_data), "pair.t"))]
  all_data = apply(all_data, 2, as.integer) %>% as.data.frame()
  rownames(all_data) = wj_data_group$associated_gene
  all_data =  all_data[-which(rowSums(all_data) == 0),]
  colnames(all_data) = str_replace_all(colnames(all_data) , '\\.', "_")
  
  group_list = c(rep("pt", ncol(all_data)/2),
                 rep("t",  ncol(all_data)/2))  
  group_list <- factor(group_list,levels = c("pt","t"))
  condition = group_list
  
  subject <- factor(c(1:(ncol(all_data)/2), 1:(ncol(all_data)/2)))
  coldata <- data.frame(row.names = colnames(all_data), condition, subject)
  
  dds <- DESeqDataSetFromMatrix(countData = all_data,
                                colData = coldata,
                                design = ~ subject + condition)
  
  #dds$condition<- relevel(dds$condition, ref = "pt") # 指定哪一组作为对照组
  dds <- DESeq(dds)
  contrast = c("condition","t","pt")
  res = results(dds, contrast)
  
  table(res$padj<0.05)
  
  res$gene= rownames(res)
  
  write.table(res, paste0( "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/","wj_gene_DESeq2Paire_DEG.txt"),
              quote = F, sep = "\t", row.names = F)
}


####zq data
{
  rm(list = ls())
  iso_anno = fread("/picb/sysgenomics/gaoli/cell_data/cell_pb_trans_merge/9_squanti_rmPT/GBC_filter_rules_RulesFilter_result_classification.txt")
  iso_anno = filter(iso_anno, filter_result == "Isoform")
  
  Indir = "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/04_stringtie/prepDE"
  zq_data = read.csv(paste0(Indir, '/', 'zq', '/', "transcript_count_matrix.csv"))
  zq_data = merge(zq_data, iso_anno, by.x = "transcript_id", by.y = "isoform", all.x = T)
  zq_data_group = group_by(zq_data, associated_gene) %>% 
    summarise(
    zq_pair_pt1 = sum(zq.pair.pt1), zq_pair_pt10 = sum(zq.pair.pt10), zq_pair_pt11 = sum(zq.pair.pt11), zq_pair_pt12 = sum(zq.pair.pt12), 
    zq_pair_pt13 = sum(zq.pair.pt13), zq_pair_pt14 = sum(zq.pair.pt14), zq_pair_pt15 = sum(zq.pair.pt15), zq_pair_pt16 = sum(zq.pair.pt16),
    zq_pair_pt17 = sum(zq.pair.pt17), zq_pair_pt18 = sum(zq.pair.pt18), zq_pair_pt19 = sum(zq.pair.pt19), zq_pair_pt2 = sum(zq.pair.pt2),
    zq_pair_pt20 = sum(zq.pair.pt20), zq_pair_pt21 = sum(zq.pair.pt21), zq_pair_pt22 = sum(zq.pair.pt22), zq_pair_pt23 = sum(zq.pair.pt23),
    zq_pair_pt24 = sum(zq.pair.pt24), zq_pair_pt25 = sum(zq.pair.pt25), zq_pair_pt26 = sum(zq.pair.pt26), zq_pair_pt27 = sum(zq.pair.pt27),
    zq_pair_pt28 = sum(zq.pair.pt28), zq_pair_pt29 = sum(zq.pair.pt29), zq_pair_pt3 = sum(zq.pair.pt3), zq_pair_pt30 = sum(zq.pair.pt30),
    zq_pair_pt31 = sum(zq.pair.pt31), zq_pair_pt32 = sum(zq.pair.pt32), zq_pair_pt33 = sum(zq.pair.pt33), zq_pair_pt34 = sum(zq.pair.pt34),
    zq_pair_pt35 = sum(zq.pair.pt35), zq_pair_pt36 = sum(zq.pair.pt36), zq_pair_pt37 = sum(zq.pair.pt37), zq_pair_pt38 = sum(zq.pair.pt38),        
    zq_pair_pt39 = sum(zq.pair.pt39), zq_pair_pt4 = sum(zq.pair.pt4), zq_pair_pt40 = sum(zq.pair.pt40), zq_pair_pt41 = sum(zq.pair.pt41), 
    zq_pair_pt5 = sum(zq.pair.pt5), zq_pair_pt6 = sum(zq.pair.pt6), zq_pair_pt7 = sum(zq.pair.pt7), zq_pair_pt8 = sum(zq.pair.pt8), 
    zq_pair_pt9 = sum(zq.pair.pt9),
    
    zq_pair_t1 = sum(zq.pair.t1), zq_pair_t10 = sum(zq.pair.t10), zq_pair_t11 = sum(zq.pair.t11), zq_pair_t12 = sum(zq.pair.t12), 
    zq_pair_t13 = sum(zq.pair.t13), zq_pair_t14 = sum(zq.pair.t14), zq_pair_t15 = sum(zq.pair.t15), zq_pair_t16 = sum(zq.pair.t16),
    zq_pair_t17 = sum(zq.pair.t17), zq_pair_t18 = sum(zq.pair.t18), zq_pair_t19 = sum(zq.pair.t19), zq_pair_t2 = sum(zq.pair.t2),
    zq_pair_t20 = sum(zq.pair.t20), zq_pair_t21 = sum(zq.pair.t21), zq_pair_t22 = sum(zq.pair.t22), zq_pair_t23 = sum(zq.pair.t23),
    zq_pair_t24 = sum(zq.pair.t24), zq_pair_t25 = sum(zq.pair.t25), zq_pair_t26 = sum(zq.pair.t26), zq_pair_t27 = sum(zq.pair.t27),
    zq_pair_t28 = sum(zq.pair.t28), zq_pair_t29 = sum(zq.pair.t29), zq_pair_t3 = sum(zq.pair.t3), zq_pair_t30 = sum(zq.pair.t30),
    zq_pair_t31 = sum(zq.pair.t31), zq_pair_t32 = sum(zq.pair.t32), zq_pair_t33 = sum(zq.pair.t33), zq_pair_t34 = sum(zq.pair.t34),
    zq_pair_t35 = sum(zq.pair.t35), zq_pair_t36 = sum(zq.pair.t36), zq_pair_t37 = sum(zq.pair.t37), zq_pair_t38 = sum(zq.pair.t38),        
    zq_pair_t39 = sum(zq.pair.t39), zq_pair_t4 = sum(zq.pair.t4), zq_pair_t40 = sum(zq.pair.t40), zq_pair_t41 = sum(zq.pair.t41), 
    zq_pair_t5 = sum(zq.pair.t5), zq_pair_t6 = sum(zq.pair.t6), zq_pair_t7 = sum(zq.pair.t7), zq_pair_t8 = sum(zq.pair.t8),
    zq_pair_t9 = sum(zq.pair.t9))
  
  all_data = zq_data_group
  all_data = all_data[, c(str_which(colnames(all_data), "pair.pt"),
                          str_which(colnames(all_data), "pair.t"))]
  all_data = apply(all_data, 2, as.integer) %>% as.data.frame()
  rownames(all_data) = zq_data_group$associated_gene
  all_data =  all_data[-which(rowSums(all_data) == 0),]
  colnames(all_data) = str_replace_all(colnames(all_data) , '\\.', "_")
  
  group_list = c(rep("pt", ncol(all_data)/2),
                 rep("t",  ncol(all_data)/2))
  
  group_list <- factor(group_list,levels = c("pt","t"))
  
  condition = group_list
  subject <- factor(c(1:(ncol(all_data)/2), 1:(ncol(all_data)/2)))
  
  coldata <- data.frame(row.names = colnames(all_data), condition, subject)
  
  dds <- DESeqDataSetFromMatrix(countData = all_data,
                                colData = coldata,
                                design = ~ subject + condition)
  
  #dds$condition<- relevel(dds$condition, ref = "pt") # 指定哪一组作为对照组
  
  dds <- DESeq(dds)
  
  contrast = c("condition","t","pt")
  res = results(dds, contrast)
  
  table(res$padj<0.05)
  
  res$gene= rownames(res)
  
  write.table(res, paste0( "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/","zq_gene_DESeq2Paire_DEG.txt"),
              quote = F, sep = "\t", row.names = F)
}


#####分别读入差异分析文件 然后对splicing factor进行标注
zq_data = read_tsv("/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/zq_gene_DESeq2Paire_DEG.txt")
colnames(zq_data) = paste0("zq_",colnames(zq_data))

LZN_data = read_tsv("/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/LZN_gene_DESeq2Paire_DEG.txt")
colnames(LZN_data) = paste0("LZN_",colnames(LZN_data))

wj_data = read_tsv("/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/wj_gene_DESeq2Paire_DEG.txt")
colnames(wj_data) = paste0("wj_", colnames(wj_data))

jxq_data = read_tsv("/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/jxq_gene_DESeq2Paire_DEG.txt")
colnames(jxq_data) = paste0("jxq_", colnames(jxq_data))

all_data = merge(zq_data, LZN_data, by.x = "zq_gene", by.y = "LZN_gene", all.x = T) %>% 
  merge(., wj_data, by.x = "zq_gene", by.y = "wj_gene", all.x = T) %>% 
  merge(., jxq_data, by.x = "zq_gene", by.y = "jxq_gene", all.x = T)

ref_data = import("gencode.v39.annotation.gtf")

ref_data_gene = subset(ref_data, type == "gene") %>% as.data.frame()
sf_gene_name = read_tsv("sf_gene_name.txt")
sf_gene_data = data.frame(sf_gene = sf_gene_name, sf_gene_factor = 1)

all_data_sf = merge(all_data, sf_gene_data, by.x = "zq_gene", by.y = "sf_gene", all.x = T)

write.table(all_data_sf, paste0( "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/","all_gene_DESeq2Paire_DEG_SF.txt"),
            quote = F, sep = "\t", row.names = F)


{
  rm(list = ls())
  all_data_sf = read_tsv(paste0( "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/","all_gene_DESeq2Paire_DEG_SF.txt"))
  colnames(all_data_sf)[1] = "gene" 
  
  zq_sf_data = all_data_sf[,c(str_which(colnames(all_data_sf), "zq"), 1, 26)] %>% as.data.frame() %>% 
    subset(., sf_gene_factor == "1")
  
  zq_sf_data$logP = -log10(zq_sf_data$zq_padj)
  
  zq_sf_data$Group = "not-significant"
  
  zq_sf_data$Group[which(zq_sf_data$zq_padj < 0.05 & zq_sf_data$zq_log2FoldChange > 1)] = "up_regulated"
  zq_sf_data$Group[which(zq_sf_data$zq_padj < 0.05 & zq_sf_data$zq_log2FoldChange < -1)] = "down-regulated"
  
  zq_sf_data$Label = ""
  
  zq_sf_data = zq_sf_data[order(zq_sf_data$zq_padj),]
  up.genes = head(zq_sf_data$gene[which(zq_sf_data$Group == "up_regulated")], 10)
  down.genes = head(zq_sf_data$gene[which(zq_sf_data$Group == "down-regulated")], 10)
  
  topGene = c(as.character(up.genes), as.character(down.genes))
  zq_sf_data$Label[match(topGene, zq_sf_data$gene)] = topGene
  
  p1 = ggscatter(zq_sf_data, x = "zq_log2FoldChange", y = "logP",
                 color = "Group",
                 palette = c("#2f5688", "#BBBBBB", "#CC0000"),
                 size = 4,
                 label = zq_sf_data$Label,
                 font.label = 8,
                 repel = T,
                 xlab = "log2FoldChange",
                 ylab = '-log10(Adjust P-value)',) +
    theme_base() +
    geom_hline(yintercept = 1.3, linetype = "dashed") +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") + ggtitle("zq_data")
  
}


####提取LZN data 并画火山图
{
  #rm(list = ls())
  all_data_sf = read_tsv(paste0( "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/","all_gene_DESeq2Paire_DEG_SF.txt"))
  colnames(all_data_sf)[1] = "gene" 
  
  LZN_sf_data = all_data_sf[,c(str_which(colnames(all_data_sf), "LZN"), 1, 26)] %>% as.data.frame() %>% 
    subset(., sf_gene_factor == "1")
  
  LZN_sf_data$logP = -log10(LZN_sf_data$LZN_padj)
  
  LZN_sf_data$Group = "not-significant"
  
  LZN_sf_data$Group[which(LZN_sf_data$LZN_padj < 0.05 & LZN_sf_data$LZN_log2FoldChange > 1)] = "up_regulated"
  LZN_sf_data$Group[which(LZN_sf_data$LZN_padj < 0.05 & LZN_sf_data$LZN_log2FoldChange < -1)] = "down-regulated"
  
  LZN_sf_data$Label = ""
  
  LZN_sf_data = LZN_sf_data[order(LZN_sf_data$LZN_padj),]
  up.genes = head(LZN_sf_data$gene[which(LZN_sf_data$Group == "up_regulated")], 10)
  down.genes = head(LZN_sf_data$gene[which(LZN_sf_data$Group == "down-regulated")], 10)
  
  topGene = c(as.character(up.genes), as.character(down.genes))
  LZN_sf_data$Label[match(topGene, LZN_sf_data$gene)] = topGene
  
  p2 = ggscatter(LZN_sf_data, x = "LZN_log2FoldChange", y = "logP",
                 color = "Group",
                 palette = c("#2f5688", "#BBBBBB", "#CC0000"),
                 size = 4,
                 label = LZN_sf_data$Label,
                 font.label = 8,
                 repel = T,
                 xlab = "log2FoldChange",
                 ylab = '-log10(Adjust P-value)',) +
    theme_base() +
    geom_hline(yintercept = 1.30, linetype = "dashed") +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") + ggtitle("LZN_data")
  
}


{
  #rm(list = ls())
  all_data_sf = read_tsv(paste0( "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/","all_gene_DESeq2Paire_DEG_SF.txt"))
  colnames(all_data_sf)[1] = "gene" 
  
  wj_sf_data = all_data_sf[,c(str_which(colnames(all_data_sf), "wj"), 1, 26)] %>% as.data.frame() %>% 
    subset(., sf_gene_factor == "1")
  
  wj_sf_data$logP = -log10(wj_sf_data$wj_padj)
  
  wj_sf_data$Group = "not-significant"
  
  wj_sf_data$Group[which(wj_sf_data$wj_padj < 0.05 & wj_sf_data$wj_log2FoldChange > 1)] = "up_regulated"
  wj_sf_data$Group[which(wj_sf_data$wj_padj < 0.05 & wj_sf_data$wj_log2FoldChange < -1)] = "down-regulated"
  
  wj_sf_data$Label = ""
  
  wj_sf_data = wj_sf_data[order(wj_sf_data$wj_padj),]
  up.genes = head(wj_sf_data$gene[which(wj_sf_data$Group == "up_regulated")], 10)
  down.genes = head(wj_sf_data$gene[which(wj_sf_data$Group == "down-regulated")], 10)
  
  topGene = c(as.character(up.genes), as.character(down.genes))
  wj_sf_data$Label[match(topGene, wj_sf_data$gene)] = topGene
  
  p3 = ggscatter(wj_sf_data, x = "wj_log2FoldChange", y = "logP",
                 color = "Group",
                 palette = c("#2f5688", "#BBBBBB", "#CC0000"),
                 size = 4,
                 label = wj_sf_data$Label,
                 font.label = 8,
                 repel = T,
                 xlab = "log2FoldChange",
                 ylab = '-log10(Adjust P-value)',) +
    theme_base() +
    geom_hline(yintercept = 1.30, linetype = "dashed") +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") + ggtitle("wj_data")
  
}

{
  #rm(list = ls())
  all_data_sf = read_tsv(paste0( "/picb/sysgenomics/gaoli/cell_data/linGang_RNAseq/","all_gene_DESeq2Paire_DEG_SF.txt"))
  colnames(all_data_sf)[1] = "gene" 
  
  jxq_sf_data = all_data_sf[,c(str_which(colnames(all_data_sf), "jxq"), 1, 26)] %>% as.data.frame() %>% 
    subset(., sf_gene_factor == "1")
  
  jxq_sf_data$logP = -log10(jxq_sf_data$jxq_padj)
  
  jxq_sf_data$Group = "not-significant"
  
  jxq_sf_data$Group[which(jxq_sf_data$jxq_padj < 0.05 & jxq_sf_data$jxq_log2FoldChange > 1)] = "up_regulated"
  jxq_sf_data$Group[which(jxq_sf_data$jxq_padj < 0.05 & jxq_sf_data$jxq_log2FoldChange < -1)] = "down-regulated"
  
  jxq_sf_data$Label = ""
  
  jxq_sf_data = jxq_sf_data[order(jxq_sf_data$jxq_padj),]
  up.genes = head(jxq_sf_data$gene[which(jxq_sf_data$Group == "up_regulated")], 10)
  down.genes = head(jxq_sf_data$gene[which(jxq_sf_data$Group == "down-regulated")], 10)
  
  topGene = c(as.character(up.genes), as.character(down.genes))
  jxq_sf_data$Label[match(topGene, jxq_sf_data$gene)] = topGene
  
  p4 = ggscatter(jxq_sf_data, x = "jxq_log2FoldChange", y = "logP",
                 color = "Group",
                 palette = c("#2f5688", "#BBBBBB", "#CC0000"),
                 size = 4,
                 label = jxq_sf_data$Label,
                 font.label = 8,
                 repel = T,
                 xlab = "log2FoldChange",
                 ylab = '-log10(Adjust P-value)',) +
    theme_base() +
    geom_hline(yintercept = 1.3, linetype = "dashed") +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") + ggtitle("jxq_data")
  
}
