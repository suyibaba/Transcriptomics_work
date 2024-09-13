control1 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/counts8.txt", header = T) 
control2 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/counts9.txt",header = T) 
### subset the df to have only the gene id and counts as column
a = c(2:6)
control1_new = control1[,-(a)]
View(control1_new)
control2_new = control2[,-(a)]
View(control2_new)

########
###load libraries
library(tidyr)
library(tidyverse)
library(dplyr)

###Rename alignment.sorted.bam to control_1
names(control1_new)[names(control1_new) == "alignment.sorted.bam"] <- "control_1"
View(control1_new)
dim(control1_new)
#######Rename alignment.sorted.bam to control_2
names(control2_new)[names(control2_new) == "alignment.sorted.bam"] <- "control_2"
View(control2_new)
dim(control2_new)

#####
##### Remove duplicates and sum the counts of these duplicates from control1 sample
control1_new_df = aggregate(control_1~Geneid, data = control1_new, FUN = sum)
dim(control1_new_df)
View(control1_new_df)


##### Remove duplicates and sum the counts of these duplicates from control2 sample
control2_new_df = aggregate(control_2~Geneid, data = control2_new, FUN = sum)
dim(control2_new_df)
View(control2_new_df)

#Merge control 1 and control 2 by common column
mergedfc_1 = merge(control2_new_df, control1_new_df, by = "Geneid")
dim(mergedfc_1)

###############
control3 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/counts10.txt", header = T)
control4 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/counts.txt", header = T) 
### subset the df to have only the gene id and counts as column
a = c(2:6)
control3_new = control3[,-(a)]
View(control3_new)
control4_new = control4[,-(a)]
View(control4_new)

############
###Rename alignment.sorted.bam to control_3
names(control3_new)[names(control3_new) == "alignment.sorted.bam"] <- "control_3"
View(control3_new)
dim(control3_new)

names(control4_new)[names(control4_new) == "alignment.sorted.bam"] <- "control_4"
View(control4_new)
dim(control4_new)

##### Remove duplicates and sum the counts of these duplicates from control1 sample
control3_new_df = aggregate(control_3~Geneid, data = control3_new, FUN = sum)
dim(control3_new_df)
View(control3_new_df)

##### Remove duplicates and sum the counts of these duplicates from control1 sample
control4_new_df = aggregate(control_4~Geneid, data = control4_new, FUN = sum)
dim(control4_new_df)
View(control4_new_df)

######## #Merge control 1 and control 2 by common column
mergedfc_1s = merge(control3_new_df, control4_new_df, by = "Geneid")
dim(mergedfc_1s)
View(mergedfc_1s)

#Merge control 1 and control 2 by common column
mergedfc_2 = merge(mergedfc_1,mergedfc_1s, by = "Geneid")
dim(mergedfc_2)
View(mergedfc_2)

#####
#control5 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/Nasa_paper/counts.txt", header = T) 
### subset the df to have only the gene id and counts as column
a =  c(2:6)
control5_new = control5[,-(a)]
View(control5_new)
control5_new = control5[,-(a)]
View(control5_new)


###Rename alignment.sorted.bam to control_3
names(control5_new)[names(control5_new) == "alignment.sorted.bam"] <- "control_5"
View(control5_new)
dim(control5_new)

##### Remove duplicates and sum the counts of these duplicates from control1 sample
control5_new_df = aggregate(control_5~Geneid, data = control5_new, FUN = sum)
dim(control5_new_df)
View(control5_new_df)

#Merge control 1 and control 2 by common column
mergedfc_2s = merge(mergedfc_2, control5_new_df, by = "Geneid")
dim(mergedfc_2s)
View(mergedfc_2s)






#####flight
flight1 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/counts4.txt", header = T)  
flight2 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/counts5.txt", header = T) 


### subset the df to have only the gene id and counts as column
a = c(2:6)
flight1_new = flight1[,-(a)]
View(flight1_new)
flight2_new = flight2[,-(a)]
View(flight2_new)

###Rename alignment.sorted.bam to fligth_1
names(flight1_new)[names(flight1_new) == "alignment.sorted.bam"] <- "flight_1"
View(flight1_new)
dim(flight1_new)
#######Rename alignment.sorted.bam to  fligth_2
names(flight2_new)[names(flight2_new) == "alignment.sorted.bam"] <- "flight_2"
View(flight2_new)
dim(flight2_new)

#####
##### Remove duplicates and sum the counts of these duplicates from control1 sample
flight1_new_df = aggregate(flight_1~Geneid, data = flight1_new, FUN = sum)
dim(flight1_new_df)
View(flight1_new_df)


flight2_new_df = aggregate(flight_2~Geneid, data = flight2_new, FUN = sum)
dim(flight2_new_df)
View(flight2_new_df)

#Merge control 1 and control 2 by common column
mergedfc_3 = merge(flight1_new_df,flight2_new_df, by = "Geneid")
dim(mergedfc_3)
View(mergedfc_3)


###########

#####flight
flight3 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/counts6.txt", header = T)  
flight4 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/counts7.txt", header = T) 


### subset the df to have only the gene id and counts as column
a = c(2:6)
flight3_new = flight3[,-(a)]
View(flight3_new)
flight4_new = flight4[,-(a)]
View(flight4_new)

###Rename alignment.sorted.bam to fligth_3
names(flight3_new)[names(flight3_new) == "alignment.sorted.bam"] <- "flight_3"
View(flight3_new)
dim(flight3_new)
#######Rename alignment.sorted.bam to  fligth_4
names(flight4_new)[names(flight4_new) == "alignment.sorted.bam"] <- "flight_4"
View(flight4_new)
dim(flight4_new)

#####
##### Remove duplicates and sum the counts of these duplicates from control1 sample
flight3_new_df = aggregate(flight_3~Geneid, data = flight3_new, FUN = sum)
dim(flight3_new_df)
View(flight3_new_df)


flight4_new_df = aggregate(flight_4~Geneid, data = flight4_new, FUN = sum)
dim(flight4_new_df)
View(flight4_new_df)

#Merge flight3 and flight4 by common column
mergedfc_4 = merge(flight3_new_df,flight4_new_df, by = "Geneid")
dim(mergedfc_4)
View(mergedfc_4)


#Merge flight by common column
mergedfc_5 = merge(mergedfc_3, mergedfc_4,by = "Geneid")
dim(mergedfc_5)
View(mergedfc_5)
mergedfc_2

#Merge n
mergedfc_all = merge(mergedfc_5,mergedfc_2,by = "Geneid")
dim(mergedfc_all)
View(mergedfc_all)

######
mergedfc_all[is.na(mergedfc_all)] <- 0

### Convert the Geneid column tom row_names
rownames(mergedfc_all) <- mergedfc_all$Geneid
View(mergedfc_all)

## remove the Gene id column
merge_fc_all_n = mergedfc_all[,-1]
View(merge_fc_all_n)

library(DESeq2)
##################    
gene_list <- c("Muc1", "Muc2", "Cldn1", "Cldn2", "Tlr4", "St6galnac1", "St3gal2", "St3gal1", "GAL3ST1","B3gnt7","Galnt5", "Spin4",
               "Atoh1","Notch2","Tjp2","Tjp1","Tjp3","Cdx1", "Cdx2", "Hes1","Hes2","Gfl1","Tlr8", "Klf4", "Elf3", "Sox9","Stk11",
               "Mex3a", "Dock4","Mmp9", "Pak1", "Nfkb1", "Nfkb2", "Ptpn11", "Sirt2", "Nfat5",
               "Foxo1","Foxo3","Kcnip3","Lrrc26", "Tgm1", "Tgm2","Tgm3", "Tgm4", "Tgm6","Tgm7", "Clca1",
               "Marcks", "Nlrp6", "Gsdmd", "Myd88", "Il4", "Il13", "Il4i1", "Il4ra",
               "Il33", "Il17a","Il6","Il18","Il1b","Il22ra1", "Il22ra2", "Il25",
               "Il10","Chst5", "Chst6", "Ahr","Agr2", "Atf4", "Hspa5", "Ddit3", "Gcnt1", "Galnt1",
               "Galnt2", "Galnt3","Galnt4", "Galnt5", "Galnt6", "Galnt7", "Large", "B3gnt6",
               "B3gnt6","B3gnt7", "St3gal3", "St3gal4", "St3gal5", "St3gal6", "St8sia1",  "St8sia2",  "Gal3st1", "Adamts1", "Adamts4",
               "Adamts5", "Mmp2", "Mmp14", "Mmp28", "Fut1","Fut2", "Tgfb1", "Tgfb2")

all <- merge_fc_all_n[rownames(merge_fc_all_n) %in% gene_list, ]
dim(all)
View(all)

data = read.csv("/Users/oluwamayowaakinsuyi/Desktop/Nasa_paper/data.csv")
View(data)
###Make the deseq2 object
dds <- DESeqDataSetFromMatrix(countData = merge_fc_all_n,
                              colData = data,
                              design = ~ Group)
#######
dds$Group <- relevel(dds$Group, ref = "ground")
#Perform geometric mean normalization
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
psfdds2 = estimateSizeFactors(dds, geoMeans = geoMeans)
psfdds2

dds2 = DESeq(dds, fitType="local")


##Investigate test results table
res = DESeq(dds2)
res = results(dds2)
res
res_df <- as.data.frame(res)
View(res_df)
dim(res)
res@elementMetadata
res@priorInfo
res_df$genes <- rownames(res_df) 
####################################









write_csv(res_df,"/Users/oluwamayowaakinsuyi/Desktop/rnaseq/deseq_fvsbsl.csv")


####### run edge r
library(edgeR)
dge <- DGEList(counts = all, group = data$Group)

keep = filterByExpr(y = dge)
dge = dge[keep, , keep.lib.sizes=FALSE]

dge$counts

dge =calcNormFactors(object = dge)

dge =estimateDisp(y = dge)

et = exactTest(object = dge)
top_degs = topTags(object = et, n = "Inf")
View(top_degs$table)

dim(merge_fc_all_n)
####wilcox
View(data)
conditions <- data$Group
conditions
gene_names <- rownames(merge_fc_all_n)
wilcox_results <- apply(merge_fc_all_n, 1, function(gene_expr) {
  wilcox.test(gene_expr ~ conditions)$p.value
})

# Store the results in a data frame 
results_df <- data.frame(Gene = gene_names, p_value = wilcox_results)

# Adjust p-values for multiple testing if needed (e.g., using the Benjamini-Hochberg procedure)
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")
View(results_df)
results_df$genes <- rownames(results_df) 
write_csv(results_df,"/Users/oluwamayowaakinsuyi/Desktop/rnaseq/wilcox_f_bsl.csv")


#########
# subset the data for upregulated, downregulated, and not significant genes
upregulated <- subset(res_df, log2FoldChange >= 1 & pvalue < 0.05)
View(upregulated)

upregulated
downregulated <- subset(res_df, log2FoldChange< -1 & pvalue  < 0.05)
View(downregulated)
not_significant <- subset(res_df, !(log2FoldChange > 1 | log2FoldChange < -1) | pvalue > 0.05)
not_significant

library(ggrepel)

#######3

tiff("volcano30days.tiff", units = "in", width = 7, height = 7.5, res = 300)
volcano_plot <- ggplot() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(data = upregulated, aes(x = log2FoldChange, y = -log10(pvalue), color = "Upregulated"), alpha = 0.8) +
  geom_point(data = downregulated, aes(x = log2FoldChange, y = -log10(pvalue), color = "Downregulated"), alpha = 0.8) +
  geom_point(data = not_significant, aes(x = log2FoldChange, y = -log10(pvalue), color = "Not Significant"), alpha = 1, size = 2) +
  labs(x = expression("log"[2]*"Fold Change"), y = "-Log10(p-value)", color = "Gene Status") +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray"), 
                     labels = c("Upregulated", "Downregulated", "Not Significant")) +
  theme_minimal() + 
  theme(
    legend.position = "right",  # Ensure 'bottom' is lowercase
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
  )
volcano_plot
dev.off()
   






