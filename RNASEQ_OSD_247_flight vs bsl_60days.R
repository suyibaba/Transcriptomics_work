control1 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/data_rnaseq/counts1.txt", header = T) 
control2 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/data_rnaseq/counts2.txt",header = T) 
control3 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/data_rnaseq/counts3.txt",header = T)

### subset the df to have only the gene id and counts as column
a = c(2:6)
control1_new = control1[,-(a)]
View(control1_new)
control2_new = control2[,-(a)]
View(control2_new)
control3_new = control3[,-(a)]
View(control3_new)
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
#######
#######Rename alignment.sorted.bam to control_3
names(control3_new)[names(control3_new) == "alignment.sorted.bam"] <- "control_3"
View(control3_new)
dim(control3_new)

#####
##### Remove duplicates and sum the counts of these duplicates from control1 sample
control1_new_df = aggregate(control_1~Geneid, data = control1_new, FUN = sum)
dim(control1_new_df)
View(control1_new_df)


##### Remove duplicates and sum the counts of these duplicates from control2 sample
control2_new_df = aggregate(control_2~Geneid, data = control2_new, FUN = sum)
dim(control2_new_df)
View(control2_new_df)

##################3
##### Remove duplicates and sum the counts of these duplicates from control3 sample
control3_new_df = aggregate(control_3~Geneid, data = control3_new, FUN = sum)
dim(control3_new_df)
View(control3_new_df)


#Merge control 1 and control 2 by common column
mergedfc_1 = merge(control2_new_df, control1_new_df, by = "Geneid")
dim(mergedfc_1)
View(mergedfc_1)
##Merge mergedfc_1 amd ccontrol3_new_df by common column
mergedfc_1s = merge(mergedfc_1, control3_new_df, by = "Geneid")
View(mergedfc_1s)

############### load count 4 to 6
control4 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/data_rnaseq/counts4.txt", header = T)
control5 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/data_rnaseq/counts5.txt", header = T)
control6 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/data_rnaseq/counts6.txt", header = T)

### subset the df to have only the gene id and counts as column
a = c(2:6)
control4_new = control4[,-(a)]
View(control4_new)
control5_new = control5[,-(a)]
View(control5_new)
control6_new = control6[,-(a)]
View(control6_new)
############
###Rename alignment.sorted.bam to control_4
names(control4_new)[names(control4_new) == "alignment.sorted.bam"] <- "control_4"
View(control4_new)
dim(control4_new)
###Rename alignment.sorted.bam to control_5
names(control5_new)[names(control5_new) == "alignment.sorted.bam"] <- "control_5"
View(control5_new)
dim(control5_new)
###Rename alignment.sorted.bam to control_6
names(control6_new)[names(control6_new) == "alignment.sorted.bam"] <- "control_6"
View(control6_new)
dim(control6_new)


##### Remove duplicates and sum the counts of these duplicates from control4 sample
control4_new_df = aggregate(control_4~Geneid, data = control4_new, FUN = sum)
dim(control4_new_df)
View(control4_new_df)

##### Remove duplicates and sum the counts of these duplicates from control5 sample
control5_new_df = aggregate(control_5~Geneid, data = control5_new, FUN = sum)
dim(control5_new_df)
View(control5_new_df)


##### Remove duplicates and sum the counts of these duplicates from control6 sample
control6_new_df = aggregate(control_6~Geneid, data = control6_new, FUN = sum)
dim(control6_new_df)
View(control6_new_df)


######## #Merge control 4 and control 5 by common column
mergedfc_1ss = merge(control4_new_df, control5_new_df, by = "Geneid")
dim(mergedfc_1ss)
View(mergedfc_1ss)

#Merge control 1 and control 2 by common column
mergedfc_2 = merge(mergedfc_1ss,control6_new_df, by = "Geneid")
dim(mergedfc_2)
View(mergedfc_2)


#####
############### load count 7 to 9
control7 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/data_rnaseq/counts7.txt", header = T)
control8 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/data_rnaseq/counts8.txt", header = T)
control9 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/data_rnaseq/counts9.txt", header = T)

### subset the df to have only the gene id and counts as column
a = c(2:6)
control7_new = control7[,-(a)]
View(control7_new)
control8_new = control8[,-(a)]
View(control8_new)
control9_new = control9[,-(a)]
View(control9_new)

############
###Rename alignment.sorted.bam to control_7
names(control7_new)[names(control7_new) == "alignment.sorted.bam"] <- "control_7"
View(control7_new)
dim(control7_new)
###Rename alignment.sorted.bam to control_8
names(control8_new)[names(control8_new) == "alignment.sorted.bam"] <- "control_8"
View(control8_new)
dim(control8_new)
###Rename alignment.sorted.bam to control_9
names(control9_new)[names(control9_new) == "alignment.sorted.bam"] <- "control_9"
View(control9_new)
dim(control9_new)


##### Remove duplicates and sum the counts of these duplicates from control7 sample
control7_new_df = aggregate(control_7~Geneid, data = control7_new, FUN = sum)
dim(control7_new_df)
View(control7_new_df)

##### Remove duplicates and sum the counts of these duplicates from control8 sample
control8_new_df = aggregate(control_8~Geneid, data = control8_new, FUN = sum)
dim(control8_new_df)
View(control8_new_df)


##### Remove duplicates and sum the counts of these duplicates from control9 sample
control9_new_df = aggregate(control_9~Geneid, data = control9_new, FUN = sum)
dim(control9_new_df)
View(control9_new_df)


######## #Merge control 7 and control 8 by common column
mergedfc_2s = merge(control7_new_df, control8_new_df, by = "Geneid")
dim(mergedfc_2s)
View(mergedfc_2s)

#Merge control mergedfc_2ss  and control 9 by common column
mergedfc_2ss = merge(mergedfc_2s,control9_new_df, by = "Geneid")
dim(mergedfc_2ss)
View(mergedfc_2ss)


##### merge all the control 
merge_c = merge(mergedfc_1s ,mergedfc_2, by = "Geneid" )
dim(merge_c)
View(merge_c)

merge_cs = merge(merge_c,mergedfc_2ss, by = "Geneid" )
dim(merge_cs)
View(merge_cs)

######################3
flight1 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/ISS/flight1_counts.txt", header = T) 
flight2 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/ISS/flight2_counts.txt",header = T) 
flight3 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/ISS/flight3_counts.txt",header = T)

View(flight1)
### subset the df to have only the gene id and counts as column
a = c(2:6)
flight1_new = flight1[,-(a)]
View(flight1_new)
flight2_new= flight2[,-(a)]
View(flight2_new)
flight3_new = flight3[,-(a)]
View(flight3_new)
########
###load libraries
library(tidyr)
library(tidyverse)
library(dplyr)

###Rename alignment.sorted.bam to flight1
names(flight1_new)[names(flight1_new) == "alignment.sorted.bam"] <- "flight1"
View(flight1_new)
dim(flight1_new)
#######Rename alignment.sorted.bam to flight2
names(flight2_new)[names(flight2_new) == "alignment.sorted.bam"] <- "flight2"
View(flight2_new)
dim(flight2_new)
#######
#######Rename alignment.sorted.bam to flight3
names(flight3_new)[names(flight3_new) == "alignment.sorted.bam"] <- "flight3"
View(flight3_new)
dim(flight3_new)

#####
##### Remove duplicates and sum the counts of these duplicates from flight1 sample
flight1_new_df = aggregate(flight1~Geneid, data = flight1_new, FUN = sum)
dim(flight1_new_df)
View(flight1_new_df)


##### Remove duplicates and sum the counts of these duplicates from flight1 sample
flight2_new_df = aggregate(flight2~Geneid, data = flight2_new, FUN = sum)
dim(flight2_new_df)
View(flight2_new_df)

##################3
##### Remove duplicates and sum the counts of these duplicates from control3 sample
flight3_new_df = aggregate(flight3~Geneid, data = flight3_new, FUN = sum)
dim(flight2_new_df)
View(flight2_new_df)


#Merge flight1 and flight2 by common column
mergedf1 = merge(flight1_new_df,flight2_new_df, by = "Geneid")
dim(mergedf1)
View(mergedf1)

##Merge mergedfc_1 amd ccontrol3_new_df by common column
mergedfc_1s = merge(mergedf1,flight3_new_df , by = "Geneid")
View(mergedfc_1s)
dim(mergedfc_1s)


####################################### 4 to 6 
flight4 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/ISS/flight4_counts.txt", header = T) 
flight5 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/ISS/flight5_counts.txt",header = T) 
flight6 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/ISS/flight6_counts.txt",header = T)

### subset the df to have only the gene id and counts as column
a = c(2:6)
flight4_new = flight4[,-(a)]
View(flight4_new)
flight5_new= flight5[,-(a)]
View(flight5_new)
flight6_new = flight6[,-(a)]
View(flight6_new)
########
###load libraries
library(tidyr)
library(tidyverse)
library(dplyr)

###Rename alignment.sorted.bam to flight4
names(flight4_new)[names(flight4_new) == "alignment.sorted.bam"] <- "flight4"
View(flight4_new)
dim(flight4_new)
#######Rename alignment.sorted.bam to flight5
names(flight5_new)[names(flight5_new) == "alignment.sorted.bam"] <- "flight5"
View(flight5_new)
dim(flight5_new)
#######
#######Rename alignment.sorted.bam to flight6
names(flight6_new)[names(flight6_new) == "alignment.sorted.bam"] <- "flight6"
View(flight6_new)
dim(flight6_new)



##### Remove duplicates and sum the counts of these duplicates from flight4 sample
flight4_new_df = aggregate(flight4~Geneid, data = flight4_new, FUN = sum)
dim(flight4_new_df)
View(flight4_new_df)


##### Remove duplicates and sum the counts of these duplicates from flight5 sample
flight5_new_df = aggregate(flight5~Geneid, data = flight5_new, FUN = sum)
dim(flight5_new_df)
View(flight5_new_df)

##################3
##### Remove duplicates and sum the counts of these duplicates from control3 sample
flight6_new_df = aggregate(flight6~Geneid, data = flight6_new, FUN = sum)
dim(flight6_new_df)
View(flight6_new_df)

#Merge flight1 and flight2 by common column
mergedf2 = merge(flight4_new_df,flight5_new_df, by = "Geneid")
dim(mergedf2)
View(mergedf2)

##Merge mergedfc_1 amd ccontrol3_new_df by common column
mergedf2s = merge(mergedf2,flight6_new_df , by = "Geneid")
View(mergedf2s)
dim(mergedf2s)

############
####################################### 7 to 9
flight7 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/ISS/flight7_counts.txt", header = T) 
flight8 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/ISS/flight8_counts.txt",header = T) 
flight9 =  read.table("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/ISS/flight9_counts.txt",header = T)

### subset the df to have only the gene id and counts as column
a = c(2:6)
flight7_new = flight7[,-(a)]
View(flight7_new)
flight8_new= flight8[,-(a)]
View(flight8_new)
flight9_new = flight9[,-(a)]
View(flight9_new)
########
###load libraries
library(tidyr)
library(tidyverse)
library(dplyr)

###Rename alignment.sorted.bam to flight7
names(flight7_new)[names(flight7_new) == "alignment.sorted.bam"] <- "flight7"
View(flight7_new)
dim(flight7_new)
#######Rename alignment.sorted.bam to flight8
names(flight8_new)[names(flight8_new) == "alignment.sorted.bam"] <- "flight8"
View(flight8_new)
dim(flight8_new)
#######
#######Rename alignment.sorted.bam to flight9
names(flight9_new)[names(flight9_new) == "alignment.sorted.bam"] <- "flight9"
View(flight9_new)
dim(flight6_new)



##### Remove duplicates and sum the counts of these duplicates from flight7 sample
flight7_new_df = aggregate(flight7~Geneid, data = flight7_new, FUN = sum)
dim(flight7_new_df)
View(flight7_new_df)


##### Remove duplicates and sum the counts of these duplicates from flight8 sample
flight8_new_df = aggregate(flight8~Geneid, data = flight8_new, FUN = sum)
dim(flight8_new_df)
View(flight8_new_df)

##################3
##### Remove duplicates and sum the counts of these duplicates from control3 sample
flight9_new_df = aggregate(flight9~Geneid, data = flight9_new, FUN = sum)
dim(flight9_new_df)
View(flight9_new_df)

#Merge flight1 and flight2 by common column
mergedff = merge(flight7_new_df,flight8_new_df, by = "Geneid")
dim(mergedf2)
View(mergedff)

##Merge mergedfc_1 amd ccontrol3_new_df by common column
mergedffs = merge(mergedff,flight9_new_df , by = "Geneid")
View(mergedffs)

########
mergefffs = merge(mergedf2s, mergedffs, by = "Geneid")
View(mergefffs)

mergefall = merge(mergefffs, mergedfc_1s, by = "Geneid")
  View(mergefall)
 

  ########## Merge flight and control 
  all = merge(mergefall,merge_cs, by = "Geneid")
  View(all)
  
all[is.na(all)] <- 0

### Convert the Geneid column tom row_names
rownames(all) <- all$Geneid
View(all)

## remove the Gene id column
all = all[,-1]
View(all)

library(DESeq2)
##################    
gene_list <- c("Muc1", "Muc2","Muc3", "Cldn1", "Cldn2","Cldn2","Cldn4","Cldn5","Cldn6","Cldn7", "Cldn8", "Tlr4", "St6galnac1", "St3gal2", "St3gal1", "GAL3ST1","B3gnt7","Galnt5", "Spin4",
               "Atoh1","Notch2","Tjp2","Tjp1","Tjp3","Cdx1", "Cdx2", "Hes1","Hes2","Gfl1","Tlr8", "Klf4", "Elf3", "Sox9","Stk11",
               "Mex3a", "Dock4","Mmp9", "Pak1", "Nfkb1", "Nfkb2", "Ptpn11", "Sirt2", "Nfat5",
               "Foxo1","Foxo3","Kcnip3","Lrrc26", "Tgm1", "Tgm2","Tgm3", "Tgm4", "Tgm6","Tgm7", "Clca1",
               "Marcks", "Nlrp6", "Gsdmd", "Myd88", "Il4", "Il13", "Il4i1", "Il4ra",
               "Il33", "Il17a","Il6","Il18","Il1b","Il22ra1", "Il22ra2", "Il25",
               "Il10","Chst5", "Chst6", "Ahr","Agr2", "Atf4", "Hspa5", "Ddit3", "Gcnt1", "Galnt1",
               "Galnt2", "Galnt3","Galnt4", "Galnt5", "Galnt6", "Galnt7", "Large", "B3gnt6",
               "B3gnt6","B3gnt7", "St3gal3", "St3gal4", "St3gal5", "St3gal6", "St8sia1",  "St8sia2",  "Gal3st1", "Adamts1", "Adamts4",
               "Adamts5", "Mmp2", "Mmp14", "Mmp28", "Fut1","Fut2", "Tgfb1", "Tgfb2","Ccl2", "Ccl4", "Cxcl1", "Cdh1","Cxcl8", "Cxcl5", "Lps")

all_n <- all[rownames(all) %in% gene_list, ]
dim(all_n)
View(all_n)

library(DESeq2)
data = read.csv("/Users/oluwamayowaakinsuyi/Desktop/Nasa_paper/datan.csv")
View(data)

###Make the deseq2 object
dds <- DESeqDataSetFromMatrix(countData = all,
                              colData = data,
                              design = ~group)

#Set reference
dds$group<- relevel(dds$group, ref = "control")


#Perform geometric mean normalization
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
View(geoMeans)
psfdds2 = estimateSizeFactors(dds, geoMeans = geoMeans)
psfdds2

dds2 = DESeq(dds, fitType="local")
dds2
#########
geometric_mean_counts <- counts(dds2, normalized = TRUE)
View(geometric_mean_counts)
#######
geometric_mean_counts <- getVarianceStabilizedData(dds2)


##Investigate test results table
res = results(dds2, cooksCutoff = FALSE)
res

res_df <- as.data.frame(res)
View(res_df)
dim(res)
res_df$genes <- rownames(res_df) 
write.csv(res_df,"/Users/oluwamayowaakinsuyi/Desktop/rnaseq/ISS/res_df60days.csv")
res_df = read.csv("/Users/oluwamayowaakinsuyi/Desktop/rnaseq/ISS/res_df60days.csv")


####wilcox test
View(data)
conditions <- data$group
gene_names <- rownames(all)
wilcox_results <- apply(all, 1, function(gene_expr) {
  wilcox.test(gene_expr ~ conditions)$p.value
})

# Store the results in a data frame
results_df <- data.frame(Gene = gene_names, p_value = wilcox_results)

# Adjust p-values for multiple testing if needed (e.g., using the Benjamini-Hochberg procedure)
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")
View(results_df)
results_df$genes <- rownames(results_df) 

write_csv(results_df,"/Users/oluwamayowaakinsuyi/Desktop/rnaseq/wilcox_f_bsl.csv")

######Run Edge r 
dge <- DGEList(counts = all_n, group = data$group)
keep = filterByExpr(y = dge)
dge = dge[keep, , keep.lib.sizes=FALSE]
dge$counts
dge =calcNormFactors(object = dge)
dge =estimateDisp(y = dge)
et = exactTest(object = dge)
top_degs = topTags(object = et, n = "Inf")
View(top_degs$table)

#########
# subset the data for upregulated, downregulated, and not significant genes
upregulated <- subset(res_df, log2FoldChange >= 1 & pvalue < 0.05)
View(upregulated)

upregulated
downregulated <- subset(res_df, log2FoldChange< -1 & pvalue  < 0.05)
View(downregulated)
not_significant <- subset(res_df, !(log2FoldChange > 1 | log2FoldChange < -1) | pvalue > 0.05)
not_significant
View(not_significant)
library(ggrepel)

#######3
tiff("volcano60days.tiff", units = "in", width =5, height = 6, res = 300)
volcano_plot3 <- ggplot() +  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(data = upregulated, aes(x = log2FoldChange, y = -log10(pvalue), color = "Upregulated"), alpha = 0.8) +
  geom_point(data = downregulated, aes(x = log2FoldChange,y = -log10(pvalue), color = "Downregulated"), alpha = 0.8) +
  geom_point(data = not_significant, aes(x = log2FoldChange, y = -log10(pvalue), color = "Not Significant"), alpha = 0.8)  +
  labs(x = expression("Log"[2]*"Fold Change"), y = "-Log10(p-value)", color = "Gene Status") +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal()  +  theme(legend.position = "right", legend.text = element_text(size = 12),
                            legend.title = element_text(size = 14),  axis.title.x = element_text(size = 16),
                            axis.title.y = element_text(size = 16)
  )
volcano_plot3
dev.off()


