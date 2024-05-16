#Start with MZmine 2.53 exported files (.csv)
#Set working directory for your files
setwd("C:\\Users\\romme\\OneDrive\\Documents\\test")

#Clean and build your files
x1 <- read.csv("MS1.csv")
x1$row.m.z <- NULL
x1$row.retention.time <- NULL
colnames(x1)[colnames(x1) == "row.ID"] ="sample"
names(x1) <- gsub(x = names(x1), pattern = ".Peak.area", replacement = "")  
headers <- colnames(x1)[-1]
metadata <- data.frame(matrix(unlist(headers), byrow=TRUE),stringsAsFactors=FALSE)
colnames(metadata)[1] = "sample"
metadata[ , 'group'] = ""
write.csv(x1, "matrix.csv", row.names = FALSE)
write.csv(metadata, "metadata.csv", row.names = FALSE)
#Fill the "group" column in the "metadata" file with the corresponding experimental groups
#Save both files as "Text (tab delimited)" in Excel

#1 - Normalization
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("NormalyzerDE")

library(NormalyzerDE)
normalyzer(
  jobName = "Project",
  designPath = "metadata.txt",
  dataPath = "matrix.txt",
  experimentObj = NULL,
  outputDir = ".",
  forceAllMethods = TRUE,
  omitLowAbundSamples = TRUE,
  sampleAbundThres = 5,
  tinyRunThres = 50,
  requireReplicates = TRUE,
  normalizeRetentionTime = TRUE,
  plotRows = 3,
  plotCols = 4,
  zeroToNA = FALSE,
  sampleColName = "sample",
  groupColName = "group",
  inputFormat = "default",
  skipAnalysis = FALSE,
  quiet = FALSE,
  noLogTransform = FALSE,
  writeReportAsPngs = FALSE,
  rtStepSizeMinutes = 1,
  rtWindowMinCount = 100,
  rtWindowShifts = 1,
  rtWindowMergeMethod = "mean"
)
#Move selected normalized file to the working directory

#2 - Imputation using Random Forest algorithm
#install.packages("missForest")
#install.packages("doParallel")
library(missForest)
library(doParallel)
x2 <- as.data.frame(read.table("VSN-normalized.txt", sep="\t", header=TRUE))
x3 <- as.data.frame(apply(x2[, -1], 2, function(x2) 2^x2))
tx3 <- as.data.frame(t(x3))
#set cores for paralellization
registerDoParallel(4)
imputed_results_tx3 <- missForest(tx3, parallelize = "variables")
imp_tx3 <- imputed_results_tx3$ximp
x4 <- as.data.frame(t(imp_tx3))
x4 <- cbind(x2[1], x4)

#3 - Prepare for MetaboAnalyst
#install.packages("dplyr")
library(dplyr)
metadata <- read.csv("metadata.csv")
tmetadata <- as.data.frame(t(metadata))
names <- rownames(tmetadata)
tmetadata <- cbind(names, tmetadata)
colnames(tmetadata) <- tmetadata[1,]
tmetadata <- tmetadata[-1,]
x5 <- rbind(tmetadata, x4)
names(x5) <- gsub(x = names(x5), pattern = "sample", replacement = "Sample")  
x5$Sample[x5$Sample == 'group'] <- 'Label'
#Export file for "Statistical Analysis [one factor]"
write.csv(x5, "MetaboAnalyst.csv", row.names = FALSE)

#4 - Differential Quantification
write.csv(x4, "diffexp.csv", row.names = FALSE)
#Save "diffexp.csv" as "Text (tab delimited)" in Excel
normalyzerDE(
  jobName = "DiffExp",
  comparisons = c("T2DM-Aged"), #Set the groups you want to compare. The second one is the reference.
  designPath = "metadata.txt",
  dataPath = "diffexp.txt",
  experimentObj = NULL,
  outputDir = ".",
  logTrans = TRUE,
  type = "limma",
  sampleCol = "sample",
  condCol = "group",
  batchCol = NULL,
  techRepCol = NULL,
  leastRepCount = 1,
  quiet = FALSE,
  sigThres = 0.05,
  sigThresType = "p",
  log2FoldThres = 0.58,
  writeReportAsPngs = FALSE
)

#5 - Volcano plotting
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("EnhancedVolcano")
diffexp <- as.data.frame(read.delim("DiffExp\\DiffExp_stats.tsv" ,sep="\t"))
#Choose the comparison you want to plot 
names(diffexp) <- gsub(x = names(diffexp), pattern = "T2DM.Aged_PValue", replacement = "pvalue")  
names(diffexp) <- gsub(x = names(diffexp), pattern = "T2DM.Aged_log2FoldChange", replacement = "log2FoldChange")  
diffexp2 = subset(diffexp, select = c("log2FoldChange", "pvalue"))
diffexp2$log2FoldChange <- as.numeric(diffexp2$log2FoldChange)
diffexp2$pvalue <- as.numeric(diffexp2$pvalue)
maxp <- -log10(min(c(diffexp2$pvalue)))
minFC <- abs(min(c(diffexp2$log2FoldChange)))
maxFC <- max(c(diffexp2$log2FoldChange))
up_count <- as.numeric(length(which(diffexp2$log2FoldChange > 0.58 & diffexp2$pvalue < 0.05)))
down_count <- as.numeric(length(which(diffexp2$log2FoldChange < -0.58 & diffexp2$pvalue < 0.05)))
peakcount <- as.numeric(length(diffexp$sample))
if (maxFC > minFC) {
  loglim <- maxFC
} else if (maxFC < minFC) {
  loglim <- minFC
} else {
  loglim <- maxFC
}

#Set colors and cutoff values as needed
keyvals <- ifelse(
  diffexp2$log2FoldChange < -0.58 & diffexp2$pvalue < 0.05, 'blue',
  ifelse(diffexp2$log2FoldChange > 0.58 & diffexp2$pvalue < 0.05, 'red',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'Up'
names(keyvals)[keyvals == 'blue'] <- 'Down'
names(keyvals)[keyvals == 'grey'] <- 'No Change'

#Set "title", "pCutoff", "FCcutoff", and other aesthetics as needed
library(EnhancedVolcano)
volcano <- EnhancedVolcano(diffexp2,
                lab = rownames(diffexp2),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = rownames(diffexp2)[which(names(keyvals) %in% c('Down', 'Up'))],
                xlim = c(-loglim - 0.2, loglim + 0.2),
                ylim = c(0, maxp + 0.5),
                title = 'T2DM vs Aged',
                subtitle = paste0("n = ", peakcount),
                subtitleLabSize = 14,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 3.0,
                labSize = 0,
                colCustom = keyvals,
                cutoffLineType = 'blank',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                colAlpha = 1,
                legendPosition = 'right',
                caption = bquote("Cut-offs:"~Log[2]~"fold change ± 0.58; p < 0.05")
)

plot(volcano) +
  annotate(geom = 'text',
           fontface = 'bold',
           color = 'black',
           size  = 6,
           label = up_count, x =  -loglim*0.5, y = floor(maxp)) +
  annotate(geom = 'text',
           fontface = 'bold',
           color = 'black',
           size  = 6,
           label = down_count, x = loglim*0.5, y = floor(maxp))

           