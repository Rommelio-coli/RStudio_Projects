#Include the SIRIUS results inside the working directory
#Set working directory for your files
setwd("C:\\Users\\romme\\OneDrive\\Documents\\test")

canopus <- read.delim("sirius_project\\canopus_formula_summary.tsv", sep="\t", na.strings=c("","NA"))
diffexp <- as.data.frame(read.delim("DiffExp\\DiffExp_stats.tsv" ,sep="\t"))

#1 - Chemical Class Up/Down Plot
#For Chemical Class of ClassyFire
chemicalclass = subset(canopus, select = c("featureId", "ClassyFire.class"))
chemicalclass <- na.omit(chemicalclass)
chemicalclass_diffexp <- merge(chemicalclass, diffexp,
                               by.x=c("featureId"), by.y=c("sample"))
#Splits Up and Down regulated classes
#Indicate comparison
upclasses <- subset(chemicalclass_diffexp[which(chemicalclass_diffexp$T2DM.Aged_PValue < 0.05
                                                & chemicalclass_diffexp$T2DM.Aged_log2FoldChange > 0.58),], 
                    select = c("ClassyFire.class"))
colnames(upclasses) <- c("Class")

downclasses <- subset(chemicalclass_diffexp[which(chemicalclass_diffexp$T2DM.Aged_PValue < 0.05
                                                  & chemicalclass_diffexp$T2DM.Aged_log2FoldChange < -0.58),], 
                      select = c("ClassyFire.class"))
colnames(downclasses) <- c("Class")

#For assigning color to classes
allchanged <- as.data.frame(unique(rbind(upclasses, downclasses)[ , c("Class")]))
colnames(allchanged) <- c("Class")
allchanged[ , 'Color'] = ""
write.csv(allchanged, "ClassColor.csv", row.names = FALSE)
#Fill the table with the colors for plotting the classes

classcolors <- read.csv("ClassColor.csv")

#Frequency tables
up_freq <- as.data.frame(table(c(upclasses$Class)))
down_freq <- as.data.frame(table(c(downclasses$Class)))

#Color and Frequency tables
up_freq <- merge(up_freq, classcolors, by.x=c("Var1"), by.y=c("Class"))
down_freq <- merge(down_freq, classcolors, by.x=c("Var1"), by.y=c("Class"))
up_freq[ , 'Group'] = "T2DM vs Aged Up"
up_freq[ , 'Sum'] = sum(up_freq$Freq)
down_freq[ , 'Group'] = "T2DM vs Aged Down"
down_freq[ , 'Sum'] = sum(down_freq$Freq)
barplot_df <- rbind(up_freq, down_freq)
barplot_df <- with(barplot_df,  barplot_df[order(barplot_df$Var1) , ])
colnames(barplot_df)[colnames(barplot_df) == "Var1"] ="Class"

#Bar Plots
#install.packages("ggplot2")
library(ggplot2)
ggplot(barplot_df, aes(x = Group, y = Freq, fill = Class)) +
  geom_bar(position="fill", stat = "identity", width = 0.8) +
  scale_fill_manual(values = unique(c(barplot_df$Color)))+
  labs(x = "Samples", y = "Relative count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = Sum),
            x = barplot_df$Group, 
            y = 1.01, vjust = 0, hjust = 0.5)
                
#2 - Chemical Diversity Contrast
#Reuse the "metadata.csv" file for grouping
group_table <- read.csv("metadata.csv")
#Reuse the "MetaboAnalyst.csv" file for abundance
abundance <- read.csv("MetaboAnalyst.csv")[-1,]
chemicalclass_abundance <- merge(chemicalclass, abundance,
                               by.x=c("featureId"), by.y=c("Sample"))
