#################################################################
#
# date: November 26, 2016
# platform: Debian 3.16.74-1
# R.version :  3.6.1
# author: Villemin Jean-Philippe
# team: Machine Learning and Gene Regulation - IGH
#
# Purpose : Differential Expression of CEL Affymetrix 2.1 ST (HG19)
#
# Usage : Do Rscript AnalyseMicroArray.R in the current directory.
# 
# Design : 
# https://www.thermofisher.com/order/catalog/product/902136#/902136
#
# Docs Pipeline : 
# http://homer.ucsd.edu/homer/basicTutorial/affymetrix.html
# https://bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html
# https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
# Bug :
# https://support.bioconductor.org/p/117119/
# https://support.bioconductor.org/p/122925/
#
# Note : You can change the cutoff for fold change line 380.
#################################################################


#BiocManager::install(c("hugene21sttranscriptcluster.db","oligo","pd.hugene.2.1.st","Biobase","RColorBrewer"))
#BiocManager::install(c("Biobase","limma","AnnotationDbi","ArrayExpress","dplyr","tidyr"))
#BiocManager::install(c("arrayQualityMetrics")) # Need to install apt-get install libcairo2-dev

library(hugene21sttranscriptcluster.db)
library(pd.hugene.2.1.st)
library(Biostrings)
library(gtools)

#https://support.bioconductor.org/p/112949/
#The difference between the pdInfo package pd.hugene.2.0.st) and the ChipDb package #hugene20sttranscriptcluster.db) is that the former contains the relatively unprocessed #annotation data from Affymetrix, and the latter contains just the annotations that could be #mapped to NCBI Gene IDs from whatever is in the geneassignment column of the Affymetrix #annotation column.

library(AnnotationDbi)
library(oligo)
library(tidyr)
library(dplyr)
library(ggplot2)
library(limma)
library(Biobase)
library(oligoClasses)
library(arrayQualityMetrics)
library(RColorBrewer)

#################################################################################################################
##########################			Read Input 	################################################################
#################################################################################################################
print("Read Input :")

cels = list.files("data", pattern = "cel",full.names = TRUE)
print(cels)

affyRaw <- oligo::read.celfiles(filenames=cels)

##########################		One Shot to retrieve the design of the prob	#################################
#################################################################################################################
field <- c('fid', 'fsetid', 'level', 'type', 'x', 'y', 'chrom', 'start', 'stop')
pInfo <- getProbeInfo(affyRaw, field=field, sortBy='fid', target='core')
data(pmSequence, package=annotation(pd.hugene.2.1.st))
idx <- match(pInfo[["fid"]], pmSequence[["fid"]])
pmSequence <- pmSequence[idx,]
pmSequence <- as.data.frame(pmSequence)
final_annotated_table_for_probeset <- dplyr::full_join(pmSequence, pInfo, by = "fid")
write.table(final_annotated_table_for_probeset,row.names=FALSE,col.names=TRUE,file="Probes.tsv" ,quote=FALSE)
#head(final_annotated_table_for_probeset)

filenames <- sampleNames(affyRaw)
pData(affyRaw)$filenames <- filenames

sampleNames <- sub(".cel", "", filenames)
sampleNames(affyRaw) <- sampleNames

pData(affyRaw)$group <- ifelse(grepl("Foyer", sampleNames(affyRaw)),"FOYER", "TUMOR")
#pData(affyRaw)$patient <- c("1","1","2","2","3","3","4","4","5","5","6","6","7","7","8","8") # NOT VERY CLEAN
pData(affyRaw)$patient <- c("1","1","2","2","3","3","4","4","5","5","6","6","7","7","8","8") # NOT VERY CLEAN

head(Biobase::pData(affyRaw))

#################################################################################################################
##########################			Quality Control of The raw data	##############################
#################################################################################################################
print("Quality Control :")

exp_raw <- log2(Biobase::exprs(affyRaw))
#dim(exp_raw)
#dim(exp_raw[complete.cases(exp_raw ), ])
#dim(t(exp_raw))
#dim(na.omit(t(exp_raw)))





#PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

#percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
#sd_ratio <- sqrt(percentVar[2] / percentVar[1])

#dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                  #  GROUP = pData(affyRaw)$group,
               #     PATIENT = pData(affyRaw)$patient)

#png("PCA_rawData.png")
#ggplot(dataGG, aes(PC1, PC2)) +
 #     geom_point(aes(shape = GROUP, colour = PATIENT)) +
 # ggtitle("PCA plot of the log2-transformed raw expression data") +
 # xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
#  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
#  theme(plot.title = element_text(hjust = 0.5))+
#  coord_fixed(ratio = sd_ratio)
 
#dev.off()


png("Boxplot_Intensities_rawData.png",width=1000,height=1000)
oligo::boxplot(affyRaw, target = "core", main = "Boxplot of log2-intensitites for the raw data")
dev.off()


#################################################################################################################
##########################			arrayQualityMetrics								#############################
#################################################################################################################
print("ArrayQualityMetrics :")

#Until now, we have only performed a very basic quality control; more elaborate quality control plots are available in the #package arrayQualityMetrics (5). The package produces an html report, containing the quality control plots together with #a description of their aims and an identification of possible outliers. We do not discuss this tool in detail here, but #simply provide the code below, which creates a report for our raw data.


#arrayQualityMetrics(expressionset = affyRaw, outdir = "./report/",force = TRUE, do.logtransform = TRUE,intgroup = c("group", "patient"))

#################################################################################################################
##########################			Relative Log Expression data quality analysis	##############################
#################################################################################################################

print("Relative Log Expression data quality analysis")

#Before calibrating and evaluating the data, we want to perform another quality control procedure, namely Relative Log #Expression (RLE), as described in the article by Gandolfo et al (9). To this end, we first perform an RMA without prior #normalization
#The RLE is performed by calculating the median log2 intensity of every transcript across all arrays.
#Note that we do not have to apply the log2 manually, as the output data of the RMA function is in log2 scale by default.
#We then substract this transcript median intensity from every transcript intensity via the sweep function.


affyRaw_eset_notNormalize <- oligo::rma(affyRaw, target = "core", normalize = FALSE)
#summarised to Exons
row_medians_assayData <- Biobase::rowMedians(as.matrix(Biobase::exprs(affyRaw_eset_notNormalize)))

RLE_data <- sweep(Biobase::exprs(affyRaw_eset_notNormalize), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <-  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

png("RLE.png",width=1000,height=1000 )
ggplot2::ggplot(RLE_data_gathered, aes(patient_array, log2_expression_deviation)) + geom_boxplot(outlier.shape = NA) + ylim(c(-2, 2)) + theme(axis.text.x = element_text(colour = "aquamarine4",  angle = 60, size = 6.5, hjust = 1 ,face = "bold"))
dev.off()

#################################################################################################################
##########################			Quality assessment of the calibrated data PCA  	##############################
#################################################################################################################

print("Normalisation by RMA & Quality assessment of the calibrated data PCA")

affyRaw_eset_norm <- oligo::rma(affyRaw, target = "core")


png("Boxplot_Intensities_aterRMA.png",width=1000,height=1000)
oligo::boxplot(affyRaw_eset_norm, target = "core", main = "Boxplot of log2-intensitites for the normalised data")
dev.off()

exp <- Biobase::exprs(affyRaw_eset_norm)

PCA <- prcomp(t(exp), scale. = TRUE)
#PCA <- prcomp(~ ., data=t(exp),na.action=na.omit, scale. = TRUE)
#
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],GROUP = Biobase::pData(affyRaw_eset_norm)$group,PATIENT = Biobase::pData(affyRaw_eset_norm)$patient)

png("PCA_afterRMA.png")
ggplot(dataGG, aes(PC1, PC2)) +
      geom_point(aes(shape = GROUP, colour = PATIENT,size=4)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) 
dev.off()


######################################################################
####### FilterOut Low Expressed ####################
#######################################################################

print(" FilterOut Low Expressed : ")

man_threshold <- 4

png("HistogramMedianNormIntensities.png",width=1000,height=1000)
medians <- rowMedians(Biobase::exprs(affyRaw_eset_norm))
hist(medians, 100, col = "cornsilk", freq = FALSE, 
            main = "Histogram of the median intensities",
            border = "antiquewhite4",
            xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)
dev.off()

#Transcripts that do not have intensities larger than the threshold in at least as many arrays as the smallest #experimental group are excluded.

no_of_samples <- table(paste0(pData(affyRaw_eset_norm)$group, "_", pData(affyRaw_eset_norm)$patient))
print(no_of_samples)
samples_cutoff <- min(no_of_samples)

idx_man_threshold <- apply(Biobase::exprs(affyRaw_eset_norm), 1,
                           function(x){
                          sum(x > man_threshold) >= samples_cutoff})
                          table(idx_man_threshold)

print(idx_man_threshold)
affyRaw_eset_norm_man_filtered <- subset(affyRaw_eset_norm, idx_man_threshold)

#affyRaw_eset_norm_man_filtered <- affyRaw_eset_norm
######################################################################
####### Annotation of the transcript clusters ##########################################
#######################################################################
print(" Annotation of the transcript clusters : ")

anno_subset <- AnnotationDbi::select(hugene21sttranscriptcluster.db,
                                  keys = (featureNames(affyRaw_eset_norm_man_filtered)),
                                  columns = c("CHR","CHRLOC","CHRLOCEND","SYMBOL", "GENENAME","ENSEMBL"),
                                  keytype = "PROBEID")
#      PROBEID CHR     CHRLOC CHRLOCCHR  CHRLOCEND       SYMBOL
#2245 16657436   1      11873         1      14409      DDX11L1
#2246 16657436   1     182387         1     184878 LOC102725121
#Chromosomal locations on both the sense and antisense strands are measured as the number of base
#pairs from the p (5’ end of the sense strand) to q (3’ end of the sense strand) arms. Chromosomal
#locations on the antisense strand have a leading "-" sign (e. g. -1234567).
#Since some genes have multiple start sites, this field can map to multiple locat
print(subset(anno_subset, SYMBOL=="ZDHHC14"))#CST3

print("anno_subset")
dim(anno_subset)
anno_subset <- subset(anno_subset, !is.na(SYMBOL))
dim(anno_subset)
###############################################################################################
#### Removing probeset mappings to several Genes ####################################################
###############################################################################################
print(" Removing multiple mappings : ")

anno_grouped    <- group_by(anno_subset, PROBEID)
print(subset(anno_grouped, PROBEID=="16917939"))#CST3
print(subset(anno_grouped, PROBEID=="17008341"))#DAAM2
print(subset(anno_grouped, PROBEID=="17014064"))#ZDHHC14 # probeset annotated on 2 places

anno_summarized <-  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
anno_filtered <- filter(anno_summarized, no_of_matches > 1)

probe_stats <- anno_filtered 
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMax(probe_stats)#35
head(probe_stats)
print(subset(probe_stats, PROBEID=="16917939"))#CST3
print(subset(probe_stats, PROBEID=="17008341"))#DAAM2
print(subset(probe_stats, PROBEID=="17014064"))#ZDHHC14 # probeset annotated on 2 places


nrow(probe_stats)


ids_to_exlude <- (featureNames(affyRaw_eset_norm_man_filtered) %in% probe_stats$PROBEID)
print("ids_to_exlude: ")
table(ids_to_exlude)
print("before / after : ")
dim(affyRaw_eset_norm_man_filtered)
affy_final <- subset(affyRaw_eset_norm_man_filtered, !ids_to_exlude)
#affy_final <- affyRaw_eset_norm_man_filtered

dim(affy_final)

validObject(affy_final)
#As we have just excluded probe IDs from the assay data, we now have to also exclude them from the feature data anno_palmieri:

# I remove every annot that can be annotated on different stuff
print("Remove annot on differrent place : ")

head("anno_subset")
dim(anno_subset)
ids_to_exlude_for_anno_subset <- subset(anno_subset, anno_subset$PROBEID %in% probe_stats$PROBEID)
head("id to excludefor_annot")
head(ids_to_exlude_for_anno_subset,2)

# WTF
#anno_subset <-  subset(anno_subset, !ids_to_exlude_for_anno_subset)
anno_subset <- anno_subset[!duplicated(anno_subset$PROBEID),]

print("anno_subset")
dim(anno_subset)
head(anno_subset,2)

fData(affy_final)$PROBEID   <- rownames(fData(affy_final))
fData(affy_final)           <- left_join(fData(affy_final), anno_subset)#,by = "PROBEID"
# https://support.bioconductor.org/p/90481/ still there is NA Values
dim(affy_final)

rownames(fData(affy_final)) <- fData(affy_final)$PROBEID 
    
validObject(affy_final)


print("Remove NA via CHRLOC")
index.where.na <- is.na(fData(affy_final)$CHRLOC )
head(index.where.na)

# REMOVE NA CHRLOC
affy_final<- affy_final[!index.where.na,]


dim(final_annotated_table_for_probeset)
final_annotated_table_for_probeset<- filter(final_annotated_table_for_probeset, type == "main")
chromsClean<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
final_annotated_table_for_probeset<- final_annotated_table_for_probeset[final_annotated_table_for_probeset$chrom %in% chromsClean ,]
dim(final_annotated_table_for_probeset)

print("Remove probeID not Main")
index.where.probe.main <-  (fData(affy_final)$PROBEID %in% final_annotated_table_for_probeset$fsetid)
affy_final<- affy_final[index.where.probe.main,]

#######################################################################
#### Differential Analysis     ########################################
#######################################################################

print("Differential Analysis :")

#However, since we have two arrays per individual patient, we have a “Paired Samples” design (see section 9.4 of the limma #user guide). #This means that the samples might be biased by the person they come from. Whenever a feature in an #experimental setup is expected to #have a systematic influence on the result, blocking factors on these features should be #introduced.

#Thus, the first factor we need is a blocking factor for the individuals that will absorb differences in expression #between them. 
#Therefore, we block on patients, which means that the patient IDs become variables of the linear model.

#Then we create factors that give us the grouping for the tissue types (non-inflamed and inflamed).
print(Biobase::pData(affy_final))
#https://support.bioconductor.org/p/68307/
# https://bioconductor.org/packages/3.10/bioc/vignettes/limma/inst/doc/usersguide.pdf 9.4.2
patients <-factor(as.character(Biobase::pData(affy_final)$patient))
#files<- as.character(Biobase::pData(affy_final)$filenames)
groups <- factor(as.character(Biobase::pData(affy_final)$group),levels=c("FOYER","TUMOR"))


#design_affy <- model.matrix(~ groups + patients)
design_affy <- model.matrix(~  patients + groups )
print("design_affy")
head(design_affy)

#colnames(design_affy)[1:2] <- c("FOYER", "TUMOR")
rownames(design_affy) <-  patients


#contrast_matrix <- makeContrasts(FOYER-TUMOR, levels = design_affy)
#print("contrast")
#print(contrast_matrix)
fit<-lmFit(affy_final,design = design_affy)
#affy_fit <- eBayes(contrasts.fit(fit,contrast_matrix))
affy_fit <- eBayes(fit)#,trend=TRUE

print("design_affy")
head(design_affy)


table <- topTable(affy_fit, number = Inf, coef="groupsTUMOR")

png("HistogramPvalue.png")
hist(table$P.Value, col = brewer.pal(3, name = "Set2")[1],main = "Foyer vs Tumor", xlab = "p-values")
dev.off()


table <- topTable(affy_fit, number = Inf, coef="groupsTUMOR",lfc=log2(1.3))
table$FC <- logratio2foldchange(table$logFC, base=2)

head(table)

nrow(subset(table, P.Value < 0.05))
nrow(subset(table))
toPrint <- subset(table,P.Value < 0.05)# P.Value <0.001
# Several Probset we keep the main 

write.csv(toPrint,file="GeneDiffExpressed.csv",row.names=FALSE)

print("volcano")
#print(affy_fit)
volcano_names <- ifelse(abs(affy_fit$coefficients)>=1,affy_fit$genes$SYMBOL, NA)
#print(volcano_names)
png("VolcanoPlotDifftop100.png",width=1000,height=1000)
           
volcanoplot(affy_fit, coef = "groupsTUMOR", style = "p-value", highlight = 100, 
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

dev.off() 
