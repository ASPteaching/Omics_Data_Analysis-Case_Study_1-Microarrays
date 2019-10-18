## ----setup, include=FALSE------------------------------------------------
require(knitr)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, 
                      comment = NA, prompt = TRUE, tidy = FALSE, 
                      fig.width = 7, fig.height = 7, fig_caption = TRUE,
                      cache=FALSE)
Sys.setlocale("LC_TIME", "C")


## ----echo=FALSE----------------------------------------------------------
if(!(require(printr))) {
  install.packages(
    'printr',
    type = 'source',
    repos = c('http://yihui.name/xran', 'http://cran.rstudio.com')
  )
}


## ----dataset, fig.cap="A simplified view of a gene expression matrix", echo=FALSE----
knitr::include_graphics("figures/Figure1.jpg")


## ----MDAProcess, fig.cap="The microarray data analysis process", echo=FALSE----
knitr::include_graphics("figures/Figure2.png")


## ----CreateFolders, warning=FALSE, eval=FALSE----------------------------
## setwd(".")
## dir.create("data")
## dir.create("results")


## ----ReadTargets---------------------------------------------------------
targets <- read.csv2("./data/targets.csv", header = TRUE, sep = ";") 
knitr::kable(
  targets, booktabs = TRUE,
  caption = 'Content of the targets file used for the current analysis')


## ----installBioC, message=FALSE, warning=FALSE, eval=FALSE---------------
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## BiocManager::install()


## ----installPackages, message=FALSE, warning=FALSE, eval=FALSE-----------
## install.packages("knitr")
## install.packages("colorspace")
## install.packages("gplots")
## install.packages("ggplot2")
## install.packages("ggrepel")
## install.packages("htmlTable")
## install.packages("prettydoc")
## install.packages("devtools")
## install.packages("BiocManager")
## BiocManager::install("oligo")
## BiocManager::install("pd.mogene.2.1.st")
## BiocManager::install("arrayQualityMetrics")
## BiocManager::install("pvca")
## # NOT NEEDED UNTIL ANALYSES ARE PERFORMED
## # BiocManager::install("limma")
## # BiocManager::install("genefilter")
## # BiocManager::install("mogene21sttranscriptcluster.db")
## # BiocManager::install("annotate")
## # BiocManager::install("org.Mm.eg.db")
## # BiocManager::install("ReactomePA")


## ----ReadCELfiles, message=FALSE, results='hide', warning=FALSE----------
require(oligo)
celFiles <- list.celfiles("./data", full.names = TRUE)
require(Biobase)
my.targets <-read.AnnotatedDataFrame(file.path("./data","targets.csv"), 
                                     header = TRUE, row.names = 1, 
                                     sep=";") 
rawData <- read.celfiles(celFiles, phenoData = my.targets)


## ----ChangeName----------------------------------------------------------
my.targets@data$ShortName->rownames(pData(rawData))
colnames(rawData) <-rownames(pData(rawData)) 

head(rawData)


## ----QCRaw, message=FALSE, warning=FALSE, eval=FALSE---------------------
## library(arrayQualityMetrics)
## arrayQualityMetrics(rawData)


## ----QCRawDataRes, fig.cap="Aspect of the summary table, in the index.html file, produced by the arrayQualityMetrics package on the raw data", echo=FALSE----
knitr::include_graphics("figures/Figure3.png")


## ------------------------------------------------------------------------
require(ggplot2)
require(ggrepel)
plotPCA3 <- function (datos, labels, factor, title, scale,colores, size = 1.5, glineas = 0.25) {
  data <- prcomp(t(datos),scale=scale)
  # plot adjustments
  dataDf <- data.frame(data$x)
  Group <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
  # main plot
  p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Group), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
    scale_fill_discrete(name = "Group")
  # avoiding labels superposition
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels),segment.size = 0.25, size = size) + 
    labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +  
    ggtitle(paste("Principal Component Analysis for: ",title,sep=" "))+ 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colores)
  }


## ----PCARaw, message=FALSE, fig.cap="Visualization of the two first Principal Components for raw data"----
plotPCA3(exprs(rawData), labels = targets$ShortName, factor = targets$Group, 
         title="Raw data", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow"))


## ----savePCAraw, echo=TRUE, results='hide'-------------------------------
tiff("figures/PCA_RawData.tiff", res = 200, width = 4.5, height = 4, units = 'in')
plotPCA3(exprs(rawData), labels = targets$ShortName, factor = targets$Group, 
         title="Raw data", scale = FALSE, size = 2, 
         colores = c("red", "blue", "green", "yellow"))
dev.off()


## ----BoxplotRaw, message=FALSE, fig.cap="Boxplot for arrays intensities (Raw Data)"----
boxplot(rawData, cex.axis=0.5, las=2,  which="all", 
         col = c(rep("red", 3), rep("blue", 3), rep("green", 3), rep("yellow", 3)),
         main="Distribution of raw intensity values")


## ----saveIntensRaw, echo=FALSE, results='hide'---------------------------
tiff("figures/Intensity_RawData.tiff", res = 200, width = 4, height = 4, units = 'in')
boxplot(rawData, cex.axis=0.5, las=2,  which="all", 
         col = c(rep("red", 3), rep("blue", 3), rep("green", 3), rep("yellow", 3)),
         main="Distribution of raw intensity values")
dev.off()


## ----Normalization-------------------------------------------------------
eset_rma <- rma(rawData)


## ----QCNorm, message=FALSE, warning=FALSE, eval=FALSE--------------------
## arrayQualityMetrics(eset_rma, outdir = file.path("./results", "QCDir.Norm"), force=TRUE)


## ----QCNormDataRes, fig.cap="Aspect of the summary table, in the index.html file, produced by the arrayQualityMetrics package on normalized data", echo=FALSE----
knitr::include_graphics("figures/Figure6.png")


## ----PCANorm, message=FALSE, fig.cap="Visualization of first two principal components for normalized data"----
plotPCA3(exprs(eset_rma), labels = targets$ShortName, factor = targets$Group, 
         title="Normalized data", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow"))


## ----savePCAnorm, echo=FALSE, results='hide'-----------------------------
tiff("figures/PCA_NormData.tiff", res = 150, width = 5, height = 5, units = 'in')
plotPCA3(exprs(eset_rma), labels = targets$ShortName, factor = targets$Group, 
         title="Normalized data", scale = FALSE, size = 2, 
         colores = c("red", "blue", "green", "yellow"))
dev.off()


## ----BoxplotNorm, message=FALSE, fig.cap="Distribution of  intensities for normalized data"----
boxplot(eset_rma, cex.axis=0.5, las=2,  which="all", 
         col = c(rep("red", 3), rep("blue", 3), rep("green", 3), rep("yellow", 3)),
         main="Boxplot for arrays intensity: Normalized Data")


## ----saveIntensNorm, echo=FALSE, results='hide'--------------------------
tiff("figures/Intensity_NormData.tiff", res = 150, width = 5, height = 5, units = 'in')
boxplot(eset_rma, cex.axis=0.5, las=2,  which="all", 
         col = c(rep("red", 3), rep("blue", 3), rep("green", 3), rep("yellow", 3)),
         main="Boxplot for arrays intensity: Normalized Data")
dev.off()


## ----BatchDetection, message=FALSE, warning=FALSE------------------------
#load the library
require(pvca)
pData(eset_rma) <- targets
#select the threshold
pct_threshold <- 0.6
#select the factors to analyze
batch.factors <- c("Genotype", "Temperature")
#run the analysis
pvcaObj <- pvcaBatchAssess (eset_rma, batch.factors, pct_threshold)


## ----plotPVCA, fig.cap="Relative importance of the different factors -genotype, temperature and interaction- affecting gene expression"----
#plot the results
bp <- barplot(pvcaObj$dat, xlab = "Effects",
  ylab = "Weighted average proportion variance",
  ylim= c(0,1.1),col = c("mediumorchid"), las=2,
  main="PVCA estimation")
axis(1, at = bp, labels = pvcaObj$label, cex.axis = 0.55, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.5)


## ----savePVCAplot, echo=FALSE, results='hide'----------------------------
tiff("figures/PVCAplot.tiff", res = 150, width = 5, height = 5, units = 'in')
bp <- barplot(pvcaObj$dat, xlab = "Effects",
  ylab = "Weighted average proportion variance",
  ylim= c(0,1.1),col = c("mediumorchid"), las=2,
  main="PVCA estimation")
axis(1, at = bp, labels = pvcaObj$label, cex.axis = 0.45, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.5)
dev.off()

