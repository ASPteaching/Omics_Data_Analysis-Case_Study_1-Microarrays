## ----setup, include=FALSE------------------------------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, 
                      comment = NA, prompt = TRUE, tidy = FALSE, 
                      fig.width = 7, fig.height = 7, fig_caption = TRUE,
                      cache=FALSE)
Sys.setlocale("LC_TIME", "C")


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------
if(!(require(printr))) {
  install.packages(
    'printr',
    type = 'source',
    repos = c('http://yihui.name/xran', 'http://cran.rstudio.com')
  )
}


## ----dataset, fig.cap="A simplified view of a gene expression matrix", echo=FALSE----------------------------------------------
knitr::include_graphics("figures/Figure1.jpg")


## ----MDAProcess, fig.cap="The microarray data analysis process", echo=FALSE----------------------------------------------------
knitr::include_graphics("figures/Figure2.png")


## ----CreateFolders, warning=FALSE, eval=FALSE----------------------------------------------------------------------------------
## setwd(".")
## dir.create("data")
## dir.create("results")


## ----ReadTargets---------------------------------------------------------------------------------------------------------------
targets <- read.csv2("./data/targets.csv", header = TRUE, sep = ";") 
knitr::kable(
  targets, booktabs = TRUE,
  caption = 'Content of the targets file used for the current analysis')


## ----installBioC, message=FALSE, warning=FALSE, eval=FALSE---------------------------------------------------------------------
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## BiocManager::install()


## ----installPackages, message=FALSE, warning=FALSE, eval=FALSE-----------------------------------------------------------------
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
## BiocManager::install("limma")
## BiocManager::install("genefilter")
## BiocManager::install("mogene21sttranscriptcluster.db")
## BiocManager::install("annotate")
## BiocManager::install("org.Mm.eg.db")
## BiocManager::install("ReactomePA")
## BiocManager::install("reactome.db")


## ----ReadCELfiles, message=FALSE, results='hide', warning=FALSE----------------------------------------------------------------
library(oligo)
celFiles <- list.celfiles("./data", full.names = TRUE)
library(Biobase)
my.targets <-read.AnnotatedDataFrame(file.path("./data","targets.csv"), 
                                     header = TRUE, row.names = 1, 
                                     sep=";") 
rawData <- read.celfiles(celFiles, phenoData = my.targets)


## ----ChangeName----------------------------------------------------------------------------------------------------------------
my.targets@data$ShortName->rownames(pData(rawData))
colnames(rawData) <-rownames(pData(rawData)) 

head(rawData)


## ----QCRaw, message=FALSE, warning=FALSE, eval=FALSE---------------------------------------------------------------------------
## library(arrayQualityMetrics)
## arrayQualityMetrics(rawData)


## ----QCRawDataRes, fig.cap="Aspect of the summary table, in the index.html file, produced by the arrayQualityMetrics package on the raw data", echo=FALSE----
knitr::include_graphics("figures/Figure3.png")


## ------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
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


## ----PCARaw, message=FALSE, fig.cap="Visualization of the two first Principal Components for raw data"-------------------------
plotPCA3(exprs(rawData), labels = targets$ShortName, factor = targets$Group, 
         title="Raw data", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow"))


## ----BoxplotRaw, message=FALSE, fig.cap="Boxplot for arrays intensities (Raw Data)"--------------------------------------------
boxplot(rawData, cex.axis=0.5, las=2,  which="all", 
         col = c(rep("red", 3), rep("blue", 3), rep("green", 3), rep("yellow", 3)),
         main="Distribution of raw intensity values")


## ----Normalization-------------------------------------------------------------------------------------------------------------
eset_rma <- rma(rawData)


## ----QCNorm, message=FALSE, warning=FALSE, eval=FALSE--------------------------------------------------------------------------
## arrayQualityMetrics(eset_rma, outdir = file.path("./results", "QCDir.Norm"), force=TRUE)


## ----QCNormDataRes, fig.cap="Aspect of the summary table, in the index.html file, produced by the arrayQualityMetrics package on normalized data", echo=FALSE----
knitr::include_graphics("figures/Figure6.png")


## ----fig:PCANorm, message=FALSE, fig.cap="Visualization of first two principal components for normalized data"-----------------
plotPCA3(exprs(eset_rma), labels = targets$ShortName, factor = targets$Group, 
         title="Normalized data", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow"))


## ----BoxplotNorm, message=FALSE, fig.cap="Distribution of  intensities for normalized data"------------------------------------
boxplot(eset_rma, cex.axis=0.5, las=2,  which="all", 
         col = c(rep("red", 3), rep("blue", 3), rep("green", 3), rep("yellow", 3)),
         main="Boxplot for arrays intensity: Normalized Data")


## ----BatchDetection, message=FALSE, warning=FALSE------------------------------------------------------------------------------
#load the library
library(pvca)
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
axis(1, at = bp, labels = pvcaObj$label, cex.axis = 0.75, las=2)
values = pvcaObj$dat
new_values = round(values , 3)
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.7)


## ----SDplot, fig.cap="Values of standard deviations allong all samples for all genes ordered from smallest to biggest"---------
sds <- apply (exprs(eset_rma), 1, sd)
sdsO<- sort(sds)
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))


## ----Filtering1, results='hide', message=FALSE---------------------------------------------------------------------------------
library(genefilter)
library(mogene21sttranscriptcluster.db)
annotation(eset_rma) <- "mogene21sttranscriptcluster.db"
filtered <- nsFilter(eset_rma, 
                     require.entrez = TRUE, remove.dupEntrez = TRUE,
                     var.filter=TRUE, var.func=IQR, var.cutoff=0.75, 
                     filterByQuantile=TRUE, feature.exclude = "^AFFX")


## ----FilterResults1, results='hide', echo=FALSE--------------------------------------------------------------------------------
names(filtered)
class(filtered$eset)


## ----FilterResults2------------------------------------------------------------------------------------------------------------
print(filtered$filter.log)
eset_filtered <-filtered$eset


## ----SaveData1, results='hide', message=FALSE----------------------------------------------------------------------------------
write.csv(exprs(eset_rma), file="./results/normalized.Data.csv")
write.csv(exprs(eset_filtered), file="./results/normalized.Filtered.Data.csv")
save(eset_rma, eset_filtered, file="./results/normalized.Data.Rda")


## ----LoadSavedData-------------------------------------------------------------------------------------------------------------
if (!exists("eset_filtered")) load (file="./results/normalized.Data.Rda")


## ----DesignMatrix, message=FALSE-----------------------------------------------------------------------------------------------
library(limma)
designMat<- model.matrix(~0+Group, pData(eset_filtered))
colnames(designMat) <- c("KO.COLD", "KO.RT", "WT.COLD", "WT.RT")
print(designMat)


## ----setContrasts--------------------------------------------------------------------------------------------------------------
cont.matrix <- makeContrasts (KOvsWT.COLD = KO.COLD-WT.COLD,
                              KOvsWT.RT = KO.RT-WT.RT,
                              INT = (KO.COLD-WT.COLD) - (KO.RT-WT.RT),
                              levels=designMat)
print(cont.matrix)


## ---- linearmodelfit-----------------------------------------------------------------------------------------------------------
library(limma)
fit<-lmFit(eset_filtered, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)


## ---- topTabs1-----------------------------------------------------------------------------------------------------------------
topTab_KOvsWT.COLD <- topTable (fit.main, number=nrow(fit.main), coef="KOvsWT.COLD", adjust="fdr") 
head(topTab_KOvsWT.COLD)


## ---- topTabs2-----------------------------------------------------------------------------------------------------------------
topTab_KOvsWT.RT <- topTable (fit.main, number=nrow(fit.main), coef="KOvsWT.RT", adjust="fdr") 
head(topTab_KOvsWT.RT)


## ---- topTabs3-----------------------------------------------------------------------------------------------------------------
topTab_INT  <- topTable (fit.main, number=nrow(fit.main), coef="INT", adjust="fdr") 
head(topTab_INT)


## ----GeneAnnotation, message=FALSE, warning=FALSE------------------------------------------------------------------------------
annotatedTopTable <- function(topTab, anotPackage)
{
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  thePackage <- eval(parse(text = anotPackage))
  geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
  annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
return(annotatedTopTab)
}


## ----annotateTopTables---------------------------------------------------------------------------------------------------------
topAnnotated_KOvsWT.COLD <- annotatedTopTable(topTab_KOvsWT.COLD,
anotPackage="mogene21sttranscriptcluster.db")
topAnnotated_KOvsWT.RT <- annotatedTopTable(topTab_KOvsWT.RT,
anotPackage="mogene21sttranscriptcluster.db")
topAnnotated_INT <- annotatedTopTable(topTab_INT,
anotPackage="mogene21sttranscriptcluster.db")
write.csv(topAnnotated_KOvsWT.COLD, file="./results/topAnnotated_KOvsWT_COLD.csv")
write.csv(topAnnotated_KOvsWT.RT, file="./results/topAnnotated_KOvsWT_RT.csv")
write.csv(topAnnotated_INT, file="./results/topAnnotated_INT.csv")


## ----annotatedTop, echo=FALSE--------------------------------------------------------------------------------------------------
short<- head(topAnnotated_KOvsWT.COLD[1:5,1:4])
# library(kableExtra)
# knitr::kable(
#   short, booktabs = TRUE,
#   caption = 'Annotations added to results "topTable" for the comparison "KOvsWT.COLD"'
# )
show(short)


## ----volcanoPlot, fig.cap="Volcano plot for the comparison between KO and WT in COLD temperature. The names of the top 4 genes (i.e. the first four genes in the topTable) are shown in the plot"----
library(mogene21sttranscriptcluster.db)
geneSymbols <- select(mogene21sttranscriptcluster.db, rownames(fit.main), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
volcanoplot(fit.main, coef=1, highlight=4, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
  abline(v=c(-1,1))


## ----saveVolcanos, echo=FALSE, results='hide'----------------------------------------------------------------------------------
pdf("figures/Volcanos.pdf")
for (i in colnames(cont.matrix)){
  volcanoplot(fit.main, coef=i, highlight=4, names=SYMBOLS,
              main=paste("Differentially expressed genes",i, sep="\n"))
  abline(v=c(-1,1))
}
dev.off()


## ----decideTests.1-------------------------------------------------------------------------------------------------------------
library(limma)
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.1, lfc=1)


## ----resumeDecideTests---------------------------------------------------------------------------------------------------------
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))


## ---- vennDiagram, fig.cap="Venn diagram showing the genes in common between the three comparisons performed"------------------
vennDiagram (res.selected[,1:3], cex=0.9)
title("Genes in common between the three comparisons\n Genes selected with FDR < 0.1 and logFC > 1")


## ----data4Heatmap--------------------------------------------------------------------------------------------------------------
probesInHeatmap <- rownames(res.selected)
HMdata <- exprs(eset_filtered)[rownames(exprs(eset_filtered)) %in% probesInHeatmap,]

geneSymbols <- select(mogene21sttranscriptcluster.db, rownames(HMdata), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
rownames(HMdata) <- SYMBOLS
write.csv(HMdata, file = file.path("./results/data4Heatmap.csv"))


## ----heatmapNoclustering, fig.cap="Heatmap for expression data without any grouping"-------------------------------------------
my_palette <- colorRampPalette(c("blue", "red"))(n = 299)
library(gplots)

heatmap.2(HMdata,
          Rowv = FALSE,
          Colv = FALSE,
          main = "Differentially expressed genes \n FDR < 0,1, logFC >=1",
          scale = "row",
          col = my_palette,
          sepcolor = "white",
          sepwidth = c(0.05,0.05),
          cexRow = 0.5,
          cexCol = 0.9,
          key = TRUE,
          keysize = 1.5,
          density.info = "histogram",
          ColSideColors = c(rep("red",3),rep("blue",3), rep("green",3), rep("yellow",3)),
          tracecol = NULL,
          dendrogram = "none",
          srtCol = 30)


## ----heatmapClustering, fig.cap="Heatmap for expression data grouping genes (rows) and samples (columns) by their similarity"----
heatmap.2(HMdata,
          Rowv = TRUE,
          Colv = TRUE,
          dendrogram = "both",
          main = "Differentially expressed genes \n FDR < 0,1, logFC >=1",
          scale = "row",
          col = my_palette,
          sepcolor = "white",
          sepwidth = c(0.05,0.05),
          cexRow = 0.5,
          cexCol = 0.9,
          key = TRUE,
          keysize = 1.5,
          density.info = "histogram",
          ColSideColors = c(rep("red",3),rep("blue",3), rep("green",3), rep("yellow",3)),
          tracecol = NULL,
          srtCol = 30)



## ----selectGenes---------------------------------------------------------------------------------------------------------------
listOfTables <- list(KOvsWT.COLD = topTab_KOvsWT.COLD, 
                     KOvsWT.RT  = topTab_KOvsWT.RT, 
                     INT = topTab_INT)
listOfSelected <- list()
for (i in 1:length(listOfTables)){
  # select the toptable
  topTab <- listOfTables[[i]]
  # select the genes to be included in the analysis
  whichGenes<-topTab["adj.P.Val"]<0.15
  selectedIDs <- rownames(topTab)[whichGenes]
  # convert the ID to Entrez
  EntrezIDs<- select(mogene21sttranscriptcluster.db, selectedIDs, c("ENTREZID"))
  EntrezIDs <- EntrezIDs$ENTREZID
  listOfSelected[[i]] <- EntrezIDs
  names(listOfSelected)[i] <- names(listOfTables)[i]
}
sapply(listOfSelected, length)


## ------------------------------------------------------------------------------------------------------------------------------
mapped_genes2GO <- mappedkeys(org.Mm.egGO)
mapped_genes2KEGG <- mappedkeys(org.Mm.egPATH)
mapped_genes <- union(mapped_genes2GO , mapped_genes2KEGG)


## ----BiologicalSig-------------------------------------------------------------------------------------------------------------
library(ReactomePA)

listOfData <- listOfSelected[1:2]
comparisonsNames <- names(listOfData)
universe <- mapped_genes

for (i in 1:length(listOfData)){
  genesIn <- listOfData[[i]]
  comparison <- comparisonsNames[i]
  enrich.result <- enrichPathway(gene = genesIn,
                                 pvalueCutoff = 0.05,
                                 readable = T,
                                 pAdjustMethod = "BH",
                                 organism = "mouse",
                                 universe = universe)
  
  cat("##################################")
  cat("\nComparison: ", comparison,"\n")
  print(head(enrich.result))

  if (length(rownames(enrich.result@result)) != 0) {
  write.csv(as.data.frame(enrich.result), 
             file =paste0("./results/","ReactomePA.Results.",comparison,".csv"), 
             row.names = FALSE)
  
  pdf(file=paste0("./results/","ReactomePABarplot.",comparison,".pdf"))
    print(barplot(enrich.result, showCategory = 15, font.size = 4, 
            title = paste0("Reactome Pathway Analysis for ", comparison,". Barplot")))
  dev.off()
  
  pdf(file = paste0("./results/","ReactomePAcnetplot.",comparison,".pdf"))
    print(cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
         vertex.label.cex = 0.75))
  dev.off()
  }
}


## ----network, fig.cap="Network obtained from the Reactome enrichment analysis on the list obtained from the comparison between KO and WT in RT"----
  cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
         vertex.label.cex = 0.75)


## ----tableReacto, echo=FALSE---------------------------------------------------------------------------------------------------
Tab.react <- read.csv2(file.path("./results/ReactomePA.Results.KOvsWT.RT.csv"), 
                       sep = ",", header = TRUE, row.names = 1)

Tab.react <- Tab.react[1:4, 1:5]
knitr::kable(Tab.react, booktabs = TRUE, caption = "First rows and columns for Reactome results on KOvsWT.RT.csv comparison")


## ----listOfFiles, echo=FALSE---------------------------------------------------------------------------------------------------
listOfFiles <- dir("./results/") 
knitr::kable(
  listOfFiles, booktabs = TRUE,
  caption = 'List of files generated in the analysis',
  col.names="List_of_Files"
)

