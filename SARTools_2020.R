################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### November 28th, 2019
### designed to be executed with SARTools 1.7.2
################################################################################
# library(DESeq2)         ## RNA-seq differential analysis
# library(limma)          ## RNA-seq differential analysis via limma-voom
# library(edgeR)          ## RNA-seq differential analysis
# library(RColorBrewer)   ## Colors for plotting
# library(venn)           ## Venn Diagrams
# library(HTSFilter)       ## To remove not expressed genes
# library(readxl)
# library(tidyverse) 
# # library(devtools)
# # utils::download.file(source_url, destfile="aa.zip")
# # source_url <- "https://api.github.com/repos/PF2-pasteur-fr/SARTools/tarball/master"
# # getOption("download.file.method")
# # capabilities(c("libcurl", "http/ftp"))
# # install_github("PF2-pasteur-fr/SARTools")
# library(SARTools)
################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session
setwd("C:/Users/a0148148/Documents/Dossiers Olivier HARLE/INRA/EtudeConsortia/ARN/DESeq2_B") ; dir()
workDir <- "C:/Users/a0148148/Documents/Dossiers Olivier HARLE/INRA/EtudeConsortia/ARN/DESeq2_B"      # working directory for the R session
projectName <- "IBIS"                         # name of the project
author <- "OH"                                # author of the statistical analysis/report

targetFile <- "target_IBIS_clearM5045t.txt"                           # path to the design/target file
rawDir <- "perfectcount"                                      # path to the directory containing raw counts files
featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove

varInt <- "pH"                                    # factor of interest
condRef <- "50"                                    # reference biological condition
batch <- "extract"                                       # blocking factor: NULL (default) or "batch" for example

fitType <- "parametric"                              # mean-variance relationship: "parametric" (default), "local" or "mean"
cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.05                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "shorth"    #first try with median          # "median" (default) or "shorth" to estimate the size factors

colors <- c("darkorange3",# "gold")#,"gold2","gold3",        # vector of colors of each biological condition on the plots
    "darkorange2")#,"darkorange2","darkorange1","darkorange")
 # "red3","red2")#,"red3")    #"red",sans "black" --> sans témoins et sans B60b
forceCairoGraph <- FALSE

################################################################################
###                             running script                               ###
################################################################################
# setwd(workDir)
library(SARTools)
if (forceCairoGraph) options(bitmapType="cairo")

# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)
target <-target[order(target$pH),]
target$group<-as.factor(as.character(target$group))
# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove) #, featuresToRemove=featuresToRemove
# counts <- counts[-grep("MaGe:712",row.names(as.data.frame(counts))),]#for DESeq2 on Lpla777
# counts <- counts[-grep("MaGe:715",row.names(as.data.frame(counts))),]#for DESeq2 on Ldel865
cAsup <- NULL
for(c in 1:ncol(counts)){
  if(sum(counts[,c])<(100000)){
    cAsup <- c(cAsup,c)
  }
}
# names(as.data.frame(counts))
# rAsup <- NULL
# for(r in 1:nrow(counts)){
#   if((sum(counts[r,])<(ncol(counts)*10))){ 
#     rAsup <- c(rAsup,r)
#   }
# }
counts <- counts[,-cAsup] #-rAsup
# counts <- as.data.frame(counts)
target <- target[-cAsup,]
# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# save image of the R session
save.image(file=paste0(projectName,"shorth", ".RData"))

# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)

