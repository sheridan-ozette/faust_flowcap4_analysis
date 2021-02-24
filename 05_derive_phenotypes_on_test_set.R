library(survminer)
library(dplyr)
#
#bring in meta data
#
metaData <- readRDS(file.path(normalizePath("."),"flowReMi_metaData.rds"))
metaData$subject <- as.character(paste0("subject_",as.numeric(as.factor(metaData$subjectID))))
#
#get the markers in the selected phenotypes
#
acp <- colnames(readRDS(file.path(normalizePath("."),"part03","faustData","faustCountMatrix.rds")))
cnMap <- readRDS(file.path(normalizePath("."),"part03","faustData","metaData","exhaustiveColNameMap.rds"))
referencePhenotype <- acp[1]
encPheno <- cnMap[which(cnMap[,"newColNames"]==referencePhenotype),"faustColNames",drop=TRUE]
markersInPheno <- strsplit(encPheno,"~[[:digit:]]~[[:digit:]]~")[[1]]
#
#bring in the new annotation thresholds and derive counts for target phenotypes on test set
#
rl <- readRDS(file.path(normalizePath("."),"part04","faustData","gateData","root_resList.rds"))
allTargets <- readRDS(file.path(normalizePath("."),"selected_phenotypes.rds"))
firstTarget <- TRUE
for (targetPhenotype in allTargets) {
    encPheno <- cnMap[which(cnMap[,"newColNames"]==targetPhenotype),"faustColNames",drop=TRUE]
    firstSample <- TRUE
    for (fn in list.files(file.path(normalizePath("."),"part04","faustData","sampleData"))) {
        exprsMat <- readRDS(file.path(normalizePath("."),"part04","faustData","sampleData",fn,"exprsMat.rds"))
        subject <- metaData[which(metaData$name==fn),"subject"]
        countLookup <- seq(nrow(exprsMat))
        for (mn in markersInPheno) {
            if (grepl(paste0(mn,"~1"),encPheno)) {
                countLookup <- intersect(countLookup,
                                         which(exprsMat[,mn] < as.numeric(rl[[mn]][[subject]])))
            }
            else {
                countLookup <- intersect(countLookup,
                                         which(exprsMat[,mn] >= as.numeric(rl[[mn]][[subject]])))
            }
        }
        parentCount <- nrow(exprsMat)
        phenotypeDF <- data.frame(
            sampleName=fn,
            prop=(length(countLookup)/parentCount),
            stringsAsFactors=FALSE
        )
        colnames(phenotypeDF) <- c("name",targetPhenotype)
        if (firstSample) {
            allPhenotypeDF <- phenotypeDF
            firstSample <- FALSE
        }
        else {
            allPhenotypeDF <- rbind(allPhenotypeDF,phenotypeDF)
        }
    }
    if (firstTarget) {
        phenotypesOut <- allPhenotypeDF
        firstTarget <- FALSE
    }
    else {
        phenotypesOut <- inner_join(phenotypesOut,allPhenotypeDF,by=c("name"))
    }
}
#
#save the count matrix
#
saveRDS(phenotypesOut,file.path(normalizePath("."),"derived_phenotype_matrix.rds"))

