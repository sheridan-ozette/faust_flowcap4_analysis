library(faust)
library(flowCore)
library(flowWorkspace)
library(dplyr)
library(xtable)

gsIn <- load_gs(file.path(normalizePath("."),"flowReMi_gs"))
pd <- pData(gsIn)
pd$subject <- as.character(paste0("subject_",as.numeric(as.factor(pd$subjectID))))
trainStatus <- rep("Training",nrow(pd))
trainStatus[which(pd$inTraining=="0")] <- "Test"
pd$trainStatus <- trainStatus
pData(gsIn) <- pd

#
#standardize the marker names in the gating set.
#
refTable <- as.data.frame(flowCore::parameters(gh_pop_get_data(gsIn[["010_preprocessed.fcs"]],"root"))@data)
for (i in seq(length(gsIn))) {
    print(i)
    pdTable <- as.data.frame(flowCore::parameters(gh_pop_get_data(gsIn[[i]],"root"))@data)
    mn <- pdTable[,"desc"]
    names(mn) <- pdTable[,"name"]
    mn[match(refTable[,"name"],names(mn))] <- refTable[,"desc",drop=TRUE]
    markernames(gsIn[[i]]) <- mn
}

#
#set the starting node of the pre-gating.
#
startingNode <- "root"

#
#subset to training samples.
#
trainingSamples <- rownames(pd)[which(pd$trainStatus == "Training")]
gsTrain <- gsIn[which(sampleNames(gsIn) %in% trainingSamples)]

#
#set the markers used in the analysis.
#
ac <- c("CD4", "CD27", "CCR7", "CD8", "CD57", "CD45RO", "CD107A","CD154", "IFNG", "TNFA", "IL2")

#
#set channel bounds
#
cbIN <- structure(
    c(0.5, 3.407, 0.5, 2.75, 0.5, 2.5, 0.5, 3.577, 0.5, 2.75, 0.5, 2.819, 0.5, 2.779,
      0.5, 2.531, 0.5, 2.603, 0.5, 1.5, 0.5, 2.551),
    .Dim = c(2L, 11L),
    .Dimnames = list(
        c("Low", "High"),
        c("CD4", "CD27", "CD8", "CD57", "CD45RO", "CD107A", "CD154", "CCR7", "IFNG", "TNFA", "IL2")
    )
)
xtable(t(cbIN))

if (!dir.exists(file.path(normalizePath("."),"part03"))) {
    dir.create(file.path(normalizePath("."),"part03"))
}

faust(
    gatingSet=gsTrain,
    activeChannels=ac,
    channelBounds=cbIN,
    startingCellPop=startingNode,
    experimentalUnit="subject",
    projectPath = file.path(normalizePath("."),"part03"),
    depthScoreThreshold = 0.005,
    selectionQuantile = 0.75,
    threadNum = 8,
    debugFlag = TRUE,
    seedValue = 1234,
    drawAnnotationHistograms = FALSE,
    annotationsApproved = TRUE
)

selected_phenotypes <- setdiff(colnames(readRDS(file.path(normalizePath("."),
                                                          "part03",
                                                          "faustData",
                                                          "faustCountMatrix.rds"))),"0_0_0_0_0")
saveRDS(selected_phenotypes,"./selected_phenotypes.rds")
