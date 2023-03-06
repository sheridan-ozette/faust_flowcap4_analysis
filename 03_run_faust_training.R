library(faust)
library(faust2)
library(flowCore)
library(flowWorkspace)
library(dplyr)
library(xtable)

pop <- 'allcells'
gsIn <- load_gs(file.path(normalizePath("."),paste0("flowReMi_gs", "_", pop)))
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
    # print(i)
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

gsTrain <- gsTrain[[1:10]]

#
#set the markers used in the analysis.
#
ac <- c("CD4", "CD27", "CCR7", "CD8", "CD57", "CD45RO", "CD107A","CD154", "IFNG", "TNFA", "IL2")
# ac <- setdiff(names(markernames(gsIn)),
#               c('G780-A', 'V450-A'))

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

if (!dir.exists(file.path(normalizePath("."),paste0("part03_",pop)))) {
  dir.create(file.path(normalizePath("."),paste0("part03_",pop)))
  dir.create(file.path(normalizePath("."),paste0("part03_",pop),'faust2'))
}

faust2(
  gs_path = file.path(normalizePath("."),paste0("flowReMi_gs", "_", pop)),
  derived_data_path = file.path(normalizePath("."),
                                paste0("part03_",pop),
                                'faust2'),
  # active_channels = ac,
  active_channels = c('G710-A'),
  channel_bounds = cbIN,
  experimental_unit = "fileNum",
  starting_cell_pop=startingNode,
  depth_score_threshold = 0.005,
  selection_quantile = 0.75,
  thread_number = 8,
  random_seed = 1234,
)

faust_old(
    gatingSet=gsTrain,
    activeChannels=ac,
    channelBounds=cbIN,
    startingCellPop=startingNode,
    experimentalUnit="subject",
    projectPath = file.path(normalizePath("."),paste0("part03_",'test')),
    depthScoreThreshold = 0.005,
    selectionQuantile = 0.75,
    threadNum = 8,
    debugFlag = TRUE,
    seedValue = 1234,
    drawAnnotationHistograms = FALSE,
    annotationsApproved = TRUE
)

faust(
    gatingSet=gsTrain,
    activeChannels=ac,
    channelBounds=cbIN,
    startingCellPop=startingNode,
    experimentalUnit="subject",
    projectPath = file.path(normalizePath("."),paste0("part03_",pop)),
    depthScoreThreshold = 0.005,
    selectionQuantile = 0.75,
    threadNum = 8,
    debugFlag = TRUE,
    seedValue = 1234,
    drawAnnotationHistograms = FALSE,
    annotationsApproved = TRUE
)

selected_phenotypes <- setdiff(colnames(readRDS(file.path(normalizePath("."),
                                                          paste0("part03_",),
                                                          "faustData",
                                                          "faustCountMatrix.rds"))),"0_0_0_0_0")
saveRDS(selected_phenotypes,paste0("./selected_phenotypes_", 'test', ".rds"))
