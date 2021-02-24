library(stringr)
library(flowCore)
library(flowWorkspace)
library(dplyr)
projPath <- normalizePath(".")
fcsFiles <- setdiff(list.files(file.path(projPath,"FR-FCM-ZZ99","dambi")),c("MetaData_preprocessing.csv"))

#
#bring in meta data and derive file for later use.
#
patMeta<- read.csv(file.path(normalizePath("."),"FR-FCM-ZZ99","attachments","MetaDataTrain.csv"),
                   as.is=TRUE)
inTraining <- rep(1,nrow(patMeta))
inTraining[which(is.na(patMeta$Status))] <- 0
trainingInfo <- data.frame(
    subjectID=paste0(patMeta$Stim,patMeta$Unstim),
    inTraining=inTraining,
    stringsAsFactors=FALSE
)
allMeta<- read.csv(file.path(normalizePath("."),"FR-FCM-ZZ99","attachments","MetaDataFull.csv"),
                   as.is=TRUE)
allMeta$subjectID <- paste0(allMeta$Stim,allMeta$Unstim)
dfStim <- allMeta[,c("Status","Survival.Time","Stim","subjectID")]
dfStim$stimStatus <- 1
dfStim$fileNum <- dfStim$Stim
dfStim <- dfStim[,c("Status","Survival.Time","subjectID","stimStatus","fileNum")]
dfUn <- allMeta[,c("Status","Survival.Time","Unstim","subjectID")]
dfUn$stimStatus <- 0
dfUn$fileNum <- dfUn$Unstim
dfUn <- dfUn[,c("Status","Survival.Time","subjectID","stimStatus","fileNum")]
fileMeta <- rbind(dfStim,dfUn)
allMeta <- inner_join(fileMeta,trainingInfo,by=c("subjectID"))
saveRDS(allMeta,
        file.path(normalizePath("."),"allSubjectsMetaData.rds"))

#
#construct the gating set for faust analysis
#
ffList <- list()
for (ff in fcsFiles) {
    print(ff)
    ffIn <- read.FCS(file.path(projPath,"FR-FCM-ZZ99","dambi",ff))
    ffList <- append(ffList,list(ffIn))
    names(ffList)[length(ffList)] <- ff
}
fs <- as(ffList,"flowSet")
gs <- GatingSet(fs)
pd <- pData(gs)
pd$fileNum <- paste0(str_extract(pd[,"name"],"\\d\\d\\d"),".fcs")
pdm <- inner_join(pd,allMeta,by=c("fileNum"))
rownames(pdm) <- pdm[,"name"]
pData(gs) <- pdm
save_gs(gs,"./flowReMi_gs")

#
#save a copy of the gating set meta data
#
gs <- load_gs("./flowReMi_gs")
flowReMi_metaData <- as.data.frame(pData(gs))
flowReMi_metaData$Survival.Time <- as.numeric(flowReMi_metaData$Survival.Time)
flowReMi_metaData$stimStatus <- as.numeric(flowReMi_metaData$stimStatus)
flowReMi_metaData$inTraining <- as.numeric(flowReMi_metaData$inTraining)
flowReMi_metaData$Status <- as.numeric(flowReMi_metaData$Status)
saveRDS(flowReMi_metaData,"./flowReMi_metaData.rds")
