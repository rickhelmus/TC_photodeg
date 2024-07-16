library(patRoon)

# NOTE: use OpenMS 2.7 instead of bundled version from patRoonExt
options(patRoon.path.OpenMS = "C:/Program Files/OpenMS-2.7.0/bin")

future::plan("multisession", workers = 14)

for (f in list.files(c("workflow", "utils"), pattern = "\\.R$", full.names = TRUE))
    source(f)

wfData <- list()
wfData$featData <- getFeatures()

bgMSMS <- getBGMSMSPeaks(analysisInfo(wfData$featData$fGroupsRaw)[analysisInfo(wfData$featData$fGroupsRaw)$group == "MQ", ])

wfData$SuSData <- getSuspectsFromStructures(wfData$featData, bgMSMS)
wfData$SuFData <- getSuspectsFromFormulas(wfData$featData, wfData$SuSData, bgMSMS)
wfData$unkData <- getUnknownsGeneral(wfData$featData, wfData$SuSData, wfData$SuFData, bgMSMS)
wfData$unkFormData <- getUnknownsFormulas(wfData$unkData, wfData$SuSData)
wfData$unkCompData <- getUnknownsCompounds(wfData$unkData, wfData$SuSData)
wfData$confData <- getIDLConfirmations()
wfData$SQData <- getSemiQuant(wfData$featData, wfData$SuSData, wfData$unkData, wfData$unkFormData, wfData$unkCompData,
                              wfData$confData$suspList)
makeReport(wfData$featData, wfData$SuSData, wfData$SuFData, wfData$unkData, wfData$unkFormData, wfData$unkCompData,
           wfData$confData)
wfData$unkEvalData <- getUnkMetricsEval(wfData$SuSData)

# Save workflow data for further post-processing (plotting etc)
saveRDS(wfData, "output/workflow.Rds")
