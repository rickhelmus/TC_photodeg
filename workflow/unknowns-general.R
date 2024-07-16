getUnknownsGeneral <- function(featData, SuSData, SuFData, bgMSMS)
{
    fGroups <- screenSuspects(featData$fGroupsSel, getParentSuspList(), onlyHits = FALSE)
    
    # keep parents while removing all features from suspect screening
    fGroups <- fGroups[, union(screenInfo(fGroups)$group,
                               setdiff(names(fGroups), union(groupNames(SuSData$componentsTPs),
                                                             groupNames(SuFData$componentsTPs))))]
    
    mslists <- getMSPL(fGroups, bgMSMS)
    
    componentsTPs <- generateComponents(fGroups[, suspects = getParentSuspList()$name], "tp", fGroupsTPs = fGroups,
                                        MSPeakLists = mslists, ignoreParents = TRUE)
    componentsTPs <- filterComponentsForParents(fGroups, componentsTPs, TRUE)
    fGroups <- fGroups[, results = componentsTPs]
    
    # NOTE: phenazone is not found here and requires a tolerance of >=8 ppm. The high m/z deviation is most likely due to
    # the very high feature intensity of phenazone.
    formulas <- getFormulas(fGroups, mslists)
    
    # HACK: based on the candidate found for M149_R515_2917 it was found that the adduct is most likely [M]+, instead of
    # the [M+H]+ used elsewhere
    adducts(fGroups)["M149_R515_2917"] <- "[M]+"
    
    mfScorings <- c("score", "fragScore", "metFusionScore", "individualMoNAScore")
    # NOTE: disable multiprocessing, as that often results in PubChem connection troubles
    withOpt(MP.maxProcs = 1, {
        compoundsUF <- getCompounds(fGroups, mslists, formulas, database = "pubchem",
                                    maxCandidatesToStop = 25000, scoreTypes = mfScorings, topMost = 5000,
                                    timeout = 600)
    })
    compounds <- filter(compoundsUF, elements = "Cl0-0Br0-0B0-0Si0-0P0-0")
    
    mslists <- filterMSPLForAnn(mslists, fGroups, compounds)
    
    return(list(
        fGroups = fGroups,
        componentsTPs = componentsTPs,
        mslists = mslists,
        formulas = formulas,
        compoundsUF = compoundsUF, compounds = compounds
    ))
}
