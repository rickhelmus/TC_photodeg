getSemiQuant <- function(featData, SuSData, unkData, unkFormData, unkCompData, suspListConfirm)
{
    # NOTE: data was exported as follows:
    # setDAMethod(anaInfoToExport, "pathToDataAnalysisMethod.m")
    # recalibrarateDAFiles(anaInfoToExport)
    # convertMSFiles(anaInfo = anaInfoToExport, from = "bruker", to = "mzML", algorithm = "bruker", centroid = TRUE)
    
    anaInfo <- read.csv("misc/analysisInfo-semi_quant.csv")
    
    fList <- findFeatures(anaInfo, "openms", noiseThrInt = 1000, chromSNR = 3, chromFWHM = 5, minFWHM = 1, maxFWHM = 30)
    
    # Group and align features between analyses
    # NOTE: extend default group RT a bit to group samples with some drift
    fGroups <- groupFeatures(fList, "openms", rtalign = FALSE, QT = TRUE, maxGroupRT = 24)
    
    # Basic rule based filtering
    fGroups <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 1000, relMinReplicateAbundance = 1,
                      maxReplicateIntRSD = 0.75, blankThreshold = 5, removeBlanks = TRUE,
                      retentionRange = c(60 * 2, Inf), mzRange = NULL)
    
    adducts(fGroups) <- rep("[M+H]+", length(fGroups))
    
    fGroups <- doScreenSuspects(fGroups, rbind(getParentSuspList(), suspListConfirm, fill = TRUE), onlyHits = TRUE)
    
    concTable <- data.table::fread("misc/concs.csv")
    
    quantTab <- getQuantCalibFromScreening(fGroups, concTable, areas = TRUE, average = TRUE)
    quantTab <- quantTab[!is.na(conc)] # UNDONE: bug?
    
    # omit maxed out results
    upperConcs <- c(
        flecainide = 75000,
        metoprolol = 18750,
        sulfamethoxazole = 18750,
        phenazone = 4688,
        "1-amino-3-[4-(2-methoxyethyl)phenoxy]propan-2-ol" = 18750,
        Aniline = 18750,
        "3-amino-5-methylisoxazole" = 37500,
        "p-aminophenol" = 150000,
        "Sulfanilic acid" = 150000,
        "4-hydroxybenzenesulfonamide" = 150000,
        Forbisen = 9375,
        Formanilide = 37500,
        "n-ethyl-3-hydroxybenzamide" = 37500,
        "n-phenylpropanamide" = 9375
    )
    quantTab <- quantTab[conc <= upperConcs[name]]
    quantTab <- quantTab[name != "p-aminophenol"] # peak elutes prior to waste switch
    quantTab <- quantTab[intensity > 0]
    
    quantTab[, c("RSQ", "residual") := {
        model <- summary(lm(intensity ~ conc, weights = 1/conc))
        slope <- model$coefficients[2]; intercept <- model$coefficients[1]
        list(
            model$r.squared,
            (intensity - (slope * conc + intercept)) / intensity * 100
        )
    }, by = "name"]
    quantTab <- quantTab[abs(residual) <= 20] # omit outliers
    quantTab <- quantTab[!name %in% c("flecainide", "metoprolol")] # these 2 really seem to affect model performance (poor R2)
    quantTab[, N := .N, by = "name"]
    quantTab <- quantTab[N >= 5]
    
    eluentTab <- data.table::data.table(time = c(0, 13*60, 14*60), B = c(5, 100, 100))
    
    fGroupsSuSQ <- predictRespFactors(SuSData$fGroups, eluent = eluentTab, organicModifier = "MeOH", pHAq = 3.41,
                                      calibrants = quantTab, concUnit = "uM", calibConcUnit = "ngL")
    fGroupsSuSQ <- calculateConcs(fGroupsSuSQ, areas = TRUE)
    
    compoundsUnkQ <- consensus(unkCompData$compoundsP, unkFormData$compoundsUnkForm)
    compoundsUnkQ <- predictRespFactors(compoundsUnkQ, featData$fGroupsSel, eluent = eluentTab,
                                        organicModifier = "MeOH", pHAq = 3.41, calibrants = quantTab, concUnit = "uM",
                                        calibConcUnit = "ngL")
    fGroupsUnkQ <- unkData$fGroups[, groupNames(compoundsUnkQ)]
    fGroupsUnkQ <- calculateConcs(fGroupsUnkQ, compoundsUnkQ, areas = TRUE)

    # for evaluating model properties    
    fGroups <- predictRespFactors(fGroups, eluent = eluentTab, organicModifier = "MeOH", pHAq = 3.41,
                                  calibrants = quantTab,
                                  concUnit = "ngL")
    fGroups <- calculateConcs(fGroups, areas = TRUE)
    
    return(list(
        fGroupsSuSQ = fGroupsSuSQ,
        fGroupsUnkQ = fGroupsUnkQ,
        fGroupsStd = fGroups
    ))
}
