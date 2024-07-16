getFeatures <- function()
{
    # NOTE: data was exported as follows:
    # setDAMethod(anaInfoToExport, "pathToDataAnalysisMethod.m")
    # recalibrarateDAFiles(anaInfoToExport)
    # convertMSFiles(anaInfo = anaInfoToExport, from = "bruker", to = "mzML", algorithm = "bruker", centroid = TRUE)
    
    anaInfo <- read.csv("misc/analysisInfo.csv", colClasses = c(conc = "numeric"))
    
    # Find all features
    fList <- findFeatures(anaInfo, "openms", noiseThrInt = 1000, chromSNR = 3, chromFWHM = 5, minFWHM = 1, maxFWHM = 30)
    
    # Group and align features between analyses
    # NOTE: extend default group RT a bit to group samples with some drift
    fGroups <- groupFeatures(fList, "openms", rtalign = FALSE, QT = TRUE, maxGroupRT = 24)
    
    # Basic rule based filtering
    fGroupsF <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 1000, relMinReplicateAbundance = 1,
                       maxReplicateIntRSD = 0.75, blankThreshold = 5, removeBlanks = TRUE,
                       retentionRange = c(60 * 2, Inf), mzRange = NULL)
    
    
    # regression prioritization
    fGroupsScrPar <- doScreenSuspects(fGroupsF, getParentSuspList(), onlyHits = TRUE, adduct = "[M+H]+")
    regTabs <- getRegTabs(fGroupsF, minRSQ = 0.5)
    regTabsGN <- unique(unlist(lapply(regTabs, "[[", "group")))
    fGroupsP <- fGroupsF[, union(regTabsGN, screenInfo(fGroupsScrPar)$group)] # only keep prioritized/parent features
    
    # Perform automatic generation of components
    components <- generateComponents(fGroupsF, "cliquems", ppm = 7, ionization = "positive")
    # Simplify components for manual processing with checkComponents() below: remove any features without annotations and
    # remove any components that do not contain any of the prioritized features
    componentsF <- delete(components, j = function(ct, ...) is.na(ct$isonr) & is.na(ct$adduct_ion))
    componentsF <- delete(componentsF, j = function(ct, ...) !any(ct$group %in% names(fGroupsP)))
    # we may end up with 1-sized components, remove these as these cannot be used for annotation
    componentsF <- filter(componentsF, size = c(2, 1E6))
    
    # manually check/filter components
    # checkComponents(componentsF, fGroupsF, session = "misc/components_check.yml")
    componentsF <- filter(componentsF, checkComponentsSession = "misc/components_check.yml")
    
    # select/annotate features
    fGroupsIonSel <- selectIons(fGroupsP, componentsF, prefAdduct = "[M+H]+", onlyMonoIso = TRUE)
    
    # manually check/filter features
    # checkFeatures(fGroupsIonSel, session = "misc/features_check.yml")
    fGroupsSel <- filter(fGroupsIonSel, checkFeaturesSession = "misc/features_check.yml")
    
    # Remove analyses of high standards: the intensities in these analyses were severely overloaded, which therefore
    # negatively impact the mass accuracy for peak lists for the parents obtained later
    fGroupsSel <- fGroupsSel[analysisInfo(fGroupsSel)$conc <= 25 | !grepl("^MIX", analysisInfo(fGroupsSel)$group)]
 
    return(list(
        fGroupsRaw = fGroups,
        fGroupsFiltered = fGroupsF,
        fGroupsPrioReg = fGroupsP,
        fGroupsIonSel = fGroupsIonSel,
        fGroupsSel = fGroupsSel
    ))
}
