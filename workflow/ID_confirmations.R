getIDLConfirmations <- function()
{
    suspList <- getTPsConfirmSuspList()
    
    # NOTE: data was exported as follows:
    # setDAMethod(anaInfoToExport, "pathToDataAnalysisMethod.m")
    # recalibrarateDAFiles(anaInfoToExport)
    # convertMSFiles(anaInfo = anaInfoToExport, from = "bruker", to = "mzML", algorithm = "bruker", centroid = TRUE)
    
    anaInfo <- read.csv("misc/analysisInfo-confirm.csv")

    fLists <- findFeatures(anaInfo, "openms", noiseThrInt = 1000, chromSNR = 3, chromFWHM = 5, minFWHM = 1,
                           maxFWHM = 30)
    
    # Group and align features between analyses
    # NOTE: extend default group RT a bit to group samples with some drift
    fGroups <- groupFeatures(fLists, "openms", rtalign = FALSE, QT = TRUE, maxGroupRT = 24)
    
    # Basic rule based filtering
    # NOTE: removed RT filter for aminophenol
    fGroups <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 1000, relMinReplicateAbundance = 1,
                      maxReplicateIntRSD = 0.75, blankThreshold = 5, removeBlanks = TRUE)
    
    adducts(fGroups) <- rep("[M+H]+", length(fGroups))
    
    fGroups <- doScreenSuspects(fGroups, suspList, onlyHits = TRUE)
    
    bgMSMS <- getBGMSMSPeaks(anaInfo[anaInfo$group == "MQ", ])
    mslists <- getMSPL(fGroups, bgMSMS)

    # export peak lists for MB
    plPath <- file.path("output", "conf-peaklists")
    dir.create(plPath, showWarnings = FALSE)
    for (i in seq_len(nrow(screenInfo(fGroups))))
    {
        MSMS <- mslists[[screenInfo(fGroups)$group[i]]]$MSMS
        if (!is.null(MSMS))
            data.table::fwrite(MSMS, file.path(plPath, paste0(screenInfo(fGroups)$name[i], ".csv")))
    }
    
    # Phenazone repeated experiments at higher conc
    anaInfoRepeatPhena <- read.csv("misc/analysisInfo-confirm_phena_high.csv")
    
    fListRepeatPhena <- findFeatures(anaInfoRepeatPhena, "openms", noiseThrInt = 1000, chromSNR = 3,
                                     chromFWHM = 5, minFWHM = 1, maxFWHM = 30)
    
    # Group and align features between analyses
    # NOTE: extend default group RT a bit to group samples with some drift
    fGroupsRepeatPhena <- groupFeatures(fListRepeatPhena, "openms", rtalign = FALSE, QT = TRUE, maxGroupRT = 24)
    
    # Basic rule based filtering
    fGroupsRepeatPhena <- filter(fGroupsRepeatPhena, preAbsMinIntensity = 100, absMinIntensity = 1000,
                                  relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75, blankThreshold = 5,
                                  removeBlanks = TRUE, retentionRange = c(60 * 2, Inf), mzRange = NULL)
    
    adducts(fGroupsRepeatPhena) <- rep("[M+H]+", length(fGroupsRepeatPhena))
    fGroupsRepeatPhena <- doScreenSuspects(fGroupsRepeatPhena, suspList, onlyHits = TRUE)
    mslistsRepeatPhena <- getMSPL(fGroupsRepeatPhena, bgMSMS)
    
    confirmedIDLevels <- data.table::data.table(
        parent = c(
            "flecainide",
            "metoprolol",
            "sulfamethoxazole",
            "sulfamethoxazole",
            "sulfamethoxazole",
            "sulfamethoxazole",
            "sulfamethoxazole",
            "sulfamethoxazole",
            "sulfamethoxazole",
            "phenazone",
            "phenazone",
            "phenazone",
            "phenazone"
        ),
        group = c(
            "M166_R473_6824",
            "M226_R366_2919",
            "M99_R289_3167",
            "M99_R136_3553",
            "M110_R126_3769",
            "M110_R177_3768",
            "M174_R162_5656",
            "M174_R162_5656",
            "M94_R231_2478",
            "M375_R598_4727",
            "M122_R471_3775",
            "M150_R528_2607",
            "M94_R231_2478"
        ),
        SMILES = c(
            "CCNC(=O)C1=CC(=CC=C1)O",
            "NCC(COC1=CC=C(C=C1)CCOC)O",
            "Cc1cc(N)no1",
            "Cc1cc(N)no1",
            "C1(=CC=C(C=C1)N)O",
            "C1(=CC=C(C=C1)N)O",
            "Nc1ccc(cc1)S(O)(=O)=O",
            "NS(=O)(=O)c1ccc(O)cc1",
            "Nc1ccccc1",
            "CC1=C(C(=O)N(N1C)C2=CC=CC=C2)C3=C(N(N(C3=O)C4=CC=CC=C4)C)C",
            "C1=CC=C(C=C1)NC=O",
            "CCC(=O)NC1=CC=CC=C1",
            "Nc1ccccc1"
        ),
        InChIKey = c(
            "XEARKGOSBVSNJD-UHFFFAOYSA-N",
            "XJWXVDJGNOHFLR-UHFFFAOYSA-N",
            "FKPXGNGUVSHWQQ-UHFFFAOYSA-N",
            "FKPXGNGUVSHWQQ-UHFFFAOYSA-N",
            "PLIKAWJENQZMHA-UHFFFAOYSA-N",
            "PLIKAWJENQZMHA-UHFFFAOYSA-N",
            "HVBSAKJJOYLTQU-UHFFFAOYSA-N",
            "DIRCLGLKRZLKHG-UHFFFAOYSA-N",
            "PAYRUJLWNCNPSJ-UHFFFAOYSA-N",
            "ANYXUEIQQWKBQV-UHFFFAOYSA-N",
            "DYDNPESBYVVLBO-UHFFFAOYSA-N",
            "ZTHRQJQJODGZHV-UHFFFAOYSA-N",
            "PAYRUJLWNCNPSJ-UHFFFAOYSA-N"
        ),
        IDL = c(
            "disproved",
            "3a",
            "1",
            "disproved",
            "disproved",
            "disproved",
            "3a",
            "disproved",
            "1",
            "disproved",
            "1",
            "disproved",
            "1"
        )
    )
    
    return(list(
        suspList = suspList,
        confirmedIDLevels = confirmedIDLevels,
        fGroups = fGroups, fGroupsRepeatPhena = fGroupsRepeatPhena,
        mslists = mslists, mslistsRepeatPhena = mslistsRepeatPhena
    ))
}
