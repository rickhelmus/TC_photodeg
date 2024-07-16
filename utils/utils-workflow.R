getParentSuspList <- function(prep = FALSE)
{
    ret <- data.table::data.table(
        name = c("flecainide", "metoprolol", "sulfamethoxazole", "phenazone"),
        SMILES = c("C1CCNC(C1)CNC(=O)C2=C(C=CC(=C2)OCC(F)(F)F)OCC(F)(F)F",
                   "CC(C)NCC(COC1=CC=C(C=C1)CCOC)O",
                   "CC1=CC(=NO1)NS(=O)(=O)C2=CC=C(C=C2)N",
                   "CC1=CC(=O)N(N1C)C2=CC=CC=C2"),
        rt = c(8.3 * 60, 6.8 * 60, 7.5 * 60, 7.6 * 60),
        # HACK: add mz fragments for phenazone since it's not annotated due to high m/z deviation (intensity overload)
        fragments_mz = c(NA, NA, NA, "147.0915;161.1072;174.0786"),
        fragments_formula = c("C17H18F6NO3;C6H12N;C11H7F6O3",
                              "C12H15O2;C9H9O;C6H14NO",
                              "C6H6NO2S;C6H6N;C4H7N2O;C6H6NO2S;O2S",
                              # from MB.eu (others MF)
                              "C9H11N2;C10H13N2C10H10N2O"))
    if (prep)
        ret <- patRoon:::prepareSuspectList(ret, NULL, TRUE, TRUE, TRUE, FALSE)
    return(ret)
}

getTPsConfirmSuspList <- function()
{
    # Standards of TPs
    data.table::data.table(
        name = c(
            "1-amino-3-[4-(2-methoxyethyl)phenoxy]propan-2-ol",
            "Aniline",
            "3-amino-5-methylisoxazole",
            "p-aminophenol",
            "Sulfanilic acid",
            "4-hydroxybenzenesulfonamide",
            "Forbisen",
            "Formanilide",
            "n-ethyl-3-hydroxybenzamide",
            "n-phenylpropanamide"
        ),
        SMILES = c(
            "NCC(COC1=CC=C(C=C1)CCOC)O",
            "Nc1ccccc1",
            "Cc1cc(N)no1",
            "C1(=CC=C(C=C1)N)O",
            "Nc1ccc(cc1)S(O)(=O)=O",
            "NS(=O)(=O)c1ccc(O)cc1",
            "CC1=C(C(=O)N(N1C)C2=CC=CC=C2)C3=C(N(N(C3=O)C4=CC=CC=C4)C)C",
            "C1=CC=C(C=C1)NC=O",
            "CCNC(=O)C1=CC(=CC=C1)O",
            "CCC(=O)NC1=CC=CC=C1"
        ),
        rt = c(366, 235, 284, 114, 170, 291, 562, 474, 434, 534)
    )
}

doScreenSuspects <- function(fg, ...) screenSuspects(fg, rtWindow = 12, mzWindow = 0.005, ...)

getRegTabs <- function(fGroups, minRSQ, filter = TRUE, minRepl = 2)
{
    fGroups <- fGroups[, rGroups = grep("^UV", replicateGroups(fGroups), value = TRUE)]
    # NOTE: normalizing is not really needed to get the linearity properties, but getting the normalized values is handy
    # for reporting...
    fGroups <- normInts(fGroups, groupNorm = TRUE)
    
    regTabs <- sapply(c("UV", "UV_H2O2", "UV_H2O2_NOM"), function(samp)
    {
        concs <- c(25, 75, 150)
        fg <- fGroups[, rGroups = paste0(samp, "_", concs)]
        fg <- filter(fg, absMinReplicates = minRepl)
        tab <- as.data.table(fg, regression = TRUE, average = FALSE, normalize = TRUE)
        if (is.null(tab[["RSQ"]]))
            tab[, c("RSQ", "p", "slope") := NA_real_]
        if (filter)
            tab <- tab[RSQ >= minRSQ & p <= 0.05 & slope > 0]
        
        if (nrow(tab) > 0)
        {
            # ensure highest concentration is part of regression
            # NOTE: combined with the min replicates filter above, this effectively means that the feature should be in the two highest mixes
            concs <- analysisInfo(fg)$conc
            tab[, maxConc := {
                wh <- unlist(.SD) != 0
                max(concs[wh])
            }, by = seq_len(nrow(tab)), .SDcols = analyses(fg)]
            
            if (filter)
                tab <- tab[maxConc == max(concs)]
        }
        
        return(tab)
    }, simplify = FALSE)
    
    return(regTabs)
}

getMSPL <- function(fg, bgMSMS, topMSMSPeaks = 25, ...)
{
    avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
    ret <- generateMSPeakLists(fg, "mzr", maxMSRtWindow = 5, precursorMzWindow = 4, avgFeatParams = avgMSListParams,
                               avgFGroupParams = avgMSListParams, ...)
    ret <- filter(ret, relMSMSIntThr = 0.05, topMSMSPeaks = topMSMSPeaks)
    ret <- filterMSPLForBlanks(ret, bgMSMS)
    
    ret <- delete(ret, j = function(pl, grp, ana, t)
    {
        if (t == "MS")
            return(FALSE)
        if (!is.null(ana))
        {
            precMZ <- ret[[ana, grp]]$MS[precursor == TRUE]$mz
            return((pl$mz - precMZ) > 4)
        }
        precRow <- ret[[grp]]$MS[precursor == TRUE]
        return((pl$mz - precRow[precursor == TRUE]$mz) > 4)
    })
    
    return(ret)
}

getFormulas <- function(fg, mslists, elements = "CHNOSF", relMzDev = 5, ...)
{
    ret <- generateFormulas(fg, mslists, "genform", relMzDev = relMzDev, elements = elements, oc = TRUE,
                            calculateFeatures = FALSE, featThresholdAnn = 0.75, ...)
    return(ret)
}

getBGMSMSPeaks <- function(anaInfo)
{
    blSpecs <- data.table::rbindlist(future.apply::future_Map(anaInfo$analysis, anaInfo$path, f = function(ana, path)
    {
        printf("Running sample %d/%d\n", match(ana, anaInfo$analysis), nrow(anaInfo))
        specs <- patRoon:::loadSpectra(file.path(path, paste0(ana, ".mzML")))
        specsMSMSHd <- specs$header[specs$header$msLevel == 2 & specs$header$retentionTime >= 120 & specs$header$basePeakIntensity >= 5000]
        specsMSMS <- data.table::rbindlist(lapply(specs$spectra[specsMSMSHd$seqNum], function(x)
        {
            x <- as.data.table(x)
            x <- x[intensity >= 1000]
            if (nrow(x) > 25)
            {
                oint <- order(x$intensity, decreasing = TRUE)
                x <- x[oint[seq_len(25)]]
            }
            return(x)
        }))
        distm <- dist(specsMSMS$mz)
        hc <- fastcluster::hclust(distm, "ward.D2")
        specsMSMS[, cl := cutree(hc, h = 0.003)]
        rm(distm, hc)
        gc()
        specsMSMSUn <- data.table::copy(specsMSMS)
        specsMSMSUn[, c("mz", "intensity", "specs") := .(mean(mz), mean(intensity), .N), by = cl]
        specsMSMSUn[, specsRel := specs / nrow(specsMSMSHd)]
        specsMSMSUn <- unique(specsMSMSUn, by = "cl")
        return(specsMSMSUn)
    }), idcol = TRUE)
    distm <- dist(blSpecs$mz)
    hc <- fastcluster::hclust(distm, "ward.D2")
    blSpecsCL <- data.table::copy(blSpecs)
    blSpecsCL[, cl := cutree(hc, h = 0.003)]
    blSpecsCL[, clSize := .N, by = "cl"]
    blSpecsCL[, clAnaSize := data.table::uniqueN(.id), by = "cl"]
    blSpecsCL[, c("mz", "intensity", "specs", "specsRel") := .(mean(mz), mean(intensity), sum(specs),
                                                               sum(specsRel) / nrow(anaInfo)), by = cl]
    blSpecsUn <- unique(blSpecsCL, by = "cl")
    blSpecsUn <- blSpecsUn[clAnaSize >= (0.8*nrow(anaInfo))]
    blSpecsUn[, intensity := intensity / max(intensity)]
    return(blSpecsUn[specsRel >= 0.1])
}

filterMSPLForBlanks <- function(mslists, bgMSMS)
{
    delete(mslists, j = function(pl, grp, ana, t)
    {
        if (t == "MS")
            return(FALSE)
        if (!is.null(ana))
            return(sapply(pl$mz, function(mz) any(abs(mz - bgMSMS$mz) <= 0.003)))
        
        return(sapply(pl$mz, function(mz) any(abs(mz - bgMSMS$mz) <= 0.003)))
    })
}

filterMSPLForAnn <- function(mslists, fGroups, comps = NULL)
{
    formsBG <- getFormulas(fGroups, mslists, elements = "CHNOPSF")
    return(filter(mslists, annotatedBy = if (!is.null(comps)) list(formsBG, comps) else formsBG))
}

filterFormulasForSusps <- function(fg, forms)
{
    delete(forms, j = function(ann, grp)
    {
        # remove candidates with F/S unless flecainide/sulphamethoxazole
        suspLine <- screenInfo(fg)[group == grp]
        if (nrow(suspLine) == 0)
            return(FALSE)
        rmElements <- character()
        if (!any(grepl("fleca", suspLine$name)))
            rmElements <- c(rmElements, "F")
        if (!any(grepl("sulfa", suspLine$name)))
            rmElements <- c(rmElements, "S")
        if (length(rmElements) == 0)
            return(FALSE)
        return(grepl(paste0(rmElements, collapse = "|"), ann$neutral_formula))
    })
}

filterComponentsForParents <- function(fGroups, components, otherOK)
{
    delete(components, j = function(ct, cn)
    {
        parFG <- componentInfo(components)[match(cn, name)]$parent_group
        par <- screenInfo(fGroups)[match(parFG, group)]$name
        par <- strtrim(par, 4) # to match with shortened names in replicate groups
        
        getParRGs <- function(p) grep(p, replicateGroups(fGroups), fixed = TRUE, value = TRUE)
        
        parRGs <- getParRGs(par)
        otherParRGs <- unlist(lapply(c("fleca", "meto", "sulfa", "phena"), getParRGs))
        otherParRGs <- setdiff(otherParRGs, parRGs)
        
        expRGs <- setdiff(replicateGroups(fGroups), parRGs)
        expRGs <- setdiff(expRGs, otherParRGs)
        expRGs <- expRGs[!grepl("MIX", expRGs, fixed = TRUE)]
        
        maxIntForRGs <- function(gn, rgs)
        {
            ints <- as.data.table(fGroups, features = TRUE)[group == gn][replicate_group %in% rgs]$intensity
            return(if (length(ints) == 0) 0 else max(ints))
        }
        
        keep <- sapply(ct$group, function(gn)
        {
            if (maxIntForRGs(gn, parRGs) > 0)
                return(TRUE)
            if (!otherOK)
                return(FALSE)
            mopi <- maxIntForRGs(gn, otherParRGs)
            return(mopi == 0 | (maxIntForRGs(gn, expRGs) / mopi >= 5))
        })
        
        return(!keep)
    })
}

getCompounds <- function(fg, mslists, formulas, ...)
{
    ret <- generateCompounds(fg, mslists, "metfrag", dbRelMzDev = 8, fragRelMzDev = 8, fragAbsMzDev = 0.002,
                             ...)
    ret <- addFormulaScoring(ret, formulas, updateScore = TRUE)
    return(ret)
}

estIDLevels <- function(formTab, compTab, formRank, compRank, annSimForm, annSimComp)
{
    # simplified version of patRoon:::estimateIdentificationLevel()
    # ID levels are as IDLevelRules.yml
    
    compRow <- compTab[compRank]
    
    # level 2a
    ims <- compRow[["individualMoNAScore"]]
    if (compRank == 1 && !is.null(ims) && ims >= 0.9)
    {
        imsNext <- if (nrow(compTab) > 1) max(compTab[["individualMoNAScore"]][-compRank]) else 0
        if (is.na(imsNext) || imsNext == 0)
            return("2a")
    }
    
    # level 3b
    if (compRank == 1 && !is.na(annSimComp) && annSimComp >= 0.9)
    {
        fs <- compTab$fragScore
        fs <- fs / max(fs) # normalize
        fsNext <- if (nrow(compTab) > 1) max(fs[-compRank]) else 0
        if (fs[compRank] >= 0.9 && (nrow(compTab) == 1 || fsNext == 0))
            return("3b")
    }
    
    # level 3c
    if (!is.null(ims) && ims >= 0.7)
        return("3c")
    
    # 3d
    if (!is.na(annSimComp) && annSimComp >= 0.7)
        return("3d")
    
    # 4
    if (!is.na(formRank) && formRank == 1)
    {
        formRow <- formTab[formRank]
        iso <- formRow$isoScore
        if (is.null(iso) || is.na(iso))
            return("5")
        
        # modification
        if (!is.na(annSimForm) && annSimForm >= 0.9 && iso >= 0.5)
            return("4b")
        
        isoNext <- if (nrow(formTab) > 1) max(formTab$isoScore[-formRank]) else -Inf
        if (!is.na(isoNext) && (iso - isoNext) < 0.2)
            return("5")
        if (!is.na(annSimForm) && annSimForm >= 0.7 && iso >= 0.5)
            return("4a")
        if (iso >= 0.9)
            return("4c")
    }
    
    return("5")
}
