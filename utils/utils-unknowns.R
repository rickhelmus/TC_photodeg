calcAnnSim <- function(featAnn, mslists, grp, UID, formulas = NULL)
{
    if (is.null(featAnn[[grp]]) || !UID %in% featAnn[[grp]]$UID)
        return(NA_real_)
    aplArgs <- list(featAnn, groupName = grp, index = which(featAnn[[grp]]$UID == UID), MSPeakLists = mslists)
    if (!is.null(formulas))
        aplArgs <- c(aplArgs, list(formulas = formulas))
    annPL <- do.call(annotatedPeakList, aplArgs)
    if (is.null(annPL))
        return(NA_real_)
    if (nrow(annPL) == 0)
        return(NA_real_)
    patRoon:::annotatedMSMSSimilarity(annPL, getDefSpecSimParams(removePrecursor = TRUE))
}

calcFormulaFit <- function(formula1, formula2)
{
    calcOneWay <- function(f1, f2)
    {
        f1ElCount <- tryCatch(sum(patRoon:::splitFormulaToList(f1)), error = function(...) NULL)
        fSub <- patRoon:::splitFormulaToList(patRoon:::subtractFormula(f2, f1))
        f1ElMissing <- sum(abs(fSub[fSub < 0]))
        return((f1ElCount - f1ElMissing) / f1ElCount)
    }
    return(max(calcOneWay(formula1, formula2), calcOneWay(formula2, formula1)))
}

getFMCS <- function(SMILES1, SMILES2, au = 1, bu = 4, matching.mode = "aromatic", CIDs = NULL, ...)
{
    # ChemmineR cannot do SMILES conversion on Windows, so first convert to a temporary SDF file with RCDK. This also
    # makes aromatic rings look better...
    
    mols <- patRoon:::getMoleculesFromSMILES(c(SMILES1, SMILES2), doTyping = TRUE, emptyIfFails = FALSE)
    if (!patRoon:::isValidMol(mols[[1]]) || !patRoon:::isValidMol(mols[[2]]))
        return(NULL) # UNDONE: or zero?
    
    mols <- lapply(mols, rcdk::generate.2d.coordinates)
    
    sdfFile <- tempfile(fileext = ".sdf")
    rcdk::write.molecules(mols, sdfFile, together = TRUE)
    
    SDFSet <- ChemmineR::read.SDFset(sdfFile)
    if (!is.null(CIDs))
        ChemmineR::cid(SDFSet) <- CIDs
    return(fmcsR::fmcs(SDFSet[[1]], SDFSet[[2]], au = au, bu = bu, matching.mode = matching.mode, ...))
}

calcStructFitFMCS <- function(SMILES1, SMILES2, ...)
{
    return(unname(fmcsR::stats(getFMCS(SMILES1, SMILES2, ...))["Overlap_Coefficient"]))
}

calcMaxSim <- function(targetSMILES, otherSMILES, fpType = "extended", fpSimMethod = "tanimoto")
{
    dists <- sapply(otherSMILES, patRoon:::distSMILES, SMI1 = targetSMILES, fpType = fpType, fpSimMethod = fpSimMethod)
    wh <- which.max(dists)
    return(list(maxSim = dists[wh], SMI = otherSMILES[wh]))
}

calcTPScoreComp <- function(fit, maxTPSim, annSim)
{
    maxOrZero(fit, maxTPSim) + maxOrZero(annSim)
}

calcTPFormulaScores <- function(fGroups, parentGroup, parentFormula, mslists, formulas, componTab,
                                minAnnSim = 0.5, minFormulaFit = 0.7)
{
    ret <- as.data.table(formulas[componTab$group])
    
    # omit S/F unless present in parent
    if (!grepl("S", parentFormula))
        ret <- ret[!grepl("S", neutral_formula)]
    if (!grepl("F", parentFormula))
        ret <- ret[!grepl("F", neutral_formula)]
    
    printf("Processing %d formulas for parent %s\n", nrow(ret), parentGroup)
    
    cat("Doing formula fits...\n")
    ret[, formFit := sapply(neutral_formula, calcFormulaFit, parentFormula)]
    ret <- ret[formFit >= minFormulaFit]
    printf("Remaining formulas after formFit filter: %d\n", nrow(ret))
    
    cat("Calculating annotation similarities...\n")
    ret[, annSim := mapply(group, UID, FUN = calcAnnSim, MoreArgs = list(featAnn = formulas, mslists = mslists,
                                                                         formulas = NULL))]
    
    ret <- ret[is.na(annSim) | annSim >= minAnnSim]
    printf("Remaining formulas after annSim filter: %d\n", nrow(ret))
    
    cat("Estimating ID levels...\n")
    
    if (nrow(ret) == 0)
    {
        ret[, fragNLMatches := numeric()]
        ret[, TPScore := numeric()]
    }
    else
    {
        cat("Doing annotation matches...\n")
        if (parentFormula %in% formulas[[parentGroup]]$UID)
        {
            parentFI <- formulas[[parentGroup]][UID == parentFormula]$fragInfo[[1]]
            getAnnMatches <- function(fi, col) sum(parentFI[[col]] %in% fi[[col]])
            ret[, fragNLMatches := {
                uid <- UID
                fi <- formulas[[group]][UID == uid]$fragInfo[[1]]
                sum(getAnnMatches(fi, "ion_formula"), getAnnMatches(fi, "neutral_loss"))
            }, by = c("group", "UID")]
        }
        else
            ret[, fragNLMatches := NA_real_]
        ret[, TPScore := formFit + maxOrZero(annSim), by = seq_len(nrow(ret))]
    }
    
    ret[, formulaDiff := sapply(neutral_formula, patRoon:::subtractFormula, formula2 = parentFormula)]
    
    return(ret[])
}

calcTPCompoundScores <- function(fGroups, parentGroup, parentFormula, parentSMILES, mslists, formulas, compounds,
                                 TPsTab, componTab, componStructTab, minFormulaFit, minCompoundFit, minMaxTPSim,
                                 minXLogPDiff = 1, minRTDiff = 30, onlyGroup = NULL)
{
    ret <- as.data.table(compounds[componTab$group])
    if (!is.null(onlyGroup))
    {
        ret <- ret[group == onlyGroup]
        if (nrow(ret) == 0)
            return(data.table())
    }
    
    printf("Processing %d compounds for parent %s\n", nrow(ret), parentGroup)
    
    cat("Doing RT dirs...\n")
    
    parentLogP <- getXLogP(parentSMILES)
    parentCompRes <- compounds[[parentGroup]][SMILES == parentSMILES]
    gInfo <- groupInfo(fGroups)
    ret[, XLogP_RCDK := future.apply::future_sapply(SMILES, getXLogP)]
    ret[, XLogPDiffParent := XLogP_RCDK - parentLogP]
    ret[, retDiffParent := gInfo[group, "rts"] - gInfo[parentGroup, "rts"]]
    ret[, retDirExpect := data.table::fcase((XLogPDiffParent + minXLogPDiff) < 0, -1,
                                            (XLogPDiffParent - minXLogPDiff) > 0, 1,
                                            default = 0)]
    ret[, retDir := data.table::fcase((retDiffParent + minRTDiff) < 0, -1,
                                      (retDiffParent - minRTDiff) > 0, 1,
                                      default = 0)]
    ret <- ret[retDirExpect == retDir | retDirExpect == 0 | retDir == 0]
    printf("Remaining compounds after retDir filter: %d\n", nrow(ret))
    
    cat("Doing formula fits...\n")
    ret[, formFit := sapply(neutral_formula, calcFormulaFit, parentFormula)]
    ret <- ret[formFit >= minFormulaFit]
    printf("Remaining compounds after formFit filter: %d\n", nrow(ret))
    
    cat("Doing compound FMCS...\n")
    ret[, compFit := future.apply::future_sapply(SMILES, calcStructFitFMCS, parentSMILES, au = 0, bu = 0)]
    
    cat("Doing TP sims...\n")
    ret[, c("maxTPSimHits", "maxTPSimHitSMILES") := .(NA_real_, NA_real_)] # disable for now: not used and saves some time
    ret[, c("maxTPSimAllPred", "maxTPSimAllPredSMILES") :=
            data.table::rbindlist(future.apply::future_lapply(SMILES, calcMaxSim, TPsTab$SMILES))]
    ret <- ret[compFit >= minCompoundFit | maxTPSimAllPred >= minMaxTPSim]
    printf("Remaining compounds after compFit/maxTPSim filter: %d\n", nrow(ret))
    
    cat("Calculating annotation similarities...\n")
    ret[, annSim := mapply(group, UID, FUN = calcAnnSim, MoreArgs = list(featAnn = compounds, mslists = mslists))]
    ret[, annSimBoth := mapply(group, UID, FUN = calcAnnSim, MoreArgs = list(featAnn = compounds, mslists = mslists, formulas = formulas))]
    
    if (nrow(ret) == 0)
    {
        ret[, c("maxTPSimHits", "maxTPSimHitSMILES", "maxTPSimAllPred", "maxTPSimAllPredSMILES") :=
                .(numeric(), character(), numeric(), character())]
        ret[, fragNLMatches := numeric()]
        ret[, TPScore := numeric()]
    }
    else
    {
        cat("Doing annotation matches...\n")
        parentFI <- parentCompRes$fragInfo[[1]]
        getAnnMatches <- function(fi, col) sum(parentFI[[col]] %in% fi[[col]])
        ret[, fragNLMatches := {
            uid <- UID
            fi <- compounds[[group]][UID == uid]$fragInfo[[1]]
            sum(getAnnMatches(fi, "ion_formula"), getAnnMatches(fi, "neutral_loss"))
        }, by = c("group", "UID")]
        
        ret[, TPScore := calcTPScoreComp(compFit, maxTPSimAllPred, annSim), by = seq_len(nrow(ret))]
    }
    
    return(ret[])
}

embedPlot <- function(code, width = 350, height = 350, res = 72, style = "")
{
    f <- tempfile(fileext = ".png")
    withr::with_png(f, code, width = width, height = height, res = res)
    sprintf("<img src=%s style='%s'></img>", knitr::image_uri(f), style)
}

makeUnknownFormTab <- function(tabs, topMost = 25)
{
    tabs <- lapply(tabs, function(x)
    {
        x <- data.table::copy(x)
        x[, ord := frank(-TPScore, ties.method = "first"), by = "group"]
        return(x[ord <= topMost])
    })
    tab <- data.table::rbindlist(tabs, idcol = "parent_group", use.names = TRUE)
    tab <- tab[, c("parent_group", "group", "neutral_formula", "TPScore", "formulaDiff", "formFit",
                   grep("annSim", names(tab), value = TRUE), "fragNLMatches", "neutralMass"),
               with = FALSE]
    return(tab)
}

getUnknownFormTabReact <- function(tab)
{
    reactable::reactable(tab, groupBy = c("parent_group", "group"), filterable = TRUE)
}


makeUnknownCompTab <- function(fGroupsUnk, tabs, topMost = 25)
{
    tabs <- lapply(tabs, function(x)
    {
        x <- data.table::copy(x)
        x[, ord := data.table::frank(-TPScore, ties.method = "first"), by = "group"]
        return(x[ord <= topMost])
    })
    ret <- data.table::rbindlist(tabs, idcol = "parent_group", use.names = TRUE)
    ret[, parent_SMILES := screenInfo(fGroupsUnk)[match(parent_group, group)]$SMILES]
    ret[, TP := {
        embedPlot({
            fmcs <- getFMCS(parent_SMILES, SMILES[1], au = 0, bu = 0)
            if (is.null(fmcs))
                patRoon:::noDataPlot()
            else
            {
                par(mfrow = c(1, 2))
                plot(fmcsR::mcs1(fmcs)[[1]][[1]], colbonds = fmcsR::mcs1(fmcs)[[2]][[1]])
                plot(fmcsR::mcs2(fmcs)[[1]][[1]], colbonds = fmcsR::mcs2(fmcs)[[2]][[1]])
            }
        }, width = 750, height = 500, style = "width: 350px; height: auto;")
    }, by = c("parent_SMILES", "SMILES")]
    ret[, TP_MSAP := {
        embedPlot({
            plot(patRoon:::getRCDKStructurePlot(rcdk::parse.smiles(maxTPSimAllPredSMILES[1])[[1]]))
        }, width = 500, height = 500, style = "width: 200px; height: auto;")
    }, by = c("maxTPSimAllPredSMILES")]
    
    
    numCols <- names(which(sapply(ret, is.numeric)))
    ret[, (numCols) := lapply(mget(numCols), round, 2)]
    return(ret)
}

getUnknownCompTabReact <- function(tab)
{
    reactable::reactable(tab[, c("parent_group", "group", "TP", "TPScore", "InChIKey", "XLogPDiffParent",
                                 "retDiffParent", "retDirExpect", "retDir", "formFit", "compFit",
                                 "maxTPSimAllPred", "TP_MSAP",
                                 grep("annSim", names(tab), value = TRUE),
                                 "fragNLMatches", "neutralMass", "neutral_formula"),
                             with = FALSE], groupBy = c("parent_group", "group"), columns = list(
                                 TP = reactable::colDef(html = TRUE, width = 350),
                                 TP_MSAP = reactable::colDef(html = TRUE, width = 200)
                             ))
}
