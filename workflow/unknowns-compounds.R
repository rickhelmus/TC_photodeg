getUnknownsCompounds <- function(unkData, SuSData)
{
    # perform compound TP score calculations and thresholding. This is quite computationally intensive and will take up
    # several hours
    compUTabs <- Map(screenInfo(unkData$fGroups)$group, screenInfo(unkData$fGroups)$formula,
                     screenInfo(unkData$fGroups)$SMILES, products(SuSData$TPsCons)[screenInfo(unkData$fGroups)$name],
                     componentTable(unkData$componentsTPs)[componentInfo(unkData$componentsTPs)[match(screenInfo(unkData$fGroups)$group, parent_group)]$name],
                     componentTable(SuSData$componentsTPs)[componentInfo(SuSData$componentsTPs)[match(screenInfo(unkData$fGroups)$group, parent_group)]$name],
                     MoreArgs = list(fGroups = unkData$fGroups, mslists = unkData$mslists,
                                     formulas = unkData$formulas, compounds = unkData$compounds,
                                     minFormulaFit = 0.94, minCompoundFit = 0.54, minMaxTPSim = 0.65),
                     f = calcTPCompoundScores)
    
    # save data to avoid the need to repeat calculations in new R sessions
    saveRDS(compUTabs, "output/compUTabs.Rds")

    # inspect results
    if (FALSE)
    {
        compUTab <- makeUnknownCompTab(unkData$fGroups, compUTabs)
        getUnknownCompTabReact(compUTab)
    }
    
    # selected candidates
    selectedCandidates <- list(
        flecainide = data.table::data.table(
            group = character(),
            InChIKey = character()
        ),
        metoprolol = data.table::data.table(
            group = character(),
            InChIKey = character()
        ),
        sulfamethoxazole = data.table::data.table(
            group = c(
                "M140_R170_3386",
                "M190_R159_4149",
                "M190_R159_4149",
                "M190_R159_4149"
            ),
            InChIKey = c(
                "DHIMKPLQNUVPHL-UHFFFAOYSA-N",
                "HCEYSAVOFADVMD-UHFFFAOYSA-N",
                "MSYVHZXMLZBEIF-UHFFFAOYSA-N",
                "DJOHZBQQLCXFQU-UHFFFAOYSA-N"
            )
        ),
        phenazone = data.table::data.table(
            group = c(
                "M176_R533_3000",
                "M219_R515_7445",
                "M375_R598_4727"
            ),
            InChIKey = c(
                "CEULXHFVIAILQG-UHFFFAOYSA-N",
                "FACUXNKAWFZFEV-UHFFFAOYSA-N",
                "ANYXUEIQQWKBQV-UHFFFAOYSA-N"
            )
        )
    )
    # add more metadata for eg reporting
    selectedCandidates <- lapply(selectedCandidates, function(sel)
    {
        sel[, SMILES := {
            IKs <- InChIKey
            as.data.table(unkData$compounds)[match(IKs, InChIKey)]$SMILES
        }]
    })
    
    # only keep selected structs to get more representative compound ranking for identification level calculations
    compoundsP <- delete(unkData$compounds,
                         j = function(ann, grp) !ann$InChIKey %in% unlist(lapply(selectedCandidates, "[[", "InChIKey")))
    
    selectedCandidates <- Map(names(selectedCandidates), selectedCandidates, f = function(parn, sel)
    {
        if (nrow(sel) == 0)
            return(sel)
        parg <- screenInfo(unkData$fGroups)[name == parn]$group
        ct <- compUTabs[[parg]]
        sel <- merge(sel, ct[, c("InChIKey", "neutral_formula", grep("^annSim", names(ct), value = TRUE), "compFit",
                                 "formFit", "maxTPSimAllPred", "TPScore"), with = FALSE],
                     by = "InChIKey", sort = FALSE)
        
        sel[, annSimForm := mapply(group, neutral_formula, FUN = calcAnnSim,
                                   MoreArgs = list(featAnn = unkData$formulas, mslists = unkData$mslists, formulas = NULL))]
        sel[, formRank := mapply(group, neutral_formula, FUN = function(grp, form)
        {
            if (!is.null(unkData$formulas[[grp]]))
                return(match(form, unkData$formulas[[grp]]$neutral_formula, nomatch = NA_integer_))
            return(NA_integer_)
        })]
        sel[, compRank := mapply(group, InChIKey, FUN = function(grp, IK)
        {
            return(match(IK, compoundsP[[grp]]$InChIKey, nomatch = NA_integer_))
        })]
        sel[, estIDLevel := mapply(estIDLevels, annotations(unkData$formulas)[group], annotations(compoundsP)[group],
                                   formRank, compRank, annSimForm, annSim)]
        
        return(sel)
    })
    
    return(list(
        selectedCandidates = selectedCandidates,
        compoundsP = compoundsP
    ))
}
