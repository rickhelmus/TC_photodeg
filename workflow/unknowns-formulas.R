getUnknownsFormulas <- function(unkData, SuSData)
{
    # do TP score calculations and thresholding
    formUTabs <- Map(screenInfo(unkData$fGroups)$group, screenInfo(unkData$fGroups)$formula, 
                     componentTable(unkData$componentsTPs)[componentInfo(unkData$componentsTPs)[match(screenInfo(unkData$fGroups)$group, parent_group)]$name],
                     MoreArgs = list(fGroups = unkData$fGroups, mslists = unkData$mslists, formulas = unkData$formulas,
                                     minAnnSim = 0, minFormulaFit = 0.94),
                     f = calcTPFormulaScores)

    if (FALSE)
    {
        # inspect results
        formUTab <- makeUnknownFormTab(formUTabs)
        getUnknownFormTabReact(formUTab)
    }
    
    # selections
    selectedCandidates <- list(
        flecainide = data.table::data.table(
            group = c(
                "M166_R473_6824"
            ),
            formula = c(
                "C9H11NO2"
            ),
            SMILES = c(
                "CCNC(=O)C1=CC(=CC=C1)O"
            )
        ),
        metoprolol = data.table::data.table(
            group = character(), formula = character(), SMILES = character()
        ),
        sulfamethoxazole = data.table::data.table(
            group = character(), formula = character(), SMILES = character()
        ),
        phenazone = data.table::data.table(
            group = c(
                "M122_R471_3775",
                "M149_R515_2917",
                "M150_R528_2607", # in compounds but didn't make top 25
                "M165_R462_6907",
                "M166_R473_6824",
                "M208_R533_2918" # unsure what structure this could be
            ),
            formula = c(
                "C7H7NO",
                "C8H9N2O", # note that the structure is charged, see below
                "C9H11NO",
                "C8H8N2O2",
                "C9H11NO2",
                "C11H13NO3"
            ),
            SMILES = c(
                "C1=CC=C(C=C1)NC=O", # is in PubChem, but MS/MS data is poor
                # HACK: note that the structure is already charged, and therefore the adduct should be M+ instead of M+H
                # used everywhere in this workflow. The feature's adduct was already changed in getUnknownsGeneral().
                "CN=[N+](C=O)C1=CC=CC=C1",
                "CCC(=O)NC1=CC=CC=C1", # in PubChem, but comp annotation rank >25
                "O=C1CONN1c2ccccc2",
                "C1=CC=C(C=C1)NC(=O)CCO", # is in PubChem, but no we lack MS/MS
                NA
            ),
            # manually looked up PubChem identifiers
            CID = c(
                7671,
                NA,
                12107,
                NA,
                15561474,
                NA
            )
        )
    )
    # add in extra metadata for eg reporting later
    selectedCandidates <- Map(names(selectedCandidates), selectedCandidates, f = function(parn, sel)
    {
        if (nrow(sel) == 0)
            return(sel)
        parinf <- screenInfo(unkData$fGroups)[name == parn]
        ft <- formUTabs[[parinf$group]]
        ft <- merge(sel, ft[, c("group", "neutral_formula", grep("^annSim", names(ft), value = TRUE), "TPScore", "formFit"),
                            with = FALSE],
                    by.x = c("group", "formula"), by.y = c("group", "neutral_formula"), sort = FALSE, all.x = TRUE)
        # HACK: for M149_R515_2917 we adjusted the formula to one that wasn't in the annotation candidates --> this results
        # in NA formFit/TP score, so fix here
        ft[group == "M149_R515_2917", formFit := calcFormulaFit(formula, parinf$formula)]
        ft[group == "M149_R515_2917", TPScore := formFit + NAToZero(annSim)]
    })
    selectedCandidates <- lapply(selectedCandidates, function(sel)
    {
        sel <- data.table::copy(sel)
        sel[, InChIKey := patRoon:::babelConvert(SMILES, "smi", "inchikey")]
        return(sel)
    })
    
    # create specialized MetFrag compound library with formula TP candidates and perform compound annotation
    MFLibUnkForm <- data.table::rbindlist(selectedCandidates, idcol = "parent", fill = TRUE)
    MFLibUnkForm <- MFLibUnkForm[!is.na(SMILES)]
    MFLibUnkForm[, Identifier := paste0(parent, "-", group)]
    data.table::setnames(MFLibUnkForm, "formula", "MolecularFormula")
    MFLibUnkForm[, InChI := patRoon:::babelConvert(SMILES, "smi", "inchi")]
    data.table::fwrite(MFLibUnkForm, "output/unk_forms-database.csv")
    compoundsUnkForm <- getCompounds(unkData$fGroups[, MFLibUnkForm$group], unkData$mslists, unkData$formulas,
                                     database = "csv",
                                     extraOpts = list(LocalDatabasePath = "output/unk_forms-database.csv"))
    
    selectedCandidates <- Map(names(selectedCandidates), selectedCandidates, f = function(parn, self)
    {
        if (nrow(self) > 0)
        {
            self[, annSimForm := mapply(group, formula, FUN = calcAnnSim,
                                        MoreArgs = list(featAnn = unkData$formulas, mslists = unkData$mslists, formulas = NULL))]
            self[, annSim := mapply(group, patRoon:::getIKBlock1(InChIKey), FUN = calcAnnSim,
                                    MoreArgs = list(featAnn = compoundsUnkForm, mslists = unkData$mslists, formulas = NULL))]
            self[, annSimBoth := mapply(group, patRoon:::getIKBlock1(InChIKey), FUN = calcAnnSim,
                                        MoreArgs = list(featAnn = compoundsUnkForm, mslists = unkData$mslists, formulas = unkData$formulas))]
            self[, formRank := mapply(group, formula, FUN = function(grp, form)
            {
                return(match(form, unkData$formulas[[grp]]$neutral_formula, nomatch = NA_integer_))
            })]
            self[, compRank := mapply(group, InChIKey, FUN = function(grp, IK)
            {
                return(match(IK, compoundsUnkForm[[grp]]$InChIKey, nomatch = NA_integer_))
            })]
            self[, estIDLevel := mapply(estIDLevels, annotations(unkData$formulas)[group], annotations(compoundsUnkForm)[group],
                                        formRank, compRank, annSimForm, annSim)]
            
            if (!all(is.na(self$SMILES)))
            {
                self[!is.na(SMILES), c("maxTPSimAllPred", "maxTPSimAllPredSMILES") :=
                         data.table::rbindlist(future.apply::future_lapply(SMILES, calcMaxSim, as.data.table(SuSData$TPsCons)$SMILES))]
                parSMI <- getParentSuspList()[name == parn]$SMILES
                self[!is.na(SMILES), compFit := future.apply::future_sapply(SMILES, calcStructFitFMCS, parSMI, au = 0, bu = 0)]
                self[, TPScore_comp := calcTPScoreComp(compFit, maxTPSimAllPred, annSim), by = seq_len(nrow(self))]
            }
            else
                self[, c("maxTPSimAllPred", "maxTPSimAllPredSMILES", "compFit", "TPScore_comp") := NA_real_]
        }
        return(self)
    })
    
    return(list(
        selectedCandidates = selectedCandidates,
        compoundsUnkForm = compoundsUnkForm
    ))
}
