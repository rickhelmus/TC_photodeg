getSuspectsFromStructures <- function(featData, bgMSMS)
{
    TPsCTS <- generateTPs("cts", parents = getParentSuspList(), transLibrary = "combined_photolysis_abiotic_hydrolysis",
                          generations = 4, calcSims = FALSE)
    TPsBTE <- generateTPs("biotransformer", parents = getParentSuspList(), generations = 4)
    withOpt(MP.maxProcs = 4, {
        TPsBTH <- generateTPs("biotransformer", type = "allHuman", parents = getParentSuspList(), generations = 2, MP = TRUE)
    })
    TPsPC <- generateTPs("library", parents = getParentSuspList(), generations = 4)
    TPsLIT <- generateTPs("library", parents = getParentSuspList(),
                                TPLibrary = data.table::fread("misc/lib-lit-expanded.csv"))
    TPsCons <- consensus(TPsCTS, TPsBTE, TPsBTH, TPsPC, TPsLIT, labels = c("CTS", "BTE", "BTH", "PC", "LIT"))
    
    # amend suspect screening results with TPs.
    # NOTE: we add the parent to the TP suspect manually (ie includeParents = FALSE) to retain their retention times
    suspListTPs <- convertToSuspects(TPsCons, includeParents = FALSE)
    data.table::fwrite(suspListTPs, "output/suspects_TPs-structures.csv")
    suspListTPs <- rbind(getParentSuspList(), suspListTPs, fill = TRUE)
    
    fGroups <- doScreenSuspects(featData$fGroupsSel, suspListTPs, onlyHits = TRUE)
    
    # formula/compound annotation for structure suspects
    mslists <- getMSPL(fGroups, bgMSMS)
    formulas <- getFormulas(fGroups, mslists)
    formulas <- filterFormulasForSusps(fGroups, formulas)
    # make customized compound database for MetFrag, so we don't rely on the presence of the TP suspects in public databases
    # such as PubChem.
    convertToMFDB(TPsCons, "output/TP-database.csv", includeParents = TRUE)
    compounds <- getCompounds(fGroups, mslists, formulas, database = "csv",
                              extraOpts = list(LocalDatabasePath = "output/TP-database.csv"))
    
    # omit mass peaks without clear annotations
    mslists <- filterMSPLForAnn(mslists, fGroups, compounds)
    
    fGroups <- annotateSuspects(fGroups, formulas = formulas, compounds = compounds, MSPeakLists = mslists,
                                IDFile = "misc/IDLevelRules.yml")
    
    componentsTPs <- generateComponents(fGroups[, suspects = getParentSuspList()$name], "tp",
                                        fGroupsTPs = fGroups, TPs = TPsCons,
                                        MSPeakLists = mslists, formulas = formulas,
                                        compounds = compounds, ignoreParents = TRUE)
    componentsTPsF <- filterComponentsForParents(fGroups, componentsTPs, TRUE)
    
    fGroups <- fGroups[results = componentsTPsF]
    
    # manually checked TP candidates that should be removed (structure same as parent, in-source fragment etc)
    removeTPs <- data.table::data.table(
        parent = c(
            "metoprolol",
            "sulfamethoxazole",
            "sulfamethoxazole",
            "phenazone",
            "phenazone",
            "sulfamethoxazole"
        ),
        group = c(
            "M268_R384_4267", # same as parent
            "M254_R349_3557", # same as parent
            "M94_R530_4020", # in-source fragment (of M150)
            "M94_R530_4020", # idem
            "M117_R529_2999", # in-source fragment (of M150)
            "M124_R144_3393" # most likely can't ionize with LC-MS
        ),
        SMILES = c(
            "CC(C)NCC(COC1=CC=C(C=C1)CCOC)O",
            "CC1=CC(=NO1)NS(=O)(=O)C2=CC=C(C=C2)N",
            "Nc1ccccc1",
            "Nc1ccccc1",
            "C(=O)(C(=[OH])CC=C)[OH2]",
            "O=N(=O)c1ccccc1"
        )
    )
    for (i in seq_len(nrow(removeTPs)))
    {
        rmTP <- removeTPs[i]
        componentsTPsF <- delete(componentsTPsF, i = componentInfo(componentsTPsF)[parent_name == rmTP$parent]$name, j = function(ct, ...)
        {
            return(ct$group == rmTP$group & ct$SMILES == rmTP$SMILES)
        })
    }
    
    fGroups <- fGroups[results = componentsTPsF]
    
    # HACK: remove suspects hits from the removeTPs list above and those that were post-filtered from the components
    structSuspKeep <- function(name, grp, SMI)
    {
        cTab <- as.data.table(componentsTPsF)
        return(name %in% getParentSuspList()$name |
                   (grp %in% cTab$group & SMI %in% cTab$SMILES))
    }
    fGroups@screenInfo <- fGroups@screenInfo[structSuspKeep(name, group, SMILES)]
    
    return(list(
        TPsCTS = TPsCTS, TPsBTE = TPsBTE, TPsBTH = TPsBTH, TPsPC = TPsPC, TPsLIT = TPsLIT, TPsCons = TPsCons,
        fGroups = fGroups,
        mslists = mslists,
        formulas = formulas,
        compounds = compounds,
        componentsTPs = componentsTPsF
    ))
}
