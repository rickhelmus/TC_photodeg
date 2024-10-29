getSuspectsFromFormulas <- function(featData, SuSData, bgMSMS)
{
    # get default metabolic logic rules of patRoon
    trans <- patRoon:::TPsLogicTransformations
    # and add addition of H2O2
    trans <- rbind(trans, data.frame(transformation = "hydroperoxidation", add = "H2O2", sub = "", retDir = 0))
    genFormLib <- genFormulaTPLibrary(getParentSuspList(), transformations = trans, generations = 2)
    
    # add formula from literature search: https://doi.org/10.1016/j.scitotenv.2016.03.057
    genFormLib <- rbind(genFormLib, data.table::data.table(
        parent_name = "sulfamethoxazole", parent_formula = "C10H11N3O3S",
        parent_neutralMass = genFormLib[parent_name == "sulfamethoxazole"]$parent_neutralMass[1],
        TP_name = "sulfa-PoirierLarabie-1", TP_formula = "C10H10N2O4S",
        TP_neutralMass = rcdk::get.formula("C10H10N2O4S")@mass,
        generation = 1, retDir = 0
    ), fill = TRUE)
    TPsForm <- generateTPs("library_formula", parents = getParentSuspList(), TPLibrary = genFormLib, generations = 2,
                           matchGenerationsBy = "formula")
    
    # Omit formulas already present as structure
    TPsFormF <- delete(TPsForm, j = function(tab, par, ...) tab$formula %in% SuSData$TPsCons[[par]]$formula)
    
    # amend suspect screening results with TPs.
    # NOTE: we add the parent to the TP suspect manually (ie includeParents = FALSE) to retain their retention times
    suspListTPs <- convertToSuspects(TPsFormF, includeParents = FALSE)
    data.table::fwrite(suspListTPs, "output/suspects_TPs-formulas.csv")
    suspListTPs <- rbind(getParentSuspList(), suspListTPs, fill = TRUE)
    
    fGroups <- doScreenSuspects(featData$fGroupsSel, suspListTPs, onlyHits = TRUE)
    
    # the next steps largely follow that of the structure suspects annotation
    mslists <- getMSPL(fGroups, bgMSMS)
    formulas <- getFormulas(fGroups, mslists)
    formulas <- filterFormulasForSusps(fGroups, formulas)
    mslists <- filterMSPLForAnn(mslists, fGroups)
    fGroups <- annotateSuspects(fGroups, formulas = formulas, MSPeakLists = mslists, IDFile = "misc/IDLevelRules.yml")
    
    componentsTPs <- generateComponents(fGroups[, suspects = getParentSuspList()$name], "tp",
                                        fGroupsTPs = fGroups, TPs = TPsFormF,
                                        MSPeakLists = mslists, formulas = formulas, ignoreParents = TRUE)
    componentsTPsFormF <- filterComponentsForParents(fGroups, componentsTPs, TRUE)
    fGroups <- fGroups[results = componentsTPsFormF]
    
    return(list(
        TPsFormF = TPsFormF,
        fGroups = fGroups,
        mslists = mslists,
        formulas = formulas,
        componentsTPs = componentsTPs
    ))
}
