getUnkMetricsEval <- function(SuSData)
{
    return(data.table::rbindlist(Map(names(SuSData$TPsCons), parents(SuSData$TPsCons)$SMILES, parents(SuSData$TPsCons)$formula, f = function(parName, parSMI, parForm)
    {
        ret <- SuSData$TPsCons[[parName]][, .(name, SMILES, formula)]
        ret <- ret[SMILES != "[NH4]"] # ChemmineR seems to have issues with this one
        
        # fits (+)
        ret[, compFitPos := future.apply::future_sapply(SMILES, calcStructFitFMCS, parSMI, au = 0, bu = 0)]
        ret[, formFitPos := sapply(formula, calcFormulaFit, parForm)]
        
        # max sims
        ret[, maxTPSimPos := future.apply::future_sapply(SMILES, function(smi) calcMaxSim(smi, setdiff(ret$SMILES, smi))$maxSim)]
        otherTPSMILES <- as.data.table(SuSData$TPsCons)[parent != parName]$SMILES
        ret[, maxTPSimNeg := future.apply::future_sapply(SMILES, function(smi) calcMaxSim(smi, setdiff(otherTPSMILES, smi))$maxSim)]
        
        # sims for (-)
        otherParSMILES <- setdiff(parents(SuSData$TPsCons)$SMILES, parSMI)
        otherParForms <- setdiff(parents(SuSData$TPsCons)$formula, parForm)
        otherFits <- data.table::rbindlist(Map(ret$SMILES, ret$formula, f = function(smi, form)
        {
            data.table::data.table(
                SMILES = smi,
                compFitNeg = sapply(otherParSMILES, calcStructFitFMCS, smi, au = 0, bu = 0),
                formFitNeg = sapply(otherParForms, calcFormulaFit, form)
            )
        }))
        
        ret <- merge(ret, otherFits, by = "SMILES", sort = FALSE)
        
        return(ret[])
    }), idcol = "parent"))
}
