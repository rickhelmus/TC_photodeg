printf <- function(...) cat(sprintf(...), sep = "")

# short names are used for e.g. replicate groups
shortenParentName <- function(par) strtrim(par, ifelse(par == "metoprolol", 4, 5))

# and acronyms for anything in the manuscript
getShortParentAcronyms <- function(parn)
{
    acr <- c(fleca = "FLE", meto = "MET", sulfa = "SMX", phena = "PHE")
    return(acr[parn])
}
getParentAcronyms <- function(parn) getShortParentAcronyms(shortenParentName(parn))

longToShortCondName <- function(ln)
{
    # convert conditions to shortened format for manuscript
    return(list(
        UV = "U", UV_H2O2 = "UH", "UV_H2O2_NOM" = "UHN"
    )[ln])
}

makeTPNames <- function(type, parn, Ms, UIDs)
{
    # all equal UIDs (IKs/formulas) should get the same name
    # all equal rounded masses get same name but with different suffix ID
    
    type <- switch(type,
                   structure = "SuS",
                   formula = "SuF",
                   unknownForm = "UnF",
                   unknownStruct = "UnC")
    parn <- getParentAcronyms(parn)
    
    tab <- data.table::data.table(M = Ms, MRound = as.integer(round(Ms)), UID = UIDs)
    tab[, TP_name := sprintf("%s-%s-M%d-%d", type, parn, MRound[1], match(UID, unique(UID))), by = "MRound"]
    
    return(tab$TP_name)
}

getFormDiff <- function(form1, form2)
{
    sfl <- patRoon:::splitFormulaToList(patRoon:::subtractFormula(form2, form1))
    ret <- ""
    subfl <- sfl[sfl < 0]
    if (length(subfl) > 0)
        ret <- paste0("-", patRoon:::formulaListToString(abs(subfl)))
    addfl <- sfl[sfl > 0]
    if (length(addfl) > 0)
        ret <- if (nzchar(ret)) paste0(ret, " +", patRoon:::formulaListToString(addfl)) else paste0("+", patRoon:::formulaListToString(addfl))
    return(ret)
}

maxOrNA <- function(...)
{
    x <- c(...)
    x <- x[!is.na(x)]
    return(if (length(x) == 0) NA_real_ else max(x))
}

maxOrZero <- function(...)
{
    m <- maxOrNA(...)
    return(if (is.na(m)) 0 else m)
}

NAToZero <- function(x) ifelse(is.na(x), 0, x)
InfToZero <- function(x) ifelse(is.infinite(x), 0, x)

getXLogP <- function(SMILES)
{
    mol <- rcdk::parse.smiles(SMILES)[[1]]
    rcdk::set.atom.types(mol)
    rcdk::do.aromaticity(mol)
    rcdk::convert.implicit.to.explicit(mol)
    return(rcdk::get.xlogp(mol))
}

getAllRPackages <- function()
{
    utils::capture.output(allPkgs <- renv::dependencies()$Package)
    
    # additional used that are not picked, ie used through patRoon
    allPkgs <- c(allPkgs, "mzR", "cliqueMS", "MS2Quant")
    
    tab <- data.table::data.table(
        package = allPkgs,
        version = sapply(allPkgs, function(p) as.character(packageVersion(p))),
        url = sapply(allPkgs, function(p)
        {
            u <- packageDescription(p)$URL
            if (is.null(u))
                return("")
            return(unlist(strsplit(u, ",|\\n|\\s"))[[1]])
        })
    )
    tab <- unique(tab)
    data.table::setorderv(tab, "package")
    
    # no URLs reported for these
    tab[package == "depict", url := "https://github.com/CDK-R/depict"]
    tab[package == "rcdk", url := "https://cdk-r.github.io/cdkr/"]
    tab[package == "RColorBrewer", url := "https://cran.r-project.org/web/packages/RColorBrewer"]
    tab[package == "job", url := "https://lindeloev.github.io/job"]
    
    tab <- tab[nzchar(url)] # rest is from base R
    
    return(tab)
}
