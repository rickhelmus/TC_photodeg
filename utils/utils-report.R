plotFMCSDepict <- function(SMILES1, SMILES2, title1 = NULL, title2 = NULL, invert = FALSE, rotate = 0, flip = FALSE, split = TRUE,
                           out = NULL, ...)
{
    fmcs <- getFMCS(SMILES1, SMILES2, ...)
    
    # HACK: somehow depict seems to conflict with devtools::load_all(patRoon)? --> run with callr
    makePlot <- function(SMILES, hlAtoms, invert, rotate, flip, first, out = NULL, SMILES1 = NULL)
    {
        # based on discussions from https://github.com/CDK-R/depict/issues/5
        
        library(depict)
        
        color <- J("java.awt.Color")
        
        mol <- parse_smiles(SMILES)
        
        if (rotate != 0)
        {
            # From CDK depict: https://github.com/cdk/depict/blob/82a444544734567b459543df11321264465fb2b6/cdkdepict-lib/src/main/java/org/openscience/cdk/app/DepictController.java#L456
            GU <- J("org.openscience.cdk.geometry.GeometryUtil")
            GU$rotate(mol, GU$get2DCenter(mol), rotate)
        }
        if (flip)
        {
            # From CDK depict: https://github.com/cdk/depict/blob/82a444544734567b459543df11321264465fb2b6/cdkdepict-lib/src/main/java/org/openscience/cdk/app/DepictController.java#L461
            v2d <- J("javax.vecmath.Point2d")
            for (i in seq_len(mol$getAtomCount()))
            {
                a <- mol$getAtom(i - 1L)
                cur2d <- a$getPoint2d()
                a$setPoint2d(new(v2d, -cur2d$x, cur2d$y))
            }
        }
        
        highlight <- new(J("java/util/HashSet"))
        
        if (invert)
            hlAtoms <- setdiff(seq_len(mol$getAtomCount()), hlAtoms)
        
        for (i in hlAtoms)
        {
            a <- mol$getAtom(i - 1L)
            highlight$add(a)
        }
        for (i in seq_len(mol$getBondCount()))
        {
            b <- mol$getBond(i - 1L)
            if ((b$getBegin()$getIndex() + 1) %in% hlAtoms && (b$getEnd()$getIndex() + 1) %in% hlAtoms)
                highlight$add(b)
        }
        
        dep <- depiction() |>
            highlight_atoms(highlight, if (first) color$lightGray else color$orange) |>
            set_zoom(3) |>
            color_atoms()
        if (!first)
            dep <- outerglow(dep)
        md <- depict(dep, mol)
        if (is.null(out))
            return(get_image(md))
        md$writeTo("svg", out)
    }
    
    # HACK: grid.raster() doesn't work well with split.screen, so use base plot()
    if (split)
    {
        scr <- split.screen(c(2, 1))
        screen(scr[1])
        withr::with_par(list(mar = rep(0, 4)), {
            callr::r(makePlot, args = list(SMILES1, fmcsR::mcs1(fmcs)[[2]][[1]], invert, 0, FALSE, TRUE)) |>
                magick::image_read() |>
                plot()
        })
        screen(scr[2])
        withr::with_par(list(mar = rep(0, 4)), {
            callr::r(makePlot, args = list(SMILES2, fmcsR::mcs2(fmcs)[[2]][[1]], invert, rotate, flip,
                                           FALSE, SMILES1 = SMILES1)) |>
                magick::image_read() |>
                plot()
        })
        close.screen(scr)
    }
    else
    {
        pl1 <- callr::r(makePlot, args = list(SMILES1, fmcsR::mcs1(fmcs)[[2]][[1]], invert, rotate, flip,
                                              TRUE, out[1]))
        pl2 <- callr::r(makePlot, args = list(SMILES2, fmcsR::mcs2(fmcs)[[2]][[1]], invert, rotate, flip,
                                              FALSE, out[2], SMILES1))
        if (is.null(out))
        {
            grid::grid.raster(pl1)
            grid::grid.newpage()
            grid::grid.raster(pl2)
        }
    }
}

plotStruct <- function(SMILES, out)
{
    # simplified version of above, mainly to plot parents in consistent way
    
    # HACK: somehow depict seems to conflict with devtools::load_all(patRoon)? --> run with callr
    makePlot <- function(SMILES, out)
    {
        library(depict)
        
        mol <- parse_smiles(SMILES)
        dep <- depiction() |>
            set_zoom(3) |>
            depict(mol)
        dep$writeTo("svg", out)
    }
    
    callr::r(makePlot, args = list(SMILES, out))
}

plotReg <- function(fGroups, minRSQ = 0.5, ...)
{
    regTabs <- getRegTabs(fGroups, minRSQ = minRSQ)
    regTabs <- regTabs[sapply(regTabs, nrow) > 0]
    colors <- RColorBrewer::brewer.pal(3, "Dark2")
    names(colors) <- names(regTabs)
    anas <- analyses(fGroups)
    anaInfo <- analysisInfo(fGroups)
    
    plot(0, type = "n", xlab = expression(paste("Initial parent concentration (", mu, "g/l)")),
         ylab = "Normalized intensity (%)", xlim = c(0, 150), ylim = c(0, 100), bty = "l", ...)
    
    for (rn in names(regTabs))
    {
        anasRN <- intersect(anas, names(regTabs[[rn]]))
        ints <- unlist(regTabs[[rn]][, anasRN, with = FALSE])
        ints[ints == 0] <- NA
        ints <- ints * 100 # convert to %
        concs <- anaInfo[match(anasRN, anaInfo$analysis), "conc"]
        
        points(concs, ints, pch = 20, col = colors[[rn]])
        
        usr <- par("usr")
        # from https://stackoverflow.com/a/10046370
        clip(min(concs), max(concs), min(ints, na.rm = TRUE), max(ints, na.rm = TRUE))
        abline(lm(ints ~ concs), col = colors[[rn]])
        do.call("clip", as.list(usr))  # reset to plot region (from ?clip examples)
    }
    
    leg <- longToShortCondName(names(regTabs))
    
    leg <- mapply(leg, regTabs, FUN = function(rn, rt) as.expression(bquote(mix[.(rn)] ~ "(" ~ R^2 ~ .(round(rt$RSQ, 2)) ~ ")")))
    
    legend("topleft", legend = leg, col = colors, text.col = colors, bty = "n", lty = 1)
}

cropSVG <- function(file)
{
    # https://stackoverflow.com/a/64116358
    # NOTE: this seems to stall when InkScape is open ...
    processx::run("C:\\Program Files\\Inkscape\\bin\\inkscape.exe",
                  c("--export-plain-svg", sprintf("--export-filename=%s", file), "--export-area-drawing", file),
                  echo = TRUE)
}

makeSVGPlot <- function(out, code, ..., n = 1)
{
    outAll <- if (n > 1) sprintf(out, seq_len(n)) else out
    if (!all(file.exists(outAll)))
    {
        svglite::svglite(out, ..., bg = "transparent")
        # svg(out, bg = "transparent")
        tryCatch({ par(bg = "transparent", family = "sans"); force(code) }, error = function(...) NULL, finally = dev.off())
        for (f in outAll)
            cropSVG(f)
    }
    return(outAll)
}

makePDReport <- function(fGroupsPar, fGroupsStruct, mslistsStruct, formulasStruct, compoundsStruct, componentsStruct,
                         fGroupsForm, mslistsForm, formulasForm, componentsForm, fGroupsUnk, mslistsUnk, formulasUnk,
                         compoundsUnk, compoundsUnkForm, componentsUnk, TPsLibPC, TPsLibManual, selectedUnknownFormulas,
                         selectedUnknownStructs, featureRemarks, confirmedIDLevels, out, clear = FALSE, onlyTypes = NULL)
{
    if (clear && dir.exists(out))
        unlink(out, recursive = TRUE)
    
    patRoon:::mkdirp(out)
    out <- normalizePath(out)
    
    plotOut <- file.path(out, "plots")
    patRoon:::mkdirp(plotOut)
    
    allFGroups <- list(structure = fGroupsStruct, formula = fGroupsForm, unknownStruct = fGroupsUnk,
                       unknownForm = fGroupsUnk, parent = fGroupsPar)
    
    featsNormList <- sapply(allFGroups, function(fg)
    {
        fg <- fg[, rGroups = grep("^UV_", replicateGroups(fg), value = TRUE)]
        fg <- normInts(fg, groupNorm = TRUE)
        return(as.data.table(fg, normalized = TRUE, average = TRUE))
    }, simplify = FALSE)
    
    regTabsList <- Map(allFGroups, names(allFGroups), f = function(fg, t) getRegTabs(fg, 0, filter = FALSE))
    regTabsListF <- Map(allFGroups, names(allFGroups), f = function(fg, t) getRegTabs(fg, 0.5, filter = TRUE))
    
    # markup remarks
    featureRemarks <- data.table::copy(featureRemarks)
    featureRemarks[, remarks := lapply(remarks, sapply,
                                       function(x) sprintf("<sup><a href=\"#remark-%d\">%d</a></sup>", x, x))]
    
    getTPAlgosStruct <- function(mergedBy, IK, par)
    {
        mergedBy <- unlist(strsplit(mergedBy, ","))
        mergedBy <- sapply(mergedBy, function(mb)
        {
            if (mb == "PC")
                return(sprintf("%s (%s)", mb, patRoon:::makeDBIdentLink("pubchem", unique(TPsLibPC[[par]][InChIKey == IK]$CID))))
            if (mb == "LIT")
            {
                doi <- TPsLibManual[[par]][InChIKey == IK]$doi
                return(sprintf("%s (%s)", mb, paste0(sprintf("<a href=\"%s\">%s</a>", doi, letters[seq_along(doi)]),
                                                     collapse = ", ")))                
            }
            return(mb)
        })
        return(paste0(mergedBy, collapse = ", "))
    }
    getTPAlgosUnkStruct <- function(IK, group, compounds)
    {
        ct <- compounds[[group]][InChIKey == IK]
        patRoon:::makeDBIdentLink("pubchem", ct$identifier)
    }
    getTPAlgosUnkForm <- function(IK, parn)
    {
        if (is.na(IK))
            return("\\-")
        cid <- selectedUnknownFormulas[[parn]][InChIKey == IK]$CID
        return(if (is.null(cid) || is.na(cid)) "\\-" else patRoon:::makeDBIdentLink("pubchem", cid))
    }
    
    getMixInt <- function(grp, conc, exp, type)
    {
        val <- featsNormList[[type]][group == grp, paste0(exp, "_", conc), with = FALSE]
        if (val == 0)
            return("")
        return(paste0(round(val * 100, 0), "%"))
    }
    getSingleInt <- function(grp, exp, par, type)
    {
        pat <- paste0("^", exp, "_", strtrim(par, 4))
        v <- featsNormList[[type]][group == grp, grep(pat, names(featsNormList[[type]]), value = TRUE), with = FALSE]
        if (v == 0)
            return("")
        return(paste0(round(v * 100, 0), "%"))
    }
    getDarkInt <- function(grp, exp, type)
    {
        val <- featsNormList[[type]][group == grp, paste0(exp, "_150_dark"), with = FALSE]
        if (val == 0)
            return("")
        return(paste0(round(val * 100, 0), "%"))
    }
    getRSQ <- function(grp, exp, type)
    {
        v <- regTabsList[[type]][[exp]][group == grp]$RSQ
        if (length(v) == 0)
            return("")
        return(sprintf("%.2f", v))
    }
    getP <- function(grp, exp, type)
    {
        v <- regTabsList[[type]][[exp]][group == grp]$p
        if (length(v) == 0)
            return("")
        return(sprintf("%.3f", v))
    }
    getSlope <- function(grp, exp, type)
    {
        v <- regTabsList[[type]][[exp]][group == grp]$slope
        if (length(v) == 0)
            return("")
        return(paste0(round(v * 100, 2), "%"))
    } 
    getExpDT <- function(grp, parName, type)
    {
        ret <- data.table::data.table(
            Condition = c(
                "U",
                "UH",
                "UHN"
            ),
            mix_25 = c(
                getMixInt(grp, 25, "UV", type),
                getMixInt(grp, 25, "UV_H2O2", type),
                getMixInt(grp, 25, "UV_H2O2_NOM", type)
            ),
            mix_75 = c(
                getMixInt(grp, 75, "UV", type),
                getMixInt(grp, 75, "UV_H2O2", type),
                getMixInt(grp, 75, "UV_H2O2_NOM", type)
            ),
            mix_150 = c(
                getMixInt(grp, 150, "UV", type),
                getMixInt(grp, 150, "UV_H2O2", type),
                getMixInt(grp, 150, "UV_H2O2_NOM", type)
            ),
            RSQ = c(
                getRSQ(grp, "UV", type),
                getRSQ(grp, "UV_H2O2", type),
                getRSQ(grp, "UV_H2O2_NOM", type)
            ),
            p = c(
                getP(grp, "UV", type),
                getP(grp, "UV_H2O2", type),
                getP(grp, "UV_H2O2_NOM", type)
            ),
            slope = c(
                getSlope(grp, "UV", type),
                getSlope(grp, "UV_H2O2", type),
                getSlope(grp, "UV_H2O2_NOM", type)
            ),
            single = c(
                getSingleInt(grp, "UV", parName, type),
                getSingleInt(grp, "UV_H2O2", parName, type),
                getSingleInt(grp, "UV_H2O2_NOM", parName, type)
            ),
            dark = c(
                getDarkInt(grp, "UV", type),
                getDarkInt(grp, "UV_H2O2", type),
                getDarkInt(grp, "UV_H2O2_NOM", type)
            )
        )
        if (grp %in% featureRemarks$group)
        {
            fr <- featureRemarks[group == grp]
            if (fr$condition[[1]][1] != "all")
            {
                # add in remarks
                rows <- data.table::fcase(fr$condition[[1]] == "U", 1L,
                                          fr$condition[[1]] == "UH", 2L,
                                          fr$condition[[1]] == "UHN", 3L)
                cols <- ifelse(fr$sample[[1]] == "all", "Condition", fr$sample[[1]])
                for (i in seq_along(rows))
                    data.table::set(ret, rows[i], j = cols[i], value = paste0(ret[rows[i], cols[i], with = FALSE], fr$remarks[[1]][i]))
            }
        }
        return(ret)
    }
    
    doChromPlot <- function(grp, fg)
    {
        makeSVGPlot(file.path(plotOut, sprintf("chrom-%s.svg", grp)), {
            withr::with_par(list(mar = c(4, 4, 1.5, 0.2)), {
                plotChroms(fg[, grp], showFGroupRect = FALSE, showPeakArea = FALSE, retMin = TRUE,
                           bty = "l", EICParams = getDefEICParams(topMost = 1), title = "",
                           xlim = (groupInfo(fg)[grp, "rts"] / 60) + c(-0.75, 0.75), intMax = "feature")
            })
        }, pointsize = 16, width = 8, height = 6)
    }
    doRegrPlot <- function(grp, fg, minRSQ = 0.5)
    {
        makeSVGPlot(file.path(plotOut, sprintf("regr-%s.svg", grp)), {
            withr::with_par(list(mar = c(4, 4, 1.5, 0.5)), {
                plotReg(fg[, grp], minRSQ)
            })
        }, pointsize = 18, width = 8, height = 5)
    }
    
    reportCmp <- function(fGroups, mslists, formulas, compounds, components, type)
    {
        if (!is.null(onlyTypes) && !type %in% onlyTypes)
            return(list())
        
        cInfo <- componentInfo(components)
        if (grepl("^unknown", type))
        {
            # HACK: 'fix' parent name: patRoon makes it group names when componentization is used w/out TPs
            cInfo <- data.table::copy(cInfo)
            cInfo[, parent_name := screenInfo(fGroups)[match(parent_name, group)][!is.na(rt)]$name]
        }
        
        scrInfo <- list()
        scrInfo$header <- switch(type,
                                 structure = "Candidates from structure suspect screening",
                                 formula = "Candidates from formula suspect screening",
                                 unknownStruct = "Candidates for unknowns from compound annotations",
                                 unknownForm = "Candidates for unknowns from formula annotations")
        
        scrInfo$parents <- lapply(split(cInfo, seq_len(nrow(cInfo))), function(parRow)
        {
            cTab <- data.table::copy(components[[parRow$name]])
            
            tabFeatsSusp <- as.data.table(fGroups, collapseSuspects = NULL)
            # NOTE: assume there aren't any multiple parent hits, and make sure parent and not TP isomer is selected
            parSuspRow <- tabFeatsSusp[group == parRow$parent_group & !is.na(susp_rt)]
            
            # HACK: for unknowns: subset and add missing info from selections
            if (type == "unknownForm")
            {
                cTab <- merge(cTab, selectedUnknownFormulas[[parSuspRow$susp_name]], by = "group")
                if (nrow(cTab) > 0)
                {
                    cTab[, neutralMass := rcdk::get.formula(formula)@mass, by = seq_len(nrow(cTab))]
                }
            }
            else if (type == "unknownStruct")
            {
                cTab <- merge(cTab, selectedUnknownStructs[[parSuspRow$susp_name]], by = "group")
                if (nrow(cTab) > 0)
                {
                    cTab[, c("neutralMass", "formula", "SMILES") := {
                        IK <- InChIKey
                        compounds[[group]][IK == InChIKey, .(neutralMass, neutral_formula, SMILES)]
                    }, by = seq_len(nrow(cTab))]
                }
            }
            
            if (nrow(cTab) == 0)
                return(list())
            
            parentInfo <- list()
            
            parentInfo$header <- sprintf("Parent '%s'", parRow$parent_name)
            
            if (F) {
                parentInfo$struct <- file.path(plotOut, sprintf("struct_par-%s.svg", parRow$parent_name))
                plotStruct(parSuspRow$susp_SMILES, parentInfo$struct)
                parentInfo$chrom <- doChromPlot(parSuspRow$group, fGroupsPar)
                parentInfo$regr <- doRegrPlot(parSuspRow$group, fGroupsPar, 0)
                parentInfo$chromHeader <- sprintf("RT: %.1f min; m/z: %.4f", 
                                                  parRow$parent_ret / 60, parRow$parent_mz)
                parentInfo$expTab <- getExpDT(parSuspRow$group, parSuspRow$susp_name, "parent")
            }
            
            # for TP naming
            TPNMs <- if (!grepl("^unknown", type))
            {
                mapply(cTab$group, cTab$TP_name,
                       FUN = function(g, n) tabFeatsSusp[group == g & susp_name == n]$neutralMass)
            }
            else
                cTab$neutralMass
            
            # assign ID levels in component tables so they can be used for sorting and reporting
            if (is.null(cTab[["estIDLevel"]])) # already set for unknowns
            {
                cTab[, estIDLevel := {
                    grp <- group
                    tabFeatsSusp[group == grp & susp_name == TP_name]$susp_estIDLevel
                }, by = c("group", "TP_name")]
            }
            
            if (!is.null(cTab[["SMILES"]]))
            {
                cTab <- merge(cTab, confirmedIDLevels[parent == parRow$parent_name, -"InChIKey"], by = c("group", "SMILES"),
                              all.x = TRUE, sort = FALSE)
                cTab[is.na(IDL), IDL := estIDLevel]
            }
            else
                cTab[, IDL := estIDLevel]
            
            cTab[, reportName := makeTPNames(type, parRow$parent_name, cTab$mz, if (!is.null(cTab[["InChIKey"]])) cTab$InChIKey else cTab$formula)]
            
            data.table::setorderv(cTab, intersect(c("mz", "IDL"), names(cTab)))
            
            cTabFGSplit <- split(cTab, by = "group")
            
            parentInfo$features <- sapply(unique(cTab$group), function(TPGrp)
            {
                featInfo <- list()
                
                featInfo$chrom <- doChromPlot(TPGrp, fGroups)
                featInfo$regr <- doRegrPlot(TPGrp, fGroups)
                
                featInfo$remarks <- if (TPGrp %in% featureRemarks$group && featureRemarks[group == TPGrp]$condition[[1]][1] == "all")
                    featureRemarks[group == TPGrp]$remark[[1]]
                else
                    ""
                
                featInfo$header <- sprintf("Feature '%s'", TPGrp)
                
                featCRow <- cTabFGSplit[[TPGrp]][1] # should be all the same for FG info, so just take first
                featInfo$chromHeader <- sprintf("RT: %.1f (\U0394 %+.1f) min; m/z: %.4f (\U0394 %+.4f)",
                                                featCRow$ret / 60, featCRow$retDiff / 60, featCRow$mz, featCRow$mzDiff)
                
                featInfo$expTab <- getExpDT(TPGrp, parSuspRow$susp_name, type)
                # non-significant rows, will be greyed out
                featInfo$expTabRowsNSig <- which(sapply(regTabsListF[[type]], function(x) !TPGrp %in% x$group))
                
                featInfo$TPs <- Map(split(cTabFGSplit[[TPGrp]], seq_len(nrow(cTabFGSplit[[TPGrp]]))), seq_len(nrow(cTabFGSplit[[TPGrp]])), f = function(cRow, cmpInd)
                {
                    cat(sprintf("Doing %s-%s-%d (%s)\n", parRow$name, TPGrp, cmpInd, type))
                    
                    hasStruct <- !is.null(cRow[["SMILES"]]) && !is.na(cRow$SMILES)
                    
                    TPInfo <- list()
                    
                    TPInfo$arrow <- normalizePath("report_pd/arrow-right.svg")
                    
                    if (hasStruct)
                    {
                        TPInfo$structs <- file.path(plotOut, sprintf("structs-%s-%s-%s-%d-%d.svg", type, parRow$name,
                                                                     TPGrp, cmpInd, 1:2))
                        if (!all(file.exists(TPInfo$structs)))
                            plotFMCSDepict(parSuspRow$susp_SMILES, cRow$SMILES, NULL, NULL, invert = TRUE, flip = TRUE,
                                           au = 0, bu = 0, split = FALSE, out = TPInfo$structs)
                        
                    }
                    else
                    {
                        TPInfo$structs <- NULL
                        TPInfo$formulas <- patRoon:::subscriptFormulaHTML(c(parSuspRow$susp_formula, cRow$formula))
                    }
                    
                    TPInfo$spec <- makeSVGPlot(file.path(plotOut, sprintf("spec-%s-%s-%s-%d.svg", type, parRow$name, TPGrp, cmpInd)), {
                        withr::with_par(list(mar = c(4, 4, 2, 0.2)), {
                            formInd <- compInd <- NA
                            if (!is.null(mslists[[cRow$group]][["MSMS"]]))
                            {
                                if (!is.null(formulas) && !is.null(formulas[[cRow$group]]))
                                    formInd <- match(cRow$formula, formulas[[cRow$group]]$neutral_formula)
                                if (!is.null(compounds) && !is.null(compounds[[cRow$group]]))
                                    compInd <- match(cRow$InChIKey, compounds[[cRow$group]]$InChIKey)
                            }
                            
                            # this could be done to disable scaling, but this makes the annotations less readable...
                            # xlim <- if (!is.null(mslists[[cRow$group]]) && !is.null(mslists[[cRow$group]][["MSMS"]])) 
                            #     c(0, max(mslists[[cRow$group]][["MSMS"]]$mz) * 1.1)
                            xlim <- NULL
                            
                            if (!is.na(compInd))
                                plotSpectrum(compounds, index = compInd, groupName = cRow$group, MSPeakLists = mslists,
                                             formulas = formulas, xlim = xlim, title = "")
                            else if (!is.na(formInd))
                                plotSpectrum(formulas, index = formInd, groupName = cRow$group,
                                             MSPeakLists = mslists, xlim = xlim, title = "")
                            else
                                patRoon:::textPlot("No MS/MS annotations")
                        })
                    }, pointsize = 18, width = 7, height = 6)
                    
                    formDiff <- getFormDiff(parSuspRow$susp_formula, cRow$formula)
                    form <- sprintf("%s (%s)", patRoon:::subscriptFormulaHTML(cRow$formula),
                                    if (nzchar(formDiff)) paste("\U0394", patRoon:::subscriptFormulaHTML(formDiff, charges = FALSE))
                                    else "no difference")
                    
                    formFit <- calcFormulaFit(parSuspRow$susp_formula, cRow$formula)
                    compFit <- if (hasStruct)
                        calcStructFitFMCS(parSuspRow$susp_SMILES, cRow$SMILES, au = 0, bu = 0)
                    else
                        NULL
                    
                    if (hasStruct)
                    {
                        XlogPPar <- getXLogP(parSuspRow$susp_SMILES)
                        XlogPTP <- getXLogP(cRow$SMILES)
                    }
                    else
                    {
                        XlogPPar <- XlogPTP <- 0
                    }
                    
                    if (!grepl("^unknown", type))
                    {
                        suspRow <- tabFeatsSusp[group == cRow$group & susp_name == cRow$TP_name]
                        massError <- suspRow$susp_d_mz
                        annSimForm <- suspRow$susp_annSimForm
                        annSimComp <- if (hasStruct) suspRow$susp_annSimComp else NA_real_
                    }
                    else
                    {
                        massError <- annotations(fGroups)[group == TPGrp]$neutralMass - cRow$neutralMass
                        annSimForm <- cRow$annSimForm; annSimComp <- cRow$annSim
                    }
                    
                    as <- if (!is.na(annSimComp))
                        sprintf("%.2f (formula), %.2f (compound)", annSimForm, annSimComp)
                    else if (!is.na(annSimForm))
                        sprintf("%.2f (formula)", annSimForm)
                    else
                        "NA"
                    
                    fit <- sprintf("fit<sub>formula</sub>: %.2f", formFit)
                    if (!is.null(compFit))
                        fit <- paste0(fit, sprintf("<br>fit<sub>compound</sub>: %.2f", compFit))
                    
                    TPInfo$header <- sprintf("Candidate '%s'", cRow$reportName)
                    
                    IDLText <- if (cRow$IDL == "disproved") "<b>disproved by standard</b>" else cRow$IDL
                    if (cRow$IDL == "disproved")
                        TPInfo$header <- paste(TPInfo$header, "(_DISPROVED_)")
                    
                    TPInfo$TPTab <- data.table::data.table(
                        prop = c(
                            "Formula",
                            "SMILES",
                            "m/z error",
                            "XLog P",
                            "Data source(s)",
                            "<i>In silico</i> similarity",
                            "ID confidence level",
                            "Fit"
                        ),
                        value = c(
                            form,
                            # Escape unfortunate SMILES which are turned into links: [X](Y) --> [X]\(). This seems the
                            # only way to avoid triggering Markdown and its Math extensions
                            if (hasStruct) gsub("](", "]\\(", cRow$SMILES, fixed = TRUE) else "",
                            sprintf("%+.1f mDa", massError * 1000),
                            sprintf("%.1f (\U0394 %+.1f)", XlogPTP, XlogPTP - XlogPPar),
                            switch(type,
                                   structure = getTPAlgosStruct(cRow$mergedBy, cRow$InChIKey, parSuspRow$susp_name),
                                   unknownStruct = getTPAlgosUnkStruct(cRow$InChIKey, cRow$group, compounds),
                                   unknownForm = getTPAlgosUnkForm(cRow$InChIKey, parRow$parent_name),
                                   # NOTE: no LIT forms were found, otherwise code is needed for references like structs
                                   "\\-"
                            ),
                            as,
                            IDLText,
                            fit
                        )
                    )
                    if (!hasStruct)
                        TPInfo$TPTab <- TPInfo$TPTab[!prop %in% c("SMILES", "XLog P")]
                    if (type == "unknownStruct" || (type == "unknownForm" && hasStruct))
                    {
                        sc <- if (type == "unknownStruct")
                            sprintf("TP_score<sub>compound</sub>: %.2f", cRow$TPScore)
                        else
                            sprintf("TP_score<sub>formula</sub>: %.2f<br>TP_score<sub>compound</sub>: %.2f", cRow$TPScore,
                                    cRow$TPScore_comp)
                        TPInfo$TPTab <- rbind(TPInfo$TPTab, data.table::data.table(
                            prop = c(
                                "TP_sim<sub>max</sub>",
                                "TP_score"
                            ),
                            value = c(
                                sprintf("%.2f", cRow$maxTPSimAllPred),
                                sc
                            )
                        ))
                    }
                    else if (type == "unknownForm")
                    {
                        TPInfo$TPTab <- rbind(TPInfo$TPTab, data.table::data.table(
                            prop = "TP_score<sub>formula</sub>",
                            value = sprintf("%.2f", cRow$TPScore)
                        ))
                    }
                    
                    matchedFGs <- setdiff(tabFeatsSusp[susp_name == cRow$TP_name]$group, cRow$group)
                    if (length(matchedFGs) > 0)
                    {
                        TPInfo$TPTab <- rbind(TPInfo$TPTab, data.table::data.table(
                            prop = "Other matches",
                            value = paste0(matchedFGs, collapse = ", ")
                        ))
                    }
                    
                    return(TPInfo)
                })
                
                return(featInfo)
            }, simplify = FALSE)
            
            return(parentInfo)
        })
        names(scrInfo$parents) <- cInfo$parent_name
        return(scrInfo)
    }
    
    allMSPL <- list(mslistsStruct, mslistsForm, mslistsUnk, mslistsUnk)
    allForms <- list(formulasStruct, formulasForm, formulasUnk, formulasUnk)
    allCompounds <- list(compoundsStruct, NULL, compoundsUnk, compoundsUnkForm)
    allComponents <- list(componentsStruct, componentsForm, componentsUnk, componentsUnk)
    
    allFGroupsNoPar <- allFGroups[setdiff(names(allFGroups), "parent")]
    repList <- Map(allFGroupsNoPar, allMSPL, allForms, allCompounds, allComponents, names(allFGroupsNoPar),
                   f = reportCmp)
    
    reportEnv <- new.env()
    reportEnv$reportList <- repList
    
    # this doesn't really work    
    # xaringan::inf_mr(file.path("report_pd", "main.Rmd"), output_file = file.path(out, "reportPD.html"),
    #                  envir = reportEnv)
    rmarkdown::render(file.path("report_pd", "main.Rmd"), output_file = file.path(out, "reportPD.html"),
                      envir = reportEnv)
}
