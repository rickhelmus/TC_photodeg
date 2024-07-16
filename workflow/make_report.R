makeReport <- function(featData, SuSData, SuFData, unkData, unkFormData, unkCompData, confData)
{
    featureRemarks <- data.table::data.table(
        type = c(
            "structure",
            "structure",
            "structure",
            "structure",
            "structure",
            "structure",
            "structure",
            "formula",
            "unknownForm",
            "unknownForm"
        ),
        parent = c(
            "flecainide",
            "metoprolol",
            "metoprolol",
            "metoprolol",
            "sulfamethoxazole",
            "sulfamethoxazole",
            "sulfamethoxazole",
            "metoprolol",
            "sulfamethoxazole",
            "phenazone"
        ),
        group = c(
            "M147_R255_4120",
            "M284_R326_2405",
            "M268_R384_4267",
            "M74_R251_2797",
            "M174_R162_5656",
            "M272_R322_3885",
            "M110_R126_3769",
            "M302_R337_5079",
            "M190_R159_4149",
            "M122_R471_3775"
        ),
        candidate = c(
            "C(CC(=O)O)CC(N)CN",
            "CC(C)NCC(COC1=CC=C(C=C1)C(COC)O)O",
            "CC(C)NC(CCOC1=CC=C(C=C1)CCOC)O",
            "C=C(O)CN",
            "Nc1ccc(cc1)S(O)(=O)=O",
            "CC(=O)CNC(=O)NS(=O)(=O)c1ccc(N)cc1",
            "C1(=CC=C(C=C1)N)O",
            "C15H27NO5",
            "C1=CC(=CC=C1O)S(=O)(=O)NO",
            "C1=CC=C(C=C1)NC=O"
        ),
        condition = list(
            "U",
            c("UH", "UH"),
            "UHN",
            c("UH", "U"),
            "U",
            "all",
            "UHN",
            "UH",
            "U",
            "UH"
        ),
        samples = list(
            "all",
            c("single", "dark"),
            "mix_150",
            c("mix_150", "all"),
            "mix_75",
            "all",
            "mix_150",
            "mix_150",
            "single",
            "all"
        ),
        remarks = list(
            6,
            c(1, 3),
            3,
            c(2, 6),
            4,
            7,
            3,
            5,
            8,
            4
        )
    )
    
    # NOTE: get parents from unfiltered fGroups: blank subtraction removes them from experiments
    fgPars <- screenSuspects(featData$fGroupsRaw, getParentSuspList(), onlyHits = TRUE, adduct = "[M+H]+")
    
    makePDReport(fgPars, SuSData$fGroups, SuSData$mslists, SuSData$formulas, SuSData$compounds, SuSData$componentsTPs,
                 SuFData$fGroups, SuFData$mslists, SuFData$formulas, SuFData$componentsTPs, unkData$fGroups,
                 unkData$mslists, unkData$formulas, unkCompData$compoundsP, unkFormData$compoundsUnkForm,
                 unkData$componentsTPs, SuSData$TPsPC, SuSData$TPsLIT, unkFormData$selectedCandidates,
                 unkCompData$selectedCandidates, featureRemarks, confData$confirmedIDLevels, "output/repPD",
                 clear = TRUE)
}
