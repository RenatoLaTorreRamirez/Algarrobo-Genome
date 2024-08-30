### Packages
library(topGO)
### Working dir
setwd("/path/to/working_directory/")

## Upload data
Prosopis_Expanded <- "Prosopis_Expanded_geneterms.txt"
Npal_Expanded <- "Npal_Expanded_geneterms.txt"
Npal_Only <- "Npal_Only_geneterms.txt"
geneID2GO <- readMappings(file="Allgene_nameterms.txt")
GO2geneID <- inverseList(geneID2GO)
geneNames <- names(geneID2GO)

## Prosopis expanded analysis
Expanded <- read.table(file=Prosopis_Expanded, sep = "\t")
Expanded <- as.character(Expanded[,1])
geneList <- factor(as.integer(geneNames %in% Expanded))
names(geneList) <- geneNames
# BP terms
BPGOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
resultBPFisher <- runTest(BPGOdata, algorithm = "weight01", statistic = "fisher")
# resultBPFisher # Put the number of GO terms in topNodes
allBPRes <- GenTable(BPGOdata, classicFisher = resultBPFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 2758)
write.table(allBPRes, file = "BP_Table_Expanded.txt", sep = "\t", quote = FALSE)
# MF terms
MFGOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
resultMFFisher <- runTest(MFGOdata, algorithm = "weight01", statistic = "fisher")
# resultMFFisher # Put the number of GO terms in topNodes
allMFRes <- GenTable(MFGOdata, classicFisher = resultMFFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 724)
write.table(allMFRes, file = "MF_Table_Expanded.txt", sep = "\t", quote = FALSE)
# CC terms
CCGOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
resultCCFisher <- runTest(CCGOdata, algorithm = "weight01", statistic = "fisher")
# resultCCFisher # Put the number of GO terms in topNodes
allCCRes <- GenTable(CCGOdata, classicFisher = resultCCFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 401)
write.table(allCCRes, file = "CC_Table_Expanded.txt", sep = "\t", quote = FALSE)

## Expanded significant only for Neltuma pallida
ExpandedNO <- read.table(file=Npal_Expanded, sep = "\t")
ExpandedNO <- as.character(ExpandedNO[,1])
geneListNO <- factor(as.integer(geneNamesNO %in% ExpandedNO))
names(geneListNO) <- geneNamesNO
# BP terms
BPGOdataNO <- new("topGOdata", ontology = "BP", allGenes = geneListNO, annot = annFUN.gene2GO, gene2GO = geneID2GONO, nodeSize = 10)
resultBPFisherNO <- runTest(BPGOdataNO, algorithm = "weight01", statistic = "fisher")
# resultBPFisherNO # Put the number of GO terms in topNodes
allBPResNO <- GenTable(BPGOdataNO, classicFisher = resultBPFisherNO, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 2758)
write.table(allBPResNO, file = "BP_Table_Expanded_Npal_Only.txt", sep = "\t", quote = FALSE)
# MF terms
MFGOdataNO <- new("topGOdata", ontology = "MF", allGenes = geneListNO, annot = annFUN.gene2GO, gene2GO = geneID2GONO, nodeSize = 10)
resultMFFisherNO <- runTest(MFGOdataNO, algorithm = "weight01", statistic = "fisher")
# resultMFFisherNO # Put the number of GO terms in topNodes
allMFResNO <- GenTable(MFGOdataNO, classicFisher = resultMFFisherNO, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 724)
write.table(allMFResNO, file = "MF_Table_Expanded_Npal_Only.txt", sep = "\t", quote = FALSE)
# CC terms
CCGOdataNO <- new("topGOdata", ontology = "CC", allGenes = geneListNO, annot = annFUN.gene2GO, gene2GO = geneID2GONO, nodeSize = 10)
resultCCFisherNO <- runTest(CCGOdataNO, algorithm = "weight01", statistic = "fisher")
# resultCCFisherNO # Put the number of GO terms in topNodes
allCCResNO <- GenTable(CCGOdataNO, classicFisher = resultCCFisherNO, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 401)
write.table(allCCResNO, file = "CC_Table_Expanded_Npal_Only.txt", sep = "\t", quote = FALSE)

# Present only for Npallida branch
ExpandedNPO <- read.table(file=Npal_Only, sep = "\t")
ExpandedNPO <- as.character(ExpandedNPO[,1])
geneList <- factor(as.integer(geneNames %in% ExpandedNPO))
names(geneList) <- geneNames
# BP terms
BPGOdataNPO <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
resultBPFisherNPO <- runTest(BPGOdataNPO, algorithm = "weight01", statistic = "fisher")
# resultBPFisherNPO # Put the number of GO terms in topNodes
allBPResNPO <- GenTable(BPGOdataNPO, classicFisher = resultBPFisherNPO, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 2758)
write.table(allBPResNPO, file = "BP_Table_NPOnly.txt", sep = "\t", quote = FALSE)
# MF terms
MFGOdataNPO <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
resultMFFisherNPO <- runTest(MFGOdataNPO, algorithm = "weight01", statistic = "fisher")
# resultMFFisherNPO # Put the number of GO terms in topNodes
allMFResNPO <- GenTable(MFGOdataNPO, classicFisher = resultMFFisherNPO, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 724)
write.table(allMFResNPO, file = "MF_Table_NPOnly.txt", sep = "\t", quote = FALSE)
# CC terms
CCGOdataNPO <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
resultCCFisherNPO <- runTest(CCGOdataNPO, algorithm = "weight01", statistic = "fisher")
# resultCCFisherNPO # Put the number of GO terms in topNodes
allCCResNPO <- GenTable(CCGOdataNPO, classicFisher = resultCCFisherNPO, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 401)
write.table(allCCResNPO, file = "CC_Table_NPOnly.txt", sep = "\t", quote = FALSE)
