#' @include Utils.R
#'
utils::globalVariables(
	names = c('gene', 'set', 'CellType', 'features.plot', 'id', 'pct.exp'),
	package = 'RIRA',
	add = TRUE
)

#' @title PlotImmuneMarkers
#'
#' @description Generate a set of Seurat FeaturePlots for common immune cell markers
#' @param seuratObj A seurat object
#' @param reductions Vector of reduction(s) to use
#' @export
#' @importFrom dplyr %>%
#' @import Seurat
PlotImmuneMarkers <- function(seuratObj, reductions = c('tsne', 'umap')) {
	reductions <- intersect(reductions, names(seuratObj@reductions))
	if (length(reductions) == 0) {
		stop('None of the requested reductions are present!')
	}

	#ENSMMUG00000003532=CD8b
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('CD8A', 'CD8B', 'ENSMMUG00000003532', 'CD4', 'IL7R', 'CD3D', 'CD3E','CD3G'), 'CD8/CD4 Markers')

	#Eff v. Mem:
	#IL7R = CD127
	#IL2RA = CD25
	#PTPRC = CD45
	#SELL = CD62-L / CD-197
	# TNFRSF6 = Fas
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('CCR7', 'SELL', 'GZMB', 'CCR5', 'IL2RA', 'PTPRC', 'IL7R', 'CTLA4', 'FAS', 'CD28', 'CD27', 'ITGA4', 'ITGB7', 'ITGB1'), 'Effector vs. Memory')

	# Integrins:
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('ITGA2', 'ITGA2B', 'ITGA4', 'ITGAD', 'ITGAE', 'ITGAL', 'ITGAL', 'ITGAX', 'ITGB1', 'ITGB2', 'ITGB3', 'ITGB4', 'ITGB5', 'ITGB6', 'ITGB7', 'ITGB8'), 'Effector vs. Memory')

	#CD8 Activation
	# XCL1 = ENSMMUG00000060218
	# CCL4 = ENSMMUG00000008111
	# LOC100423131 = XCL1, ENSMMUG00000013779, Lymphotactin
	# LOC100430627 = CCL4, ENSMMUG00000008111
	# LOC100426632 = C-C motif chemokine 4
	# LOC100426537 = C-C motif chemokine 3-like
	# LOC114673087 = C-C motif chemokine 3-like 1
	# LOC100429751 = C-X-C chemokine receptor type 1-like
	# LOC701946 = C-X-C motif chemokine 5
	# LOC703222 = C-X-C motif chemokine 5
	# LOC100423954 = C-C motif chemokine 3-like 1
	# CD154 = CD40L
	# LAMP1 = CD107a, cytotoxic capacity
	# NR4A1 = Nur77
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('IFNG', 'CD69', 'TNF', 'NFKBID', 'LTB', 'CSF1', 'CSF2', 'CSF3', 'IL2', 'TNFRSF9', 'CCL4L1', 'NR4A3', 'TNFSF14', 'CD82', 'PIGT', 'IRF8', 'IRF4', 'IRF2', 'RGCC', 'PD1', 'PDCD1', 'TNFAIP3', 'LOC100423131', 'LOC100430627', 'ENSMMUG00000013779', 'XCL1', 'ENSMMUG00000060218', 'CCL4', 'CCL3', 'PLEK', 'NR4A2', 'LOC100426537', 'LOC114673087', 'KLF10', 'GADD45B', 'CD154', 'LAMP1', 'NR4A1', 'TIMP1', 'LOC703222', 'LOC701946', 'LOC100429751', 'LOC100423954'), 'CD8 Activation Markers')

	# LEF1 = naive
	# STAT1 = Th1 helper
	# CD40LG, GATA3 = Th2
	# FOXP3 = Treg
	# AQP3, GPR183 = Tcm
	# ANXA1, GPR183 = Tem
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('CD4', 'SELL', 'LEF1', 'STAT1', 'CD40LG', 'GATA3', 'FOXP3', 'AQP3', 'GPR183', 'HOPX', 'ITGB2', 'AHNAK', 'ANXA1'), 'CD4 Phenotypic Markers')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('PRF1', 'GNLY', 'NKG7', 'GZMA','GZMB','GZMH','GZMK','GZMM'), 'Cytotoxicity')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('MS4A1', 'CD79A', 'CD74', 'DRA'), 'B-cell Markers')

	# Also: CD19-, MS4A1-, CD79B-
	PlotMarkerSeries(seuratObj, reductions = reductions, title = 'Plasma Cells', features = c('CD79A', 'JCHAIN', 'MZB1', 'XBP1', 'CD79A'))

	PlotMarkerSeries(seuratObj, reductions = reductions, title = 'Monocyte', features = c('LYZ', 'CST3', 'S100A6', 'VIM'))

	# ITGB3 = CD61. LOC703451=PF4
	PlotMarkerSeries(seuratObj, reductions = reductions, title = 'Platelet/MK', features = c('ITGB3', 'PPBP', 'LOC703451'))

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('PPBP', 'PF4'), 'Megakaryocytes')

	PlotMarkerSeries(seuratObj, reductions = reductions, title = 'Stemness', features = c('CD34'))

	# MRC1 = CD206
	# ITGAM = CD11b, AM=CD11b-, Non-AM=CD11b+
	# Non-AMs: CD16+/CD206-/HLA-DR+/CD11b+
	# M1/M2: CCR7 (M1), CD163 (M2), https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7050058/
	# CD163, CD68 = macrophage markers
	# NOS2 = iNOS
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('CD14', 'FCGR3A', 'S100A8', 'S100A6', 'MARCO', 'MRC1', 'CD163', 'CHIT1', 'APOBEC3A', 'ITGAM', 'HLA-DRB1', 'CCR7', 'CD68', 'NOS2'), 'Myeloid')

	# IL3RA = CD123 / pDC
	# CLEC4C = CD303 / pDC
	# CD1c = mDC
	# THBD = CD141 / mDC
	# CD80, CD86 = co-stimulatory. low expression = tolerogenic
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('CD14', 'FCGR3A', 'IL3RA', 'CLEC4C', 'CD1C', 'THBD', 'CD80', 'CD86', 'TGFB1'), 'DCs')

	PlotMarkerSeries(seuratObj, reductions = reductions, title = 'Regulatory Markers', features = c('CD4', 'FOXP3', 'IL2RA', 'NR4A1'))

	PlotMarkerSeries(seuratObj, reductions = reductions, title = 'Treg17', features = c('RORA', 'RORB', 'RORC', 'IL4', 'STAT3'))

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('HAVCR2'), 'Th1')

	# ZBTB16 = PLZF
	PlotMarkerSeries(seuratObj, reductions = reductions, features = GetGeneSet('MAIT_Markers'), 'MAIT')

	# ZNF683 = HOBIT
	# LOC100423131 = XCL1, ENSMMUG00000013779, Lymphotactin
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('NCR3', 'CD69', 'KLRC2', 'XCL1', 'LOC100423131', 'ZNF683', 'CD7'), 'NKT')

	# ZNF683 = HOBIT
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('TBX21', 'GATA3', 'RORC', 'FOXP3', 'BCL6', 'EOMES', 'TOX', 'GATA2', 'TCF7', 'KLF2', 'NR4A1', 'LEF1', 'PRDM1', 'ID2', 'ID3', 'ZNF683', 'BHLHE40', 'EGR1', 'EGR2', 'EGR3'), 'Transcription Factors')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('TIGIT', 'CTLA4', 'BTLA', 'PDCD1', 'CD274'), 'Inhibitory Markers')

	# https://www.frontiersin.org/articles/10.3389/fimmu.2016.00076/full#:~:text=The%20three%20main%20pathways%20activated,activation%20(1%2C%202).

	PlotMarkerSeries(seuratObj, reductions = reductions, title = 'STATs', features = c('STAT1', 'STAT2', 'STAT3', 'STAT4'))

	PlotMarkerSeries(seuratObj, reductions = reductions, title = 'NFATs', features = c('NFATC1', 'NFATC2', 'NFAT3', 'NFATC4', 'NFAT5'))

	#DAP10/12
	# LOC707555 = ZAP70
	PlotMarkerSeries(seuratObj, reductions = reductions, title = 'Signaling', features = c('HCST', 'TYROBP', 'SYK', 'ZAP70', 'LOC707555', 'ITK'))

	#LILR/KIR:
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('LILRA5','LILRA6','LILRB4','LILRB5','KIR2DL4','KIR3DX1', 'MAMU-KIR', 'KIR2DL4', 'KIR3DL2'), 'LILR/KIR')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('FCGR1A','FCGR2A','FCGR2B','FCGR3', 'FCGR3A', 'FCGR3B', 'FCER1G'), 'FCGR')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = GetGeneSet('EffectorCytokines'), 'Effector Cytokines')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = GetGeneSet('Apoptosis'), 'Apoptosis')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = GetGeneSet('InterferonAndCytokines'), 'Interferon And Cytokines')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = GetGeneSet('Interleukins'), 'Interleukins')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = GetGeneSet('InterleukinReceptors'), 'Interleukin Receptors')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = GetGeneSet('ExhaustionOrInhibitory'), 'Exhaustion Or Inhibitory')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = GetGeneSet('MHC-II'), 'MHC-II')

	PlotMarkerSeries(seuratObj, features = c('TGFB1', 'TGFB2', 'TGFB3'), title = 'TGFB')

	PlotMarkerSeries(seuratObj, features = c('TGFBR1', 'TGFBR2', 'TGFBR3'), title = 'TGFB Receptor')

	# KLRC2 = ENSMMUG00000050862
	klrs <- c('KLRB1', 'KLRC1', 'KLRD1', 'KLRF1', 'KLRF2', 'KLRG1', 'KLRG2', 'KLRC2', 'KLRC3', 'KLRK1', 'ENSMMUG00000050862')
	PlotMarkerSeries(seuratObj, reductions = reductions, features = klrs, 'KLRs')

	# NCR1 = NKp46
	# NCR2 = NKp44
	# NCR3 = NKp30
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('NCR1', 'NCR2', 'NCR3'), 'NCRs')

	# ITGAE = CD103
	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('ITGAE', 'ITGB7', 'CD69', 'CXCR6', 'TYROBP'), 'Resident Memory')

	#chemokines
	chemokines <- c('CCL1','CCL11','CCL13','CCL16','CCL17','CCL18','CCL19','CCL2','CCL20','CCL21','CCL22','CCL24','CCL25','CCL26','CCL27','CCL28','CCL5','CCL7','CCL8')
	chemokines <- c(chemokines, c('CCR1','CCR2','CCR3','CCR4','CCR5','CCR6','CCR7','CCR8','CCR9','CCR10','CCRL2'))
	chemokines <- c(chemokines, c('CXCL1','CXCL10','CXCL11','CXCL12','CXCL13','CXCL14','CXCL16','CXCL17','CXCL5','CXCL6','CXCL8','CXCL9','CXCR1','CXCR2','CXCR3','CXCR4','CXCR5','CXCR6','XCR1'))

	PlotMarkerSeries(seuratObj, reductions = reductions, features = chemokines, 'Chemokines/Receptors')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('MKI67', 'TOP2A'), 'Cell Proliferation')

	#PlotMarkerSeries(seuratObj, reductions = reductions, features = c('EPCAM'), 'Epithelial Cells')

	# LOC710951 = TRAC
	# LOC114677140 = TRBC1
	# LOC711031 = TRDC
	# LOC720538 = TRGC1
	# LOC705095 = TRGC2
	PlotMarkerSeries(seuratObj, reductions = reductions, features = GetGeneSet('MMul10TcrConstantRegion'), 'TCR Constant Region')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('MMP2','COL1A1','COL1A2','COL5A1','LUM','PDGFRA'), 'Stromal')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('CCR9','LILRA4'), 'pDC')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('CD1C','THBD'), 'pDC')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('CDH1','FLT1'), 'Epithelial')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('LYZ','CSF1R','MSR1','MAFB','CD300E'), 'MoMacDC')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('CSF3R','FCGR3B'), 'Neutrophils')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = c('HBB','HBA2','HBA1'), 'Erythrocyte')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = GetGeneSet('MemoryAndNaive'), 'Putative Memory and Naive')

	PlotMarkerSeries(seuratObj, reductions = reductions, features = GetGeneSet('TandNK_Activation.1'), 'TandNK_Activation.1')

	# IGHD, IGHM = naive B and geminal center?
	# IGHA1 = plasma cells
	# IGKC = plasma cells
}


#' @title PlotMarkerSeries
#'
#' @description Iteratively plots a set of markers
#' @param seuratObj The seurat object
#' @param features A vector of feature names
#' @param reductions The reductions to plot
#' @param title An optional title of this plot series
#' @param setSize The maximum number of features to include per FeaturePlot
#' @export
#' @import Seurat
PlotMarkerSeries <- function(seuratObj, features, reductions = c('umap'), title = NULL, setSize = 4) {
	featuresToPlot <- .FindPlotableFeatures(seuratObj, features)
	if (is.null(featuresToPlot)) {
		print('None of the requested features were found, aborting')
		return()
	}

	steps <- ceiling(length(featuresToPlot) / setSize) - 1

	for (i in 0:steps) {
		start <- (i * 4) + 1
		end <- min((start + 3), length(featuresToPlot))
		genes <- featuresToPlot[start:end]


		PlotMarkerSet(seuratObj, reductions = reductions, title = title, features = genes)
	}
}

.FindPlotableFeatures <- function(seuratObj, features) {
	tryCatch({
		df <- Seurat::FetchData(seuratObj, vars = unique(features), cells = 1)
		return(unique(names(df)))
	}, error = function(e){
		return(NULL)
	})
}

.RemoveUnchangedOrZero <- function(seuratObj, reduction, features) {
	ret <- c()
	#Remove zeros or unchanged:
	dims <- paste0(Key(object = seuratObj[[reduction]]), c(1,2))
	data <- FetchData(object = seuratObj, vars = c(dims, features), cells = colnames(x = seuratObj))
	for (feature in features) {
		if (!feature %in% names(data)) {
			stop(paste0('Feature not found: ', feature, '. Presend: ', paste0(names(data), collapse = ',')))
		}

		if (sum(data[,feature] > 0) > 1 && length(unique(data[, feature])) > 1) {
			ret <- c(ret, feature)
		}
	}

	return(ret)
}

#' @import Seurat
#' @import patchwork
PlotMarkerSet <- function(seuratObj, reductions, title, features) {
	P <- NULL

	featuresToPlot <- .FindPlotableFeatures(seuratObj, features)
	if (is.null(featuresToPlot)) {
		print('None of the requested features were found, aborting')
		return()
	}

	featuresToPlotNonZero <- .RemoveUnchangedOrZero(seuratObj, reduction, featuresToPlot)

	if (length(featuresToPlotNonZero) != length(featuresToPlot)){
		missingFeats <- featuresToPlot[!(featuresToPlot %in% featuresToPlotNonZero)]
		print(paste0('The following features were requested, but not present: ', paste0(unique(missingFeats), collapse = ',')))
	}

	if (length(featuresToPlotNonZero) == 0){
		print('None of the requested features were present, skipping')
		return()
	}

	for (reduction in reductions) {
		P1 <- FeaturePlot(seuratObj, features = featuresToPlotNonZero, reduction = reduction, min.cutoff = 'q05', max.cutoff = 'q95')
		if (all(is.null(P))) {
			P <- P1
		} else {
			P <- P | P1
		}
	}

	if (!all(is.null(P))) {
		if (!is.null(title)) {
			P <- P + patchwork::plot_annotation(title = title)
		}

		print(P)
	}
}


pkg.env$GENE_SETS <- list()

.RegisterGeneSet <- function(name, genes) {
	if (name %in% names(pkg.env$GENE_SETS)) {
		stop(paste0('Name has already been registered: ', name))
	}

	pkg.env$GENE_SETS[[name]] <- genes
}

#' @title GetGeneSet
#'
#' @description Returns a vector with the set of genes registered under the provided name.
#' @param name The name of the gene set
#' @export
#' @import Seurat
GetGeneSet <- function(name) {
	if (! (name %in% names(pkg.env$GENE_SETS))) {
		warning(paste0('Unknown gene set: ', name))
	}

	return(pkg.env$GENE_SETS[[name]])
}

### Begin Phenotyping Gene Sets
.RegisterGeneSet('TandNK_Activation.1', c('CCL4L1','MIR155HG','RGCC','NFKBIA','IFNG','NR4A3','TNFSF14','CCL3'))

.RegisterGeneSet('TandNK_Activation.Core', c('CCL4L1','LOC100430627', 'IFNG','NR4A3','TNFSF14','CCL3', 'LOC100423131'))

# This is based on T-cell analysis from lung T/NK cells. This my not be precisely Memory/Naive.
.RegisterGeneSet('MemoryAndNaive', c('SELL', 'IL7R', 'LTB', 'SPOCK2', 'COTL1', 'JUNB', 'GPR183', 'CCR7', 'FUOM', 'CD7', 'PECAM1'))

.RegisterGeneSet('Cytotoxicity', c('PRF1', 'GNLY', 'NKG7', 'GZMA','GZMB','GZMH','GZMK','GZMM'))

.RegisterGeneSet('Cytotoxicity.GzmABH', c('GZMA','GZMB','GZMH'))
.RegisterGeneSet('Cytotoxicity.GzmKM', c('GZMK','GZMM'))
.RegisterGeneSet('Metallothionein', c('MT1A', 'MT1E', 'MT1M', 'MT1X', 'MT2A', 'MT1JP'))
.RegisterGeneSet('Metallothionein.Core', c('MT1E', 'MT1M', 'MT1X', 'MT2A'))
.RegisterGeneSet('IEGs', c('FOS', 'JUN', 'ZFP36'))
.RegisterGeneSet('TCellMemoryS100', c('S100A4', 'S100A6','S100A10','S100A11'))

.RegisterGeneSet('EffectorT', c('CCL4L1','CCL5','CCR7-','CD7-','FUOM-','GNLY','GZMB','GZMH','HOPX','LTB-','NKG7','PECAM1-','PRF1','RGS9','S100A4','SELL-','SPOCK2-','JUNB-'))

.RegisterGeneSet('CentralMemT', c('CCR7', 'CLDND1', 'GPR183'))

# Note: 'CD7', 'CCR7' excluded since they are in Tcm as well
.RegisterGeneSet('NaiveT', c('CTSH', 'CA6', 'GSTT1', 'LEF1','RGS10', 'TMIGD2'))

.RegisterGeneSet('CD8Memory', c('CD8A', 'CST7', 'CTSW', 'GNLY', 'GZMK'))

.RegisterGeneSet('InterferonAndCytokines', c('IFNA1', 'IFNA2', 'IFNB1', 'IL6', 'IFNG', 'TNF', 'CCL3', 'CCL3L1', 'CCL4', 'CCL4L1'))

.RegisterGeneSet('Interleukins', c('IL1', 'IL2', 'IL3', 'IL4', 'IL5', 'IL6', 'IL7', 'CXCL8', 'IL9', 'IL10', 'IL11', 'IL13', 'TXLNA', 'IL15', 'IL16', 'IL17A', 'IL17B', 'IL18', 'IL19', 'IL20', 'IL21', 'IL22'))

.RegisterGeneSet('InterleukinReceptors', c('IL1R1','IL1R2','IL2RA','IL2RB','IL2RG','IL3RA','IL4R','IL5RA','IL6R','IL7R','IL10RA','IL12RB1','IL12RB2','IL13RA2','IL15RA','IL17RA','IL17RB','IL17RC','IL17RD','IL17RE','IL20RA','IL20RB','IL21R','IL22RA2','IL27RA','IL31RA'))

# https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2024.1389630/full
# TNFRSF6 = Fas
.RegisterGeneSet('Apoptosis', c('FASLG', 'FAS', 'BCL2', 'IRF9', 'JAK3', 'CASP3', 'CASP8', 'BAX', 'TNFRSF1A', 'TP53', 'CYCS', 'PARP1'))

.RegisterGeneSet('EffectorCytokines', c('IFNG', 'TNF', 'IL17A', 'IL2', 'CSF1', 'CSF2', 'CCL3', 'CCL3L1', 'CCL4', 'CCL4L1'))
.RegisterGeneSet('ExhaustionOrInhibitory', c('PDCD1', 'TIGIT', 'HAVCR2', 'LAG3', 'CTLA4', 'VTCN1', 'CD244', 'KLRG1', 'TNFRSF14', 'BTLA', 'CD160'))

.RegisterGeneSet('Myelocytes', c("OLFM4", "LTF", "CAMP", "LCN2"))
.RegisterGeneSet('Pro_Myelocytes', c("BPI", "MPO", "ELANE", "RTD1A", "RTD1B", "DEFA1B", "AZU1"))
#These were found by literature investigation, and the MCMs were added via correlation with UNG - a necessary gene for somatic hypermutation.
.RegisterGeneSet('SomaticHypermutation', c("AICDA", "UNG", "MCM5", "MCM4", "MCM3", "MCM6", "MCM2"))
#These BCell signatures were compiled by Greg for use in RIRA (evidence is literature support & marker-uniqueness within RIRA)
.RegisterGeneSet("GerminalCenter", c("AICDA", "RGS13", "MEF2B", "ELL3", "FAS")) #CD40 omitted
.RegisterGeneSet("CD40_Pos_BCells", c("CD69", "AHNAK", "MAMU-DRA", "MS4A1", "CD1C")) #this is a weak signature, with some overlap between CD40- B cells. This should be supplemented with Memory vs Naive gene sets.
.RegisterGeneSet("CD40_Neg_BCells", c("ANHAK", "SOX5", "MAMU-DRA", "MS4A1", "CD1C")) #similarly, this signature has overlap with CD40+. 
.RegisterGeneSet("Pre-BCells", c("STMN1", "VPREB1", "SOX4"))
.RegisterGeneSet("PlasmaCells", c("DERL3", "TNFRSF17", "JCHAIN", "FKBP11", "XBP1")) #I omitted IGHA1, since we should probably keep this gene set class agnostic.
.RegisterGeneSet("NaiveB",  c("TCL1A", "FCER2", "CXCR4")) 
.RegisterGeneSet("MemoryB",  c("CR2", "PLAC8", "LY86", "CD44")) 
.RegisterGeneSet("InnateB",  c("CD1C", "AHNAK", "SOX5")) 
### End Phenotyping Gene Sets

# Sustaining CD4 survival?
# TNFRSF4 = OX40
# TNFSF4 = OX40L

.RegisterGeneSet('MMul10TcrConstantRegion', c('LOC710951', 'LOC114677140', 'LOC711031', 'LOC720538', 'LOC705095'))

.RegisterGeneSet('MAIT_Markers', c('KLRB1', 'CEPBD', 'NCR3', 'ZBTB16', 'RORC', 'SLC4A10', 'DPP4'))

# Dysfunction?
# 'MT1A', 'MT2A', 'MT1M'
.RegisterGeneSet('MMul10TcrGenes', c('LOC703029','LOC696306','LOC106999340','LOC106996262','LOC106999345','LOC106997707','LOC106997706','LOC710149','LOC700771','LOC699427','LOC711871','LOC709081','LOC698785','LOC114677052','LOC114676933','LOC106999353','LOC106999351','LOC106999350','LOC106999349','LOC106999348','LOC106999347','LOC106999346','LOC106999343','LOC106999341','LOC106999339','LOC106999337','LOC106999336','LOC106999335','LOC106999312','LOC106997705','LOC106997704','LOC106997703','LOC106997702','LOC106997697','LOC106997453','LOC106997452','LOC106997451','LOC106995765','LOC106992460','LOC106992446','LOC106992434','LOC106992433','LOC720456','LOC716949','LOC716866','LOC711537','LOC711386','LOC711194','LOC711141','LOC711066','LOC710821','LOC710627','LOC710455','LOC710361','LOC710183','LOC710093','LOC709531','LOC708581','LOC708328','LOC704883','LOC703153','LOC702904','LOC702550','LOC702113','LOC701992','LOC701875','LOC701745','LOC701395','LOC701262','LOC701152','LOC700224','LOC700154','LOC700105','LOC699912','LOC699790','LOC699543','LOC699298','LOC699162','LOC698913','LOC698543','LOC698289','LOC698161','LOC697792','LOC697466','LOC697234','LOC697054','LOC696752','LOC696684','LOC696557','LOC696075','LOC695943','LOC114679533','LOC114679531','LOC114677139','LOC114677137','LOC114677136','LOC114677055','LOC114677054','LOC114677050','LOC114677049','LOC114677047','LOC114675766'))
.RegisterGeneSet("MMul10_MHC", c("MAMU-A", "MAMU-A3", "MAMU-AG", "LOC719250", "LOC100426197", "LOC714466", "LOC694372", "LOC114677644", "LOC114675360", "LOC106997902", "LOC698738", "LOC106996077", "LOC106992378", "LOC100424348", "LOC106995519", "LOC106992468", "LOC114676051", "LOC106996627", "LOC106997893", "LOC106997885", "LOC114675646", "LOC114669738", "LOC720132", "LOC719702", "LOC106996676", "LOC714964", "LOC114669810", "LOC114675357", "LOC100428435", "LOC100429195", "LOC699987", "LOC114677642", "LOC723473", "LOC106995461", "LOC106995452"))
.RegisterGeneSet("MMul10_KIR", c("KIR3DL0", "KIR3DL2", "KIR2DS4", "KIR2DL4", "KIR3DS05", "KIR3DL12", "KIR3DL1", "KIR3DH", "KIR3DH5", "KIR3DL21", "KIR3DL11", "MAMU-KIR", "LOC106994859", "LOC114669735", "LOC100424026", "LOC100125572"))

# Generated using: sort(rownames(seuratObj@assays$RNA)[grepl(rownames(seuratObj@assays$RNA), pattern = '^RP[SL]')])
ribosomalGenes <- c("RPL10", "RPL10A", "RPL10L", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL22L1", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL26L1", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL3", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL36AL", "RPL37", "RPL37A", "RPL37A-1", "RPL38", "RPL39", "RPL39L", "RPL3L", "RPL4", "RPL41", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL7L1", "RPL8", "RPL9", "RPLP0", "RPLP1", "RPLP2", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS19BP1", "RPS2", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS27A-1", "RPS27L", "RPS28", "RPS29", "RPS3", "RPS4X", "RPS4Y1", "RPS4Y2", "RPS5", "RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KA6", "RPS6KB1", "RPS6KB2", "RPS6KC1", "RPS6KL1", "RPS7", "RPS8", "RPS9", "RPSA")
  
# These were identified by grep on the descriptions of all LOC* genes for ribosomal
ribosomalGenes <- c(ribosomalGenes, c("LOC717349", "LOC718596", "LOC705988", "LOC703571", "LOC100429589", "LOC114672439", "LOC711644", "LOC705840", "LOC114674534", "LOC697510", "LOC100426562", "LOC100430397", "LOC705494", "LOC718356", "LOC708673", "LOC708221", "LOC693728", "LOC701224", "LOC695949", "LOC100430766", "LOC708539", "LOC114670541", "LOC694687", "LOC696921", "LOC106999910", "LOC705621", "LOC706532", "LOC701369", "LOC106994001", "LOC114675582", "LOC706087", "LOC711318", "LOC695920", "LOC697723", "LOC700257", "LOC710331", "LOC700807", "LOC106992793", "LOC717846", "LOC694712", "LOC708815", "LOC702606", "LOC707204", "LOC106996615", "LOC711714", "LOC697734", "LOC100425632", "LOC709255", "LOC703037", "LOC710558", "LOC114670812", "LOC721996", "LOC716997", "LOC698732", "LOC702174", "LOC114679384", "LOC697431", "LOC704429", "LOC106992303", "LOC114676316", "LOC100430733", "LOC708154", "LOC114673092", "LOC697548", "LOC695895", "LOC708833", "LOC701909", "LOC715154", "LOC717289", "LOC106998255", "LOC701084", "LOC707369", "LOC697035", "LOC720167", "LOC700616", "LOC700401", "LOC106996191", "LOC106998695", "LOC707090", "LOC717450", "LOC693713", "LOC713382", "LOC716142", "LOC696119", "LOC714106", "LOC695421", "LOC709128", "LOC714782", "LOC693486", "LOC106992469", "LOC699246", "LOC708791", "LOC699630", "LOC697141", "LOC710726", "LOC107000432", "LOC100428100", "LOC718153", "LOC100427366", "LOC693343", "LOC114669718", "LOC713514", "LOC698010", "LOC715166", "LOC695327", "LOC702713", "LOC704640", "LOC114674200", "LOC700134", "LOC715041", "LOC722876", "LOC710644", "LOC705596", "LOC702961", "LOC704054", "LOC715668", "LOC717801", "LOC707085", "LOC114679603", "LOC106999518", "LOC718267", "LOC705327", "LOC106992281", "LOC695218", "LOC702957", "LOC708147", "LOC719947", "LOC701392", "LOC699166", "LOC106999718", "LOC702760", "LOC702328", "LOC704267", "LOC719289", "LOC106997304", "LOC100423207", "LOC100428069", "LOC106999080", "LOC100427696", "LOC712160", "LOC712267", "LOC100425072", "LOC100427254", "LOC696068", "LOC694967", "LOC708492", "LOC710583", "LOC717763", "LOC693820", "LOC114678655", "LOC114678366", "LOC114678384", "LOC114671497", "LOC114675033", "LOC697719", "LOC694071", "LOC695523", "LOC703455", "LOC715814", "LOC106998837", "LOC713060", "LOC693904", "LOC711043", "LOC709241", "LOC693584", "LOC698301", "LOC703631", "LOC701292", "LOC706910", "LOC712630", "LOC714801", "LOC717838", "LOC698099", "LOC694067", "LOC706082", "LOC702514", "LOC711324", "LOC703616", "LOC704396", "LOC712274", "LOC698273", "LOC718556", "LOC715043", "LOC696015", "LOC698297", "LOC700262", "LOC695670", "LOC698383", "LOC708603", "LOC701421", "LOC697065", "LOC114677128", "LOC714656", "LOC106996346", "LOC703683", "LOC701691", "LOC698792", "LOC107000151", "LOC697561", "LOC710477", "LOC100427278", "LOC720748", "LOC114672187", "LOC708140", "LOC698768", "LOC716735", "LOC114673268", "LOC114674609", "LOC114674676", "LOC114674882", "LOC114675055", "LOC114675456", "LOC114674615", "LOC114674626", "LOC114674635", "LOC114674636", "LOC114674637", "LOC114674638", "LOC114674640", "LOC114674642", "LOC114674650", "LOC114674666", "LOC114674683", "LOC114673870", "LOC114675629", "LOC114675635", "LOC114675637", "LOC114675638", "LOC114675641", "LOC114675642", "LOC114675643", "LOC114675636", "LOC114675644", "LOC114675630", "LOC114675633", "LOC114675631", "LOC114675632", "LOC114675739", "LOC114675743", "LOC114675745", "LOC114675746", "LOC114675747", "LOC114675740", "LOC114675741", "LOC114675856", "LOC114675854", "LOC114675853", "LOC114675857", "LOC114675858", "LOC114675859", "LOC114675860", "LOC114675913", "LOC114675919", "LOC114675922", "LOC114675923", "LOC114675924", "LOC114675925", "LOC114675914", "LOC114675915", "LOC114675916", "LOC114675918", "LOC114675141", "LOC718979", "LOC106993215", "LOC106995250", "LOC709583", "LOC114679608", "LOC717779", "LOC707414", "LOC709376", "LOC710038", "LOC706793", "LOC106999723", "LOC696958", "LOC106993579", "LOC106995074", "LOC694799", "LOC106995394", "LOC698530", "LOC107000213", "LOC705400", "LOC720019", "LOC698143", "LOC697585", "LOC713714", "LOC702620", "LOC715037", "LOC699544", "LOC100427498", "LOC699681", "LOC701302", "LOC700008", "LOC698713", "LOC702871", "LOC693856", "LOC701386", "LOC716126", "LOC696423", "LOC114669857", "LOC697476", "LOC100429843", "LOC696901", "LOC693272", "LOC704083", "LOC702723", "LOC703433", "LOC722545", "LOC696877", "LOC106998431", "LOC106999185", "LOC706492", "LOC699695", "LOC713984", "LOC703435", "LOC709778", "LOC694122", "LOC711449", "LOC701763", "LOC719242", "LOC715711", "LOC700872", "LOC717246", "LOC698182", "LOC701156", "LOC106996680", "LOC710859", "LOC106997591", "LOC707834", "LOC100428649", "LOC695340", "LOC718357", "LOC702542", "LOC718020", "LOC699687", "LOC711174", "LOC709462", "LOC712987", "LOC106995957", "LOC708800", "LOC717286", "LOC698999", "LOC696145", "LOC106999371", "LOC698377", "LOC709045", "LOC711760", "LOC114680233", "LOC107000909", "LOC106994077", "LOC701710", "LOC114674062", "LOC699398", "LOC693578", "LOC106992566", "LOC106997285", "LOC718737", "LOC706606", "LOC704012", "LOC702297", "LOC114680376", "LOC696154", "LOC721751", "LOC704238", "LOC709820", "LOC114676029", "LOC708995", "LOC710502", "LOC695715", "LOC698942", "LOC696134", "LOC696785", "LOC698684", "LOC695740", "LOC706177", "LOC717053", "LOC698197", "LOC714598", "LOC106998794", "LOC704510", "LOC698130", "LOC107000483", "LOC716320", "LOC106994271", "LOC708074", "LOC100425638", "LOC709201", "LOC707358", "LOC706712", "LOC701457", "LOC114676221", "LOC114677736", "LOC100423145", "LOC708782", "LOC714171", "LOC707437", "LOC701466", "LOC719029", "LOC697987", "LOC714120", "LOC708085", "LOC694937", "LOC106993116", "LOC706340", "LOC697214", "LOC106995063", "LOC699344", "LOC114676734", "LOC702875", "LOC695122", "LOC717674", "LOC712555", "LOC711491", "LOC702677", "LOC704365", "LOC708118", "LOC709678", "LOC106994904", "LOC703784", "LOC702847", "LOC709086", "LOC114674318", "LOC693848", "LOC719244", "LOC114680572", "LOC719770", "LOC710317", "LOC714620", "LOC703648", "LOC703365", "LOC706181", "LOC700425", "LOC708644", "LOC709650", "LOC705067", "LOC713497", "LOC106994132", "LOC703919", "LOC114671097", "LOC703853", "LOC114674108", "LOC114676998", "LOC712599", "LOC114678374", "LOC114678397", "LOC114679023", "LOC699392", "LOC114678167", "LOC694228", "LOC710665", "LOC107000401", "LOC100423661", "LOC106999480", "LOC106997918", "LOC107000809", "LOC114673683", "LOC106998280", "LOC106998813", "LOC100429752", "LOC100423423", "LOC106994361", "LOC106994914", "LOC705174", "LOC707072", "LOC700862", "LOC693472", "LOC106998156", "LOC721192", "LOC707117", "LOC705979", "LOC705677", "LOC710287", "LOC704345", "LOC699239", "LOC702116", "LOC706754", "LOC703797", "LOC710877", "LOC696151", "LOC702932", "LOC695743", "LOC693287", "LOC703024", "LOC698657", "LOC100429815", "LOC707342", "LOC695463", "LOC699805", "LOC722113", "LOC694734", "LOC106998116", "LOC716075", "LOC703563", "LOC713469", "LOC714185", "LOC703251", "LOC703283", "LOC698417", "LOC704062", "LOC698226", "LOC114677434", "LOC693480", "LOC114676944", "LOC702978", "LOC100427953", "LOC100427759", "LOC704256", "LOC700045", "LOC702525", "LOC710595", "LOC698574", "LOC707182", "LOC700578", "LOC698448", "LOC709444", "LOC698086", "LOC714855", "LOC708405", "LOC114675505", "LOC697178", "LOC703216", "LOC716530", "LOC716150", "LOC694093", "LOC718478", "LOC707861", "LOC106998326", "LOC695279", "LOC106995044", "LOC100423364", "LOC100423295", "LOC704394", "LOC703794", "LOC693573", "LOC699375", "LOC701250", "LOC706496", "LOC693947", "LOC713991", "LOC710642", "LOC109910387", "LOC109910386", "LOC106997745", "LOC714700", "LOC100426636", "LOC106998409"))
.RegisterGeneSet("MMul10_Ribosomal", ribosomalGenes)

# Identified by grep on the descriptions of all LOC* genes for 'mitochond'
mitochondrialGenes <- c("LOC703370", "LOC721815", "LOC720015", "LOC703345", "LOC708979", "LOC704978", "LOC717349", "LOC718596", "LOC705988", "LOC703571", "LOC100429589", "LOC114672439", "LOC711644", "LOC705840", "LOC697510", "LOC100426562", "LOC100430397", "LOC705494", "LOC718356", "LOC708673", "LOC708221", "LOC693728", "LOC701224", "LOC695949", "LOC716130", "LOC100430580", "LOC100423136", "LOC106998857", "LOC701346", "LOC695671", "LOC696244", "LOC100427705", "LOC114673694", "LOC106997368", "LOC107000220", "LOC696031", "LOC708562", "LOC716946", "LOC715928", "LOC713542", "LOC719688", "LOC699117", "LOC694182", "LOC696890", "LOC707114", "LOC704076", "LOC704440", "LOC709495", "LOC702581", "LOC711636", "LOC698075", "LOC711183", "LOC698715", "LOC694260", "LOC707476", "LOC700815", "LOC711417", "LOC706691", "LOC696144", "LOC695041", "LOC701131", "LOC114676443", "LOC712886", "LOC114669858", "LOC106996230", "LOC701964", "LOC705178", "LOC700087", "LOC695277", "LOC695732", "LOC114678889", "LOC106995294", "LOC100425230", "LOC100425700", "LOC697409", "LOC710274", "LOC707227", "LOC701661", "LOC696956", "LOC106996920", "LOC701236", "LOC106992283", "LOC717879", "LOC114674848", "LOC106998783", "LOC106998836", "LOC114678228", "LOC114673956", "LOC114677410", "LOC114678786", "LOC114679364", "LOC114679428", "LOC114680400", "LOC114671073", "LOC100423866", "LOC106993150", "LOC106996000", "LOC107000216", "LOC699137", "LOC106996051", "LOC713467", "LOC114675539", "LOC693843", "LOC697678", "LOC702427", "LOC114675509", "LOC698520", "LOC106997494", "LOC100425012", "LOC106999579", "LOC107000931", "LOC106997075", "LOC717923", "LOC106998133", "LOC712330", "LOC706686", "LOC715327", "LOC693799", "LOC693766", "LOC698522", "LOC114675319", "LOC718745", "LOC718390", "LOC114678805", "LOC706366", "LOC114674124", "LOC106996824", "LOC106994877", "LOC106992611", "LOC107000653", "LOC107000976", "LOC713028", "LOC693813", "LOC711899", "LOC718711", "LOC114669720", "LOC721639", "LOC698302", "LOC718712", "LOC698305", "LOC706710", "LOC705282", "LOC114670386", "LOC699779", "LOC100427337", "LOC100429934", "LOC704917", "LOC715730", "LOC705066", "LOC709894", "LOC711634", "LOC723069", "LOC705183", "LOC702309", "LOC701966", "LOC704336", "LOC697416", "LOC106993325", "LOC716730", "LOC717607", "LOC106997886", "LOC107000665", "LOC106999488", "LOC106998010", "LOC705476", "LOC114673096", "LOC696696", "LOC715179")
.RegisterGeneSet("MMul10_Mitochondrial", mitochondrialGenes)

igHeavyVariable <- c("LOC715358", "LOC114669771", "LOC106999620", "LOC106999619", "LOC114679735", "LOC114679700", "LOC114679702", "LOC114679730", "LOC114679696", "LOC106999621", "LOC114679699", "LOC721017", "LOC106995637", "LOC721043", "LOC721139", "LOC114679710", "LOC715165", "LOC720899", "LOC714894", "LOC106999623", "LOC106999617", "LOC714310", "LOC723533", "LOC715400", "LOC701407", "LOC720918", "LOC106999616", "LOC721119", "LOC114679704", "LOC715260", "LOC719378", "LOC106999631", "LOC721089", "LOC106999608", "LOC721132", "LOC720890", "LOC106999629", "LOC114679692", "LOC100430012", "LOC106999615", "LOC114669802", "LOC114679711", "LOC714939", "LOC715543", "LOC106999624", "LOC721104", "LOC114679713", "LOC721116", "LOC114679694", "LOC114679729", "LOC106999604", "LOC715594", "LOC106995900", "LOC106999610", "LOC106999609", "LOC721083", "LOC106999633", "LOC106996003", "LOC114679703", "LOC721002", "LOC106999628", "LOC114679693", "LOC722278", "LOC720985", "LOC720452", "LOC720904", "LOC114679728", "LOC114679707", "LOC720935", "LOC721108", "LOC720974", "LOC106999618", "LOC114679733", "LOC715733", "LOC114679736", "LOC714419")
igKappaVariable <- c("LOC106993058", "LOC721353", "LOC703375", "LOC701600", "LOC701068", "LOC106993063", "LOC707007", "LOC709162", "LOC709703", "LOC709066", "LOC106993071", "LOC709793", "LOC707181", "LOC106993142", "LOC106992444", "LOC106993055", "LOC703139", "LOC114671809", "LOC106993143", "LOC723660", "LOC703610", "LOC106993056", "LOC114671797", "LOC708976", "LOC106993144", "LOC106998680", "LOC106993140", "LOC106993064", "LOC106993141", "LOC706869", "LOC707609", "LOC708249", "LOC709890", "LOC106993048", "LOC700696")
igLambdaVariable <-c("LOC107000570", "LOC106992425", "LOC701976", "LOC708547", "LOC107000578", "LOC114670549", "LOC107000585", "LOC107000586", "LOC107000581", "LOC114670550", "LOC107000576", "LOC107000577", "LOC702421", "LOC701772", "LOC107000579", "LOC701240")
.RegisterGeneSet('MMul10_Ig_Variable', unique(c(igHeavyVariable, igKappaVariable, igLambdaVariable)))

## Human gene sets
# Generated using: sort(rownames(seuratObj[['RNA']])[grepl(rownames(seuratObj[['RNA']]), pattern = '^RP[SL]')])
humanRibosomalGenes <- read.table(system.file("extdata/humanRibosomalGenes_Human_GRCh38_p13_Ensembl.tsv", package = "RIRA"), header = T)$ribosomal_genes
.RegisterGeneSet("HumanRibosomalGenes", humanRibosomalGenes)
# Generated using: sort(rownames(seuratObj[["RNA"]])[grepl(rownames(seuratObj[["RNA"]]), pattern = '^MT-')])
humanMitochondrialGenes <- c('MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6','MT-RNR1','MT-RNR2','MT-TA','MT-TC','MT-TD','MT-TE','MT-TF','MT-TG','MT-TH','MT-TI','MT-TK','MT-TL1','MT-TL2','MT-TM','MT-TN','MT-TP','MT-TQ','MT-TR','MT-TS1','MT-TS2','MT-TT','MT-TV','MT-TW','MT-TY')
.RegisterGeneSet("HumanMitochondrialGenes", humanMitochondrialGenes)

humanTravGenes <- c('TRAV1-1','TRAV1-2','TRAV10','TRAV11','TRAV12-1','TRAV12-2','TRAV12-3','TRAV13-1','TRAV13-2','TRAV14DV4','TRAV15','TRAV16','TRAV17','TRAV18','TRAV19','TRAV2','TRAV20','TRAV21','TRAV22','TRAV23DV6','TRAV24','TRAV25','TRAV26-1','TRAV26-2','TRAV27','TRAV28','TRAV29DV5','TRAV3','TRAV30','TRAV31','TRAV32','TRAV33','TRAV34','TRAV35','TRAV36DV7','TRAV37','TRAV38-1','TRAV38-2DV8','TRAV39','TRAV4','TRAV40','TRAV41','TRAV5','TRAV6','TRAV7','TRAV8-1','TRAV8-2','TRAV8-3','TRAV8-4','TRAV8-5','TRAV8-6','TRAV8-7','TRAV9-1','TRAV9-2')
humanTrbvGenes <- c('TRBV1','TRBV10-1','TRBV10-2','TRBV10-3','TRBV11-1','TRBV11-2','TRBV11-3','TRBV12-1','TRBV12-2','TRBV12-3','TRBV12-4','TRBV12-5','TRBV13','TRBV14','TRBV15','TRBV16','TRBV17','TRBV18','TRBV19','TRBV2','TRBV20-1','TRBV20OR9-2','TRBV21-1','TRBV21OR9-2','TRBV22-1','TRBV22OR9-2','TRBV23-1','TRBV23OR9-2','TRBV24-1','TRBV24OR9-2','TRBV25-1','TRBV25OR9-2','TRBV26','TRBV26OR9-2','TRBV27','TRBV28','TRBV29-1','TRBV29OR9-2','TRBV3-1','TRBV30','TRBV4-1','TRBV4-2','TRBV5-1','TRBV5-2','TRBV5-3','TRBV5-4','TRBV5-5','TRBV5-6','TRBV5-7','TRBV6-1','TRBV6-2','TRBV6-4','TRBV6-5','TRBV6-6','TRBV6-7','TRBV6-8','TRBV7-1','TRBV7-2','TRBV7-3','TRBV7-4','TRBV7-5','TRBV7-6','TRBV7-7','TRBV7-9','TRBV8-1','TRBV8-2','TRBV9','TRBVA','TRBVB')
humanTrbdGenes <- c('TRBD1')
humanTrajGenes <- c('TRAJ1','TRAJ10','TRAJ11','TRAJ12','TRAJ13','TRAJ14','TRAJ16','TRAJ17','TRAJ18','TRAJ19','TRAJ2','TRAJ20','TRAJ21','TRAJ22','TRAJ23','TRAJ24','TRAJ25','TRAJ26','TRAJ27','TRAJ28','TRAJ29','TRAJ3','TRAJ30','TRAJ31','TRAJ32','TRAJ33','TRAJ34','TRAJ35','TRAJ36','TRAJ37','TRAJ38','TRAJ39','TRAJ4','TRAJ40','TRAJ41','TRAJ42','TRAJ43','TRAJ44','TRAJ45','TRAJ46','TRAJ47','TRAJ48','TRAJ49','TRAJ5','TRAJ50','TRAJ51','TRAJ52','TRAJ53','TRAJ54','TRAJ55','TRAJ56','TRAJ57','TRAJ58','TRAJ59','TRAJ6','TRAJ60','TRAJ61','TRAJ7','TRAJ8','TRAJ9')
humanTrbjGenes <- c('TRBJ1-1','TRBJ1-2','TRBJ1-3','TRBJ1-4','TRBJ1-5','TRBJ1-6','TRBJ2-1','TRBJ2-2','TRBJ2-2P','TRBJ2-3','TRBJ2-4','TRBJ2-5','TRBJ2-6','TRBJ2-7')
.RegisterGeneSet("HumanTcrGenes", c(humanTravGenes, humanTrbvGenes, humanTrajGenes, humanTrbjGenes, humanTrbdGenes))

# Generated via sort(rownames(seuratObj[["RNA"]])[grepl(rownames(seuratObj[["RNA"]]), pattern = '^IGHV')]).  
humanIgHeavyVariable <- c('IGHV1-12','IGHV1-14','IGHV1-17','IGHV1-18','IGHV1-2','IGHV1-24','IGHV1-3','IGHV1-45','IGHV1-46','IGHV1-58','IGHV1-67','IGHV1-68','IGHV1-69','IGHV1-69-2','IGHV1-69D','IGHV1OR15-1','IGHV1OR15-2','IGHV1OR15-3','IGHV1OR15-4','IGHV1OR15-6','IGHV1OR15-9','IGHV1OR16-1','IGHV1OR16-2','IGHV1OR16-3','IGHV1OR16-4','IGHV1OR21-1','IGHV2-26','IGHV2-5','IGHV2-70','IGHV2-70D','IGHV2OR16-5','IGHV3-11','IGHV3-13','IGHV3-15','IGHV3-16','IGHV3-19','IGHV3-20','IGHV3-21','IGHV3-22','IGHV3-23','IGHV3-25','IGHV3-29','IGHV3-30','IGHV3-30-2','IGHV3-32','IGHV3-33','IGHV3-33-2','IGHV3-35','IGHV3-36','IGHV3-37','IGHV3-38','IGHV3-41','IGHV3-42','IGHV3-43','IGHV3-47','IGHV3-48','IGHV3-49','IGHV3-50','IGHV3-52','IGHV3-53','IGHV3-54','IGHV3-57','IGHV3-6','IGHV3-60','IGHV3-62','IGHV3-63','IGHV3-64','IGHV3-64D','IGHV3-65','IGHV3-66','IGHV3-69-1','IGHV3-7','IGHV3-71','IGHV3-72','IGHV3-73','IGHV3-74','IGHV3-75','IGHV3-76','IGHV3-79','IGHV3OR15-7','IGHV3OR16-10','IGHV3OR16-11','IGHV3OR16-12','IGHV3OR16-13','IGHV3OR16-15','IGHV3OR16-16','IGHV3OR16-17','IGHV3OR16-6','IGHV3OR16-7','IGHV3OR16-8','IGHV3OR16-9','IGHV4-28','IGHV4-31','IGHV4-34','IGHV4-39','IGHV4-4','IGHV4-55','IGHV4-59','IGHV4-61','IGHV4-80','IGHV4OR15-8','IGHV5-10-1','IGHV5-51','IGHV5-78','IGHV6-1','IGHV7-27','IGHV7-34-1','IGHV7-4-1','IGHV7-40','IGHV7-56','IGHV7-81','IGHV8-51-1','IGHVII-1-1','IGHVII-15-1','IGHVII-22-1','IGHVII-26-2','IGHVII-28-1','IGHVII-30-1','IGHVII-30-21','IGHVII-33-1','IGHVII-40-1','IGHVII-43-1','IGHVII-44-2','IGHVII-46-1','IGHVII-49-1','IGHVII-51-2','IGHVII-53-1','IGHVII-60-1','IGHVII-62-1','IGHVII-65-1','IGHVII-67-1','IGHVII-74-1','IGHVII-78-1','IGHVIII-11-1','IGHVIII-13-1','IGHVIII-16-1','IGHVIII-2-1','IGHVIII-22-2','IGHVIII-25-1','IGHVIII-26-1','IGHVIII-38-1','IGHVIII-44','IGHVIII-47-1','IGHVIII-5-1','IGHVIII-5-2','IGHVIII-67-2','IGHVIII-67-3','IGHVIII-67-4','IGHVIII-76-1','IGHVIII-82','IGHVIV-44-1')
# Generated via sort(rownames(seuratObj[["RNA"]])[grepl(rownames(seuratObj[["RNA"]]), pattern = '^IGHJ')]).  
humanIgHeavyJoining <- c('IGHJ1','IGHJ1P','IGHJ2','IGHJ2P','IGHJ3','IGHJ3P','IGHJ4','IGHJ5','IGHJ6')
# Generated via sort(rownames(seuratObj[["RNA"]])[grepl(rownames(seuratObj[["RNA"]]), pattern = '^IGHD')]).  
humanIgHeavyDiversity <- c('IGHD','IGHD1-1','IGHD1-14','IGHD1-20','IGHD1-26','IGHD1-7','IGHD1OR15-1A','IGHD1OR15-1B','IGHD2-15','IGHD2-2','IGHD2-21','IGHD2-8','IGHD2OR15-2A','IGHD2OR15-2B','IGHD3-10','IGHD3-16','IGHD3-22','IGHD3-3','IGHD3-9','IGHD3OR15-3A','IGHD3OR15-3B','IGHD4-11','IGHD4-17','IGHD4-23','IGHD4-4','IGHD4OR15-4A','IGHD4OR15-4B','IGHD5-12','IGHD5-18','IGHD5-24','IGHD5-5','IGHD5OR15-5A','IGHD5OR15-5B','IGHD6-13','IGHD6-19','IGHD6-25','IGHD6-6','IGHD7-27')
# Generated via sort(rownames(seuratObj[["RNA"]])[grepl(rownames(seuratObj[["RNA"]]), pattern = '^IGL')]). I excluded IGLC genes to conserve the constant chain expression. 
humanIgLambda <- c('IGLJ1','IGLJ2','IGLJ3','IGLJ4','IGLJ5','IGLJ6','IGLJ7','IGLJCOR18','IGLL1','IGLL3P','IGLL4P','IGLL5','IGLON5','IGLV1-36','IGLV1-40','IGLV1-41','IGLV1-44','IGLV1-47','IGLV1-50','IGLV1-51','IGLV1-62','IGLV10-54','IGLV10-67','IGLV11-55','IGLV2-11','IGLV2-14','IGLV2-18','IGLV2-23','IGLV2-28','IGLV2-33','IGLV2-34','IGLV2-5','IGLV2-8','IGLV3-1','IGLV3-10','IGLV3-12','IGLV3-13','IGLV3-15','IGLV3-16','IGLV3-17','IGLV3-19','IGLV3-2','IGLV3-21','IGLV3-22','IGLV3-24','IGLV3-25','IGLV3-26','IGLV3-27','IGLV3-29','IGLV3-30','IGLV3-31','IGLV3-32','IGLV3-4','IGLV3-6','IGLV3-7','IGLV3-9','IGLV4-3','IGLV4-60','IGLV4-69','IGLV5-37','IGLV5-45','IGLV5-48','IGLV5-52','IGLV6-57','IGLV7-35','IGLV7-43','IGLV7-46','IGLV8-61','IGLV8OR8-1','IGLV9-49','IGLVI-20','IGLVI-38','IGLVI-42','IGLVI-56','IGLVI-63','IGLVI-68','IGLVI-70','IGLVIV-53','IGLVIV-59','IGLVIV-64','IGLVIV-65','IGLVIV-66-1','IGLVIVOR22-1','IGLVIVOR22-2','IGLVV-58','IGLVV-66','IGLVVI-22-1','IGLVVI-25-1','IGLVVII-41-1')
# Generated via sort(rownames(seuratObj[["RNA"]])[grepl(rownames(seuratObj[["RNA"]]), pattern = '^IGK')]). I excluded IGKC genes to conserve the constant chain expression. 
humanIgKappa <- c('IGKJ1','IGKJ2','IGKJ3','IGKJ4','IGKJ5','IGKV1-12','IGKV1-13','IGKV1-16','IGKV1-17','IGKV1-22','IGKV1-27','IGKV1-32','IGKV1-33','IGKV1-35','IGKV1-37','IGKV1-39','IGKV1-5','IGKV1-6','IGKV1-8','IGKV1-9','IGKV1D-12','IGKV1D-13','IGKV1D-16','IGKV1D-17','IGKV1D-22','IGKV1D-27','IGKV1D-32','IGKV1D-33','IGKV1D-35','IGKV1D-37','IGKV1D-39','IGKV1D-42','IGKV1D-43','IGKV1D-8','IGKV1OR-2','IGKV1OR-3','IGKV1OR1-1','IGKV1OR10-1','IGKV1OR2-1','IGKV1OR2-108','IGKV1OR2-11','IGKV1OR2-118','IGKV1OR2-2','IGKV1OR2-3','IGKV1OR2-6','IGKV1OR2-9','IGKV1OR22-1','IGKV1OR22-5','IGKV1OR9-1','IGKV1OR9-2','IGKV2-10','IGKV2-14','IGKV2-18','IGKV2-19','IGKV2-23','IGKV2-24','IGKV2-26','IGKV2-28','IGKV2-29','IGKV2-30','IGKV2-36','IGKV2-38','IGKV2-4','IGKV2-40','IGKV2D-10','IGKV2D-14','IGKV2D-18','IGKV2D-19','IGKV2D-23','IGKV2D-24','IGKV2D-26','IGKV2D-28','IGKV2D-29','IGKV2D-30','IGKV2D-36','IGKV2D-38','IGKV2D-40','IGKV2OR2-1','IGKV2OR2-10','IGKV2OR2-2','IGKV2OR2-7','IGKV2OR2-7D','IGKV2OR2-8','IGKV2OR22-3','IGKV2OR22-4','IGKV3-11','IGKV3-15','IGKV3-20','IGKV3-25','IGKV3-31','IGKV3-34','IGKV3-7','IGKV3D-11','IGKV3D-15','IGKV3D-20','IGKV3D-25','IGKV3D-31','IGKV3D-34','IGKV3D-7','IGKV3OR2-268','IGKV3OR2-5','IGKV3OR22-2','IGKV4-1','IGKV5-2','IGKV6-21','IGKV6D-21','IGKV6D-41','IGKV7-3')
.RegisterGeneSet("HumanBcrGenes", c(humanIgHeavyVariable, humanIgHeavyJoining, humanIgHeavyDiversity, humanIgLambda, humanIgKappa))
# Generated via sort(rownames(seuratObj[["RNA"]])[grepl(rownames(seuratObj[["RNA"]]), pattern = '^KIR)]). 
humanKIRGenes <- c('KIR2DL1','KIR2DL3','KIR2DL4','KIR2DP1','KIR2DS4','KIR3DL1','KIR3DL2','KIR3DL3','KIR3DP1','KIR3DX1','KIRREL1','KIRREL1-IT1','KIRREL2','KIRREL3','KIRREL3-AS1','KIRREL3-AS2','KIRREL3-AS3')
.RegisterGeneSet("HumanKIRGenes", humanKIRGenes)
# Generated via sort(rownames(seuratObj[["RNA"]])[grepl(rownames(seuratObj[["RNA"]]), pattern = '^HLA-[^D]')]). 
humanHLAGenes <-  c('HLA-A','HLA-B','HLA-C','HLA-E','HLA-F','HLA-F-AS1','HLA-G','HLA-H','HLA-J','HLA-K','HLA-L','HLA-L.1','HLA-N','HLA-P','HLA-S','HLA-T','HLA-U','HLA-V','HLA-V.1','HLA-W','HLA-Z')
.RegisterGeneSet("HumanHLAGenes", humanHLAGenes)

# NOTE: this has been deprecated in favor of .2, which omits TCR constant region
.RegisterGeneSet('VariableGenes_Exclusion.1', unique(c(
	GetGeneSet('MMul10TcrGenes'),
	GetGeneSet('MMul10TcrConstantRegion'),
	GetGeneSet('MMul10_MHC'),
	GetGeneSet('MMul10_KIR'),
	GetGeneSet('MMul10_Ribosomal'),
	GetGeneSet('MMul10_Mitochondrial'),
	GetGeneSet('MMul10_Ig_Variable')
)))

.RegisterGeneSet('VariableGenes_Exclusion.2', unique(c(
	GetGeneSet('MMul10TcrGenes'),
	GetGeneSet('MMul10_MHC'),
	GetGeneSet('MMul10_KIR'),
	GetGeneSet('MMul10_Ribosomal'),
	GetGeneSet('MMul10_Mitochondrial'),
	GetGeneSet('MMul10_Ig_Variable')
)))

.RegisterGeneSet('HumanVariableGenes_Exclusion', unique(c(
  GetGeneSet('HumanMitochondrialGenes'), 
  GetGeneSet('HumanRibosomalGenes'), 
  GetGeneSet('HumanTcrGenes'), 
  GetGeneSet('HumanBcrGenes'), 
  GetGeneSet('HumanKIRGenes'), 
  GetGeneSet('HumanHLAGenes')
)))

.RegisterGeneSet("Glycolysis", c("ALDOA", "BPGM", "ENO1", "ENO2", "GAPDH", "HK1", "HK2", "HKDC1", "PFKL", "PGAM1", "PGAM2", "PGK1", "PKLR", "PKM", "TPI1"))
.RegisterGeneSet('Interferon_Response', c('IFI6','IFI27','MX1','ISG15','STAT1','LOC114672189','MX2','IFIT3'))
.RegisterGeneSet('Interferon_Response_IFI6_correlated', c('IFI6', 'ISG15', 'MX1', 'RNF213', 'BST2', 'DDX60', 'MX2', 'LOC100427967', 'STAT1', 'OAS2', 'DHX58', 'LY6E-1', 'IFIT3', 'SP100', 'EPSTI1'))

.RegisterGeneSet('MHC-II', c('IFI30', 'CD74', 'MAMU-DRB1', 'MAMU-DRA'))

# Lists provided by Rebecca Skalsky. NOTE: CD24, CD27 and CEACAM21 not annotated in MMul10
.RegisterGeneSet('B cells', c('MS4A1', 'CD79A', 'CD79B'))
.RegisterGeneSet('Follicular B cells', c('CD19', 'MS4A1', 'CR2', 'CD22', 'FCER2', 'CXCR5', 'MAMU-DRA', 'PAX5'))  # Also 'CD24'
.RegisterGeneSet('Activated B cell', c('CD19', 'CD69', 'MS4A1', 'IL2RA', 'CD86', 'FLT3', 'PAX5', 'BACH2', 'IRF4')) # Also CD27
.RegisterGeneSet('Germinal Center B cell', c('CD38', 'CD37', 'FCER2')) # Also CD27
.RegisterGeneSet('Dark zone (centroblast)', c('CD38', 'CXCR4', 'A4GALT', 'AICDA', 'FOXO1', 'MKI67', 'CDK1', 'CCND3', 'PAX5', 'BCL6', 'MEF2B'))
.RegisterGeneSet('Light zone (centrocyte)', c('CD38', 'CXCR5', 'CD83', 'BATF', 'BCL2A1'))
.RegisterGeneSet('Memory B cell precursor in LZ', c('CCR6', 'GPR183', 'PLAC8', 'MAML2', 'IL23A'))
.RegisterGeneSet('Memory B cells', c('CD19', 'MS4A1', 'CR2', 'CD40', 'CD80', 'CD86', 'MAMU-DRA', 'POU2AF1', 'PAX5')) # Also CEACAM21, CD27
.RegisterGeneSet('Plasma cell', c('SDC1', 'CD38', 'TNFRSF17', 'CXCR4', 'TNFRSF17', 'PRDM1', 'IRF4', 'XBP1')) # Also CD27
.RegisterGeneSet('Short-lived plasmablast', c('CD19', 'CD38', 'TNFRSF17', 'CXCR4')) # Also CD27
.RegisterGeneSet('Marginal Zone B cells (spleen)', c('CD19', 'MS4A1', 'CR2', 'FCER2', 'FCRL3', 'CD1C', 'EBF1', 'TCF3', 'PAX5')) # Also CD27

# SLC2A1 = Glut-1, SLC3A2 = CD98
.RegisterGeneSet('Redox', c('SLC2A1', 'SLC3A2'))

#' @title GetMMul10TcrGenes
#' @param includeConstantRegion If true, the MMul10TcrConstantRegion GeneSet will be included.
#' @description Returns a vector with MMul10 gene IDs (NCBI build) for TCR genes.
#' @export
#' @import Seurat
GetMMul10TcrGenes <- function(includeConstantRegion = FALSE){
	ret <- GetGeneSet('MMul10TcrGenes')
	if (includeConstantRegion) {
		ret <- unique(c(ret, GetGeneSet('MMul10TcrConstantRegion')))
	}

	return(ret)
}

#' @title GetMMul10IgGenes
#'
#' @description Returns a vector with MMul10 gene IDs (NCBI build) for TCR genes.
#' @export
#' @import Seurat
GetMMul10IgGenes <- function(){
	return(c(
		'LOC720839', #immunoglobulin heavy constant alpha 1-like
		'LOC708891', #immunoglobulin heavy constant gamma 1-like
		'LOC114679689', #immunoglobulin heavy constant epsilon-like
		'LOC114679691', #immunoglobulin heavy constant gamma 2-like
		'LOC114679690', #immunoglobulin heavy constant gamma 4-like
		'LOC710905', #immunoglobulin heavy constant gamma 4-like
		'LOC711872', # immunoglobulin heavy constant mu-like
		'LOC698810', #immunoglobulin iota chain
		'LOC114679695', #immunoglobulin mu heavy chain-like
		'LOC701504', # immunoglobulin kappa light chain
		'LOC107000555', # immunoglobulin lambda constant 6-like
		'LOC106996055', #immunoglobulin lambda-1 light chain-like
		'LOC708771', #immunoglobulin lambda constant 6-like
		'LOC106992418', #immunoglobulin lambda constant 6-like
		'LOC114670307', #immunoglobulin lambda constant 6-like
		'LOC106999340' # maybe immunoglobulin kappa light chain, though also labeled TRAV25
	))
}

#' @title ExpandGeneList
#'
#' @description Takes an input gene list and identifies any entries matching registered gene sets. Those will be expanded to the full gene list.
#' @param genes A vector of genes or gene set names
#' @param verbose Whether to log information about matches
#' @export
ExpandGeneList <- function(genes, verbose = TRUE) {
	genesMatchingSets <- genes[genes %in% names(pkg.env$GENE_SETS)]
	if (verbose && length(genesMatchingSets) > 0) {
		print(paste0('The following symbols match gene sets and will be expanded: ', paste0(genesMatchingSets, collapse = ',')))
	}

	ret <- genes[!genes %in% names(pkg.env$GENE_SETS)]
	for (geneSet in genesMatchingSets) {
		ret <- unique(c(ret, pkg.env$GENE_SETS[[geneSet]]))
	}

	return(ret)
}

.FindGeneSetsContaining <- function(genes){
  
  matchingGeneSets <- unlist(pkg.env$GENE_SETS)[ unlist(pkg.env$GENE_SETS) %in% genes]
  gene_matches_df <- data.frame(gene = matchingGeneSets, set = names(matchingGeneSets))
  #strip the trailing numbers from geneset names (and occasional peroids at the end of gene set signatures)
  gene_matches_df$set <- gsub("\\.?[0-9]*$", "", gene_matches_df$set) 
  
  gene_matches_df <- gene_matches_df |> 
    dplyr::group_by(gene) |>
    dplyr::mutate(geneset_union = paste(set, collapse = ", "))
  unique_gene_set_pairs <- unique.data.frame(gene_matches_df[,c("gene", "geneset_union")])
  
  return(unique_gene_set_pairs)
}

#' @title MakePhenotypingDotPlot
#'
#' @description Creates a DotPlot using custom gene sets and attempts to coarsely group gene sets by cell type. 
#' @param seuratObj A Seurat Object storing the count matrix to be used for phenotyping. 
#' @param yField The grouping variable used to calculate the average expression of genes and the y axis of the DotPlot.
#' @param scaled A boolean defining whether to color dots by scaled expression or unscaled expression.
#' @param gene_lists A vector of gene lists (defined by .RegisterGeneSet) to be queried and their genes be plotted.
#' @param scale.by Allow different scaling methods for dot size. 'radius' will de-emphasize lower/intermediately percent expressed genes.
#' @param assay Which assay to use in the input seuratObj
#' @export
MakePhenotypingDotPlot <- function(seuratObj,
                                   yField = 'ClusterNames_0.2',
                                   scaled = T, 
                                   gene_lists = c('Cytotoxicity', 'EffectorCytokines'), 
                                   assay = "RNA", 
                                   scale.by = "size"
){
  if (!is.logical(scaled)){
    stop("Please ensure scaled is either TRUE (to use per-gene scaled expression) or FALSE (for raw gene expression)")
  }
  #Parse gene_lists and coerce into a vector of genes to be plotted
  meta_gene_vector <- unique(unlist(sapply(gene_lists, FUN = GetGeneSet)))

  # drop trailing '-', in case the signature has a negative marker
  meta_gene_vector <- gsub(meta_gene_vector, pattern = '-$', replacement = '')
  
  #Get initial plotting and expression data from Seurat's version of the DotPlot
  plt <- Seurat::DotPlot(seuratObj, features = meta_gene_vector, group.by = yField, assay = assay)
  dotplot_df <- plt$data
  
  #Begin Phenotyping
  dotplot_df$Phenotype <- 0 #pass 0 as an easy debugging test for un/underannotated gene sets
  matching_genesets <- .FindGeneSetsContaining(dotplot_df$features.plot)
  dotplot_df <- merge(dotplot_df, matching_genesets, by.x = "features.plot", by.y ='gene')
  
  #Parse Phenotype Field to determine which cell type the phenotype targets
  dotplot_df$CellType <- "Unknown_CellType" #set default celltype as unknown
  dotplot_df[grepl("MemoryAndNaive", dotplot_df$geneset_union), "CellType"] <- "Cell Type: T Cells"
  dotplot_df[grepl("Myelocytes|Pro_Myelocytes", dotplot_df$geneset_union), "CellType"] <- "Cell Type: Neutrophil Precursors"
  dotplot_df[grepl("TandNK_Activation|Cytotoxicity|EffectorCytokines|Exhaustion", dotplot_df$geneset_union), "CellType"] <- "Cell Type: T Cells / NK Cells"
  
  #Sort the dataframe for faceting
  dotplot_df$CellType <- naturalsort::naturalfactor(dotplot_df$CellType)
  dotplot_df <- dotplot_df |> dplyr::arrange(CellType)
  
  #Determine if the data should be scaled or not, then alter the color scheme accordingly
  if(scaled){
    colorField <- sym("avg.exp.scaled")
    colors <- c("blue", "white", "red")
    colorLabel <- 'Scaled Average Expression'
  } else if(!scaled) {
    colorField <- sym("avg.exp")
    colors <- c("navy", "gold", "orange", "red")
    colorLabel <- 'Average Expression'
  }

  P1 <- ggplot(dotplot_df, aes(x = features.plot, y = id, size = pct.exp, color = !!colorField)) + 
    geom_point() + 
    facet_wrap(~CellType+geneset_union, scales = "free_x") + 
    egg::theme_article() + 
    scale_color_gradientn(colors = colors) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    ylab(yField) + 
    xlab("Genes") +
    labs(color=colorLabel, size = "Percentage of Cells\nWith Gene Expression")
  
 if (scale.by == 'size'){
    P1 <- P1 + ggplot2::scale_size()
  } else if (scale.by == 'radius'){
    P1 <- P1 + ggplot2::scale_radius()
  } else {
    stop("Please specify scale.by = 'size' or scale.by = 'radius'")
  }
 
  return(P1)
}


#' @title ListGeneSets
#'
#' @description Prints a list of the gene sets registered in this package
#' @export
ListGeneSets <- function(){
	print(sort(names(pkg.env$GENE_SETS)))
}