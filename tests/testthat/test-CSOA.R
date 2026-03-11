test_that("basicHeatmap returns a ggplot object", {
    mat <- matrix(0, 10, 20)
    mat[sample(length(mat), 50)] <- runif(50, max = 2.5)
    p <- basicHeatmap(mat)
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
})


test_that("getPairs works", {
    v <- c('ASD', 'VBN', 'HJKL')
    expect_identical(getPairs(v), list(c('ASD', 'VBN'),
                                       c('ASD', 'HJKL'),
                                       c('VBN', 'HJKL')))
})

test_that("geneRadialPlot returns a ggplot object", {
    edgesDF <- data.frame(gene1 = paste0('G', c(1, 2, 3, 4, 7, 8, 10,
                                                11, 11, 10, 10, 10)),
                          gene2 = paste0('G', c(2, 5, 1, 8, 4, 9, 12,
                                                13, 14, 13, 16, 14)))
    edgesDF <- henna::connectedComponents(edgesDF, 'group')
    p <- geneRadialPlot(edgesDF, 'component', extraCircles=1)
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
})

test_that("overlapNetworkPlot returns a ggraph object", {
    df <- data.frame(gene1 = paste0('G', c(1, 2, 5, 6, 7, 17)),
                     gene2 = paste0('G', c(2, 5, 8, 11, 11, 11)),
                     rank = c(1, 1, 3, 3, 3, 3))
    p <- overlapNetworkPlot(df)
    expect_equal(is(p), 'ggraph')
})

test_that("overlapCutoffPlot returns a gg object", {
    overlapDF <- data.frame(gene1=paste0('G', c(1, 3, 7, 6, 8, 2, 4, 3, 4, 5)),
                            gene2=paste0('G', c(2, 7, 2, 5, 4, 5, 1, 2, 2, 8)),
                            rank=c(1, 2, 3, 4, 4, 6, 7, 7, 7, 10))
    p <- overlapCutoffPlot(overlapDF)
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
})

test_that("percentileSets works", {
    mat <- matrix(0, 5, 200)
    expect_error(percentileSets(mat))

    rownames(mat) <- paste0('G', seq(5))
    expect_error(percentileSets(mat))

    colnames(mat) <- paste0('C', seq(200))
    mat[1, c(3, 5, seq(11, 20))] <- c(3, 7, rep(1, 10))
    mat[3, c(4, 7, 8, seq(11, 30))] <- c(4, 5, 8, rep(1, 20))
    mat[5, c(2, 8, 9, seq(11, 20))] <- c(1, 5, 5, rep(1, 10))
    expect_warning(percentileSets(mat))

    mat[2, c(3, 10, seq(11, 20))] <- c(2, 2, rep(1, 10))
    mat[4, c(4, 14, 15, 20, seq(21, 70))] <- c(3, 3, 3, 3, rep(1, 50))
    res <- list(c('C3', 'C5'),
                c('C3', 'C10'),
                c('C4', 'C7', 'C8'),
                c('C4', 'C14', 'C15', 'C20'),
                c('C8', 'C9'))
    names(res) <- rownames(mat)
    expect_identical(percentileSets(mat), res)
})

test_that("overlapGenes works", {
    df <- data.frame(gene1 = paste0('G', c(1, 2, 7, 3, 4, 5, 1)),
                     gene2 = paste0('G', c(2, 3, 5, 4, 7, 6, 4)))
    expect_identical(overlapGenes(df), paste0('G', c(1, 2, 7, 3, 4, 5, 6)))
})

test_that("overlapPairs works", {
    df <- data.frame(gene1 = paste0('G', c(1, 2, 7, 3, 4, 5, 1)),
                     gene2 = paste0('G', c(2, 3, 5, 4, 7, 6, 4)))
    expect_identical(overlapPairs(df), list(c('G1', 'G2'),
                                            c('G2', 'G3'),
                                            c('G7', 'G5'),
                                            c('G3', 'G4'),
                                            c('G4', 'G7'),
                                            c('G5', 'G6'),
                                            c('G1', 'G4')))
})

test_that("runCSOA works", {
    sceObj <- scRNAseq::BaronPancreasData('human')
    sceObj <- scuttle::logNormCounts(sceObj)
    acinarMarkers <- c('PRSS1', 'KLK1', 'CTRC', 'PNLIP', 'AKR1C3', 'CTRB1',
                       'DUOXA2', 'ALDOB', 'REG3A', 'SERPINA3', 'PRSS3', 'REG1B',
                       'CFB',	'GDF15',	'MUC1','ANPEP', 'ANGPTL4', 'OLFM4',
                       'GSTA1', 'LGALS2', 'PDZK1IP1', 'RARRES2', 'CXCL17',
                       'UBD', 'GSTA2', 'LYZ', 'RBPJL', 'PTF1A', 'CELA3A',
                       'SPINK1', 'ZG16', 'CEL', 'CELA2A', 'CPB1', 'CELA1',
                       'PNLIPRP1', 'RNASE1', 'AMY2B', 'CPA2','CPA1', 'CELA3B',
                       'CTRB2', 'PLA2G1B', 'PRSS2', 'CLPS', 'REG1A', 'SYCN'
    )
    sceObj <- runCSOA(sceObj, list(CSOA_acinar = acinarMarkers))
    expect_equal(mean(sceObj$CSOA_acinar[sceObj$label == 'acinar']), 0.4988258,
                 tolerance=0.001)

    sceObj <- runCSOA(sceObj, list(CSOA_acinar = acinarMarkers),
                      jaccardCutoff=1/3)
    expect_equal(mean(sceObj$CSOA_acinar[sceObj$label == 'acinar']), 0.4979294,
                 tolerance=0.001)

    sceObj <- runCSOA(sceObj, list(CSOA_acinar = acinarMarkers),
                      osMethod='minmax')
    expect_equal(mean(sceObj$CSOA_acinar[sceObj$label == 'acinar']), 0.4899093,
                 tolerance=0.001)

    sceObj <- runCSOA(sceObj, list(CSOA_acinar = acinarMarkers),
                      percentile=80)
    expect_equal(mean(sceObj$CSOA_acinar[sceObj$label == 'acinar']), 0.5220412,
                 tolerance=0.001)

    expect_warning(sceObj <- runCSOA(sceObj,
                                     list(CSOA_null =
                                              rownames(sceObj)[c(3,
                                                                 17,
                                                                 210,
                                                                 333,
                                                                 422)])))
    sceObj <- runCSOA(sceObj, list(CSOA_acinar = acinarMarkers),
                      pairFileTemplate='pairs')
    pairDF <- qGrab('pairsCSOA_acinar.qs2')
    expect_equal(pairDF$pairScore[[1]], 1.643546,
                 tolerance=0.001)

    alphaMarkers <- c('GCG', 'TTR', 'PCSK2', 'FXYD5', 'LDB2', 'MAFB',
                      'CHGA', 'SCGB2A1', 'GLS', 'FAP', 'DPP4', 'GPR119',
                      'PAX6', 'NEUROD1', 'LOXL4', 'PLCE1', 'GC', 'KLHL41',
                      'FEV', 'PTGER3', 'RFX6', 'SMARCA1', 'PGR', 'IRX1',
                      'UCP2', 'RGS4', 'KCNK16', 'GLP1R', 'ARX', 'POU3F4',
                      'RESP18', 'PYY', 'SLC38A5', 'TM4SF4', 'CRYBA2', 'SH3GL2',
                      'PCSK1', 'PRRG2', 'IRX2', 'ALDH1A1','PEMT', 'SMIM24',
                      'F10', 'SCGN', 'SLC30A8')

    ductalMarkers <- c('CFTR', 'SERPINA5', 'SLPI', 'TFF1', 'CFB', 'LGALS4',
                       'CTSH',	'PERP', 'PDLIM3',	'WFDC2', 'SLC3A1', 'AQP1',
                       'ALDH1A3', 'VTCN1',	'KRT19', 'TFF2', 'KRT7', 'CLDN4',
                       'LAMB3', 'TACSTD2', 'CCL2', 'DCDC2','CXCL2', 'CLDN10',
                       'HNF1B', 'KRT20', 'MUC1', 'ONECUT1', 'AMBP', 'HHEX',
                       'ANXA4', 'SPP1', 'PDX1', 'SERPINA3', 'GDF15', 'AKR1C3',
                       'MMP7', 'DEFB1', 'SERPING1', 'TSPAN8', 'CLDN1', 'S100A10',
                       'PIGR')

    expect_warning(sceObj <- runCSOA(sceObj, list(CSOA_alpha = alphaMarkers,
                                   CSOA_ductal = ductalMarkers)))
    expect_equal(mean(sceObj$CSOA_alpha[sceObj$label == 'alpha']), 0.2896088,
                 tolerance=0.001)

    expect_equal(mean(sceObj$CSOA_ductal[sceObj$label == 'ductal']), 0.1955368,
                 tolerance=0.001)

    expect_warning(sceObj <- runCSOA(sceObj, list(CSOA_alpha = alphaMarkers,
                                   CSOA_ductal = ductalMarkers), percentile=80))
    expect_equal(mean(sceObj$CSOA_alpha[sceObj$label == 'alpha']), 0.3677292,
                 tolerance=0.001)

    expect_equal(mean(sceObj$CSOA_ductal[sceObj$label == 'ductal']), 0.1976859,
                 tolerance=0.001)
})

test_that("scoreModules works", {
    mat <- matrix(0, 500, 300)
    rownames(mat) <- paste0('G', seq(500))
    colnames(mat) <- paste0('C', seq(300))
    mat[sample(8000)] <- runif(8000, max=13)
    genes1 <- paste0('G', seq(100))
    mat[genes1, 20:50] <- matrix(runif(100 * 31, min = 14, max = 15),
    nrow = 100, ncol = 31)
    genes2 <- paste0('G', seq(101, 200))
    mat[genes2, 70:100] <- matrix(runif(100 * 31, min = 14, max = 15),
    nrow = 100, ncol = 31)
    genes <- union(genes1, genes2)
    mat <- mat[genes, ]
    overlapDF <- generateOverlaps(mat)
    overlapDF <- processOverlaps(overlapDF)
    overlapDF <- connectedComponents(overlapDF)
    components <- unique(overlapDF$component)
    df <- scoreModules(mat, overlapDF, components)[[2]]

    expect_equal(nrow(df), 300)
    expect_equal(ncol(df), length(components))
})
