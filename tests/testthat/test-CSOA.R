test_that("connectedComponents works", {
    df <- data.frame(gene1 = paste0('G', c(1, 2, 12, 3, 4, 4, 7, 8, 12, 11)),
                     gene2 = paste0('G', c(2, 3, 6, 3, 7, 8, 8, 9, 10, 12)))
    expect_equal(max(connectedComponents(df)$component), 3)

    df <- data.frame(gene1 = paste0('G', c(1, 15, 4, 1, 1, 7, 7, 22, 12, 7)),
                     gene2 = paste0('G', c(2, 7, 6, 15, 2, 4, 5, 7, 7, 6)))
    expect_equal(max(connectedComponents(df)$component), 1)

    df <- data.frame(gene1 = paste0('G', c(1, 3, 5, 7, 9)),
                     gene2 = paste0('G', c(2, 4, 6, 8, 10)))
    expect_equal(max(connectedComponents(df)$component), 5)
})

test_that("findRankCutoff works", {
    freqDF <- data.frame(rank = c(1, 2, 4, 8),
                         n = c(1, 3, 4, 2))
    expect_equal(findRankCutoff(freqDF), 4)

    freqDF <- data.frame(rank = c(1, 2, 5, 8, 11, 14),
                         n = c(1, 3, 3, 3, 3, 2))
    expect_equal(findRankCutoff(freqDF), 6.5)

    freqDF <- data.frame(apples = c(1, 2, 5, 8, 11, 14),
                         oranges = c(1, 3, 3, 3, 3, 2))
    expect_error(findRankCutoff(freqDF))

    freqDF <- data.frame(rank = c(),
                         n = c())
    expect_error(findRankCutoff(freqDF))
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
