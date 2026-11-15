# CellChat 2.1.2 itself, in R: the tutorial pipeline up to computeCommunProb, on the inputs that
# native_r.py writes. Usage: Rscript run_cellchat.R <dir> <nboot> <seed>
#
# <dir>/genes.txt, cells.txt, labels.txt, levels.txt, species.txt and the genes x cells normalised matrix as raw
# little-endian binaries (i.bin int32, p.bin int32, x.bin float64), so R reads exactly the doubles
# Python holds. Writes prob.bin and pval.bin (float64, R's column-major order, K x K x nLR),
# LRsig.txt (rownames of LRsig, in R's order), features_sig.txt and session.txt.
suppressPackageStartupMessages({library(CellChat); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; nboot <- as.integer(args[2]); seed <- as.integer(args[3])
rd <- function(f) readLines(file.path(dir, f))
bin <- function(f, what, n) readBin(file.path(dir, f), what = what, size = if (what == "double") 8 else 4,
                                    n = n, endian = "little")
genes <- rd("genes.txt"); cells <- rd("cells.txt"); labels <- rd("labels.txt"); levels <- rd("levels.txt")
species <- if (file.exists(file.path(dir, "species.txt"))) rd("species.txt")[1] else "human"
p <- bin("p.bin", "integer", length(cells) + 1)
X <- sparseMatrix(i = bin("i.bin", "integer", p[length(p)]), p = p, x = bin("x.bin", "double", p[length(p)]),
                  dims = c(length(genes), length(cells)), dimnames = list(genes, cells), index1 = FALSE)
meta <- data.frame(labels = factor(labels, levels = levels), row.names = cells)

t0 <- Sys.time()
future::plan("sequential")
cellchat <- createCellChat(object = X, meta = meta, group.by = "labels")
cellchat@DB <- subsetDB(if (species == "mouse") CellChatDB.mouse else CellChatDB.human)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean", nboot = nboot, seed.use = seed)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

net <- cellchat@net
stopifnot(identical(dimnames(net$prob)[[2]], levels))
writeBin(as.vector(net$prob), file.path(dir, "prob.bin"), size = 8, endian = "little")
writeBin(as.vector(net$pval), file.path(dir, "pval.bin"), size = 8, endian = "little")
writeLines(rownames(cellchat@LR$LRsig), file.path(dir, "LRsig.txt"))
writeLines(cellchat@var.features$features, file.path(dir, "features_sig.txt"))
writeLines(c(R.version.string, paste("CellChat", as.character(packageVersion("CellChat"))),
             paste("presto", as.character(packageVersion("presto"))), paste("nboot", nboot), paste("seed", seed),
             paste("seconds", round(elapsed, 1)), paste("dim prob", paste(dim(net$prob), collapse = " x "))),
           file.path(dir, "session.txt"))
