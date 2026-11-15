# R packages for the native CellChat column of Benchmark 1 (tools/cellchat/native_r.py runs run_cellchat.R).
# Nothing else in the project needs R. The versions the results were produced with are in R_sessionInfo.txt
# (R 4.3.2, CellChat 2.1.2, Matrix 1.6-1.1).
#
#   Rscript tools/cellchat/install_r_packages.R
options(repos = c(CRAN = "https://cloud.r-project.org"))
for (p in c("remotes", "BiocManager", "Matrix", "NMF", "circlize")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
# CellChat's Bioconductor dependencies are not installed by install_github
BiocManager::install(c("Biobase", "BiocGenerics", "BiocNeighbors", "ComplexHeatmap"), update = FALSE, ask = FALSE)
remotes::install_github("jinworks/CellChat@v2.1.2", upgrade = "never")
stopifnot(as.character(packageVersion("CellChat")) == "2.1.2")
