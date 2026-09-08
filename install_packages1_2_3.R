if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_pkgs <- c("shiny", "bslib", "DT", "dplyr", "ggplot2")
bioc_pkgs <- c(
  "Biostrings",
  "BSgenome.Hsapiens.UCSC.hg38",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "GenomicFeatures",
  "GenomicRanges"
)

install.packages(setdiff(cran_pkgs, rownames(installed.packages())))
BiocManager::install(setdiff(bioc_pkgs, rownames(installed.packages())), ask = FALSE)

cat("All packages installed. Run the app with: shiny::runApp('app.R')\n")
