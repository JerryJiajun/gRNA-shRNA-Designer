# gRNA/shRNA Designer

A suite of R/Shiny apps for designing CRISPR and RNAi reagents from a gene symbol or genomic coordinate, using the **hg38** reference genome.

**[View the screenshot demo / GitHub Pages site →](https://JerryJiajun.github.io/gRNA-shRNA-Designer/)**

## Apps

| App | File | Description |
|---|---|---|
| CRISPR/Cas9 Knockout gRNA Designer | [`app1_CRISPR_KO_gRNA.R`](app1_CRISPR_KO_gRNA.R) | Pulls a gene's merged exons from hg38 and scans selected exons for SpCas9 guides, scoring each for on-target efficiency and off-target risk to surface the strongest frameshift-KO candidates. |
| CRISPRi Knockdown gRNA Designer | [`app2_CRISPRi_KD_gRNA.R`](app2_CRISPRi_KD_gRNA.R) | Locates a gene's TSS in hg38, scans the promoter window on both strands, and reports signed distance-to-TSS for each guide — flagging the −50 to +300 bp zone where CRISPRi silencing is most effective. |
| CRISPR SNP Knock-in Designer | [`app4_CRISPR_KO_SNP_knockin.R`](app4_CRISPR_KO_SNP_knockin.R) | Verifies a SNP's reference allele against hg38, ranks nearby guides by cut-site distance, and auto-builds an ssODN repair template with an optional PAM-disrupting mutation. |
| Gene Knockdown shRNA Designer | [`app3_gene_KD_shRNA.R`](app3_gene_KD_shRNA.R) | Scans every annotated transcript of a gene for candidate siRNA target sites, scores each by Reynolds/Schwarz-Khvorova design rules, and returns a ready-to-clone U6/H1 hairpin. |

## Running the apps

Each app is a standalone R/Shiny script. Install dependencies with `install_packages1_2_3.R`, then open the desired `appN_*.R` file in RStudio and click **Run App**, or run it from the command line:

```r
shiny::runApp("app1_CRISPR_KO_gRNA.R")
```

## Demo site

`docs/` contains a static, screenshot-based walkthrough of the apps, published via GitHub Pages (no live R server required). See [docs/index.md](docs/index.md).
