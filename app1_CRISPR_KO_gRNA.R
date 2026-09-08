suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicFeatures)
  library(GenomicRanges)
})

GENOME <- BSgenome.Hsapiens.UCSC.hg38
TXDB   <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ── Helper functions ──────────────────────────────────────────────────────────

gc_pct <- function(s) {
  round(100 * nchar(gsub("[^GCgc]", "", s)) / nchar(s), 1)
}

has_poly_t <- function(s, n = 4) {
  grepl(sprintf("T{%d,}", n), s, ignore.case = TRUE)
}

# On-target efficiency score (0–100).
# Captures key predictors from Doench 2014 and CRISPRscan literature:
#   GC content, position-specific nucleotide preferences, homopolymer penalties.
on_target_score <- function(seq) {
  score <- 0

  # GC content: optimal 40-70% (40 pts)
  gc <- gc_pct(seq) / 100
  score <- score + if (gc >= 0.40 && gc <= 0.70) 40
                   else if ((gc >= 0.30 && gc < 0.40) || (gc > 0.70 && gc <= 0.80)) 20
                   else 0

  # Position 20 (PAM-proximal): G/C preferred (15 pts)
  score <- score + switch(substr(seq, 20, 20), G = 15, C = 10, T = 5, A = 5, 0)

  # Position 1 (5' end): G preferred for Pol-III initiation (10 pts)
  score <- score + switch(substr(seq, 1, 1), G = 10, A = 7, T = 5, C = 3, 0)

  # Seed region GC (positions 9-20): optimal 30-70% (20 pts)
  sgc <- gc_pct(substr(seq, 9, 20)) / 100
  score <- score + if (sgc >= 0.30 && sgc <= 0.70) 20
                   else if ((sgc >= 0.20 && sgc < 0.30) || (sgc > 0.70 && sgc <= 0.80)) 10
                   else 0

  # Penalty: any homopolymer >=4 (-15); extra -10 for poly-T (Pol-III terminator)
  if (grepl("A{4,}|C{4,}|G{4,}|T{4,}", seq))  score <- score - 15
  if (grepl("T{4,}", seq, ignore.case = TRUE))  score <- score - 10

  round(min(100, max(0, score)), 1)
}

# Off-target risk score (0–100; lower = more specific).
# Heuristic based on seed-region composition and guide complexity.
off_target_risk <- function(seq) {
  risk <- 0

  # High GC in PAM-proximal seed (positions 9-20) → more stable mismatches
  risk <- risk + round(gc_pct(substr(seq, 9, 20)) / 100 * 40, 1)

  # Low sequence complexity → matches repetitive genome regions
  n_unique <- length(unique(strsplit(seq, "")[[1]]))
  risk <- risk + (4 - n_unique) * 10

  # Dinucleotide repeats
  if      (grepl("(..)\\1{2,}", seq, perl = TRUE)) risk <- risk + 20
  else if (grepl("(..)\\1",     seq, perl = TRUE)) risk <- risk + 10

  # G-runs → G-quadruplex potential
  if      (grepl("GGGG", seq)) risk <- risk + 10
  else if (grepl("GGG",  seq)) risk <- risk + 5

  round(min(100, max(0, risk)), 1)
}

# Fetch merged exons for a gene, numbered 5'→3' along the gene.
# Returns data.frame: exon_num, chr, start, end, width, strand
get_gene_exons <- function(eid) {
  exon_list <- exonsBy(TXDB, by = "gene")
  key <- as.character(unname(eid))   # unname() strips the named-vector attribute from mapIds
  if (!key %in% names(exon_list)) return(NULL)

  ex_gr <- exon_list[[key]]
  # runValue() extracts Rle run values — safer than as.character(strand(...)[1])
  gene_strand <- as.character(runValue(strand(ex_gr)))[1]
  if (is.na(gene_strand)) gene_strand <- "+"

  ex_red <- reduce(ex_gr)

  # Build data frame directly from GRanges accessors to avoid as.data.frame/rename issues
  ex_df <- data.frame(
    chr    = as.character(seqnames(ex_red)),
    start  = as.integer(start(ex_red)),
    end    = as.integer(end(ex_red)),
    width  = as.integer(width(ex_red)),
    strand = gene_strand,
    stringsAsFactors = FALSE
  )

  ex_df <- if (gene_strand == "+") arrange(ex_df, start) else arrange(ex_df, desc(start))
  ex_df$exon_num <- seq_len(nrow(ex_df))
  ex_df
}

# Scan an exon sequence for NGG-PAM SpCas9 guides (both strands).
design_ko_guides <- function(seq_str, chr, win_start, exon_num, glen) {
  slen <- nchar(seq_str)
  rows <- list()
  add  <- function(...) rows[[length(rows) + 1]] <<- list(...)

  # + strand: look for [ACGT]GG (NGG PAM)
  pm <- gregexpr("[ACGT]GG", seq_str, perl = TRUE)[[1]]
  if (pm[1] > 0) {
    for (p in pm) {
      ge <- p - 1;  gs <- ge - glen + 1
      if (gs < 1) next
      sq <- substr(seq_str, gs, ge)
      if (grepl("[^ACGT]", sq)) next
      add(guide_sequence = sq,
          pam            = substr(seq_str, p, p + 2),
          guide_strand   = "+",
          chromosome     = chr,
          position       = win_start + gs - 1,
          exon           = exon_num)
    }
  }

  # - strand: CC[ACGT] on + strand encodes NGG on - strand
  pm <- gregexpr("CC[ACGT]", seq_str, perl = TRUE)[[1]]
  if (pm[1] > 0) {
    for (p in pm) {
      gs <- p + 3;  ge <- gs + glen - 1
      if (ge > slen) next
      sq_p <- substr(seq_str, gs, ge)
      if (grepl("[^ACGT]", sq_p)) next
      add(guide_sequence = as.character(reverseComplement(DNAString(sq_p))),
          pam            = as.character(reverseComplement(DNAString(substr(seq_str, p, p + 2)))),
          guide_strand   = "-",
          chromosome     = chr,
          position       = win_start + ge - 1,
          exon           = exon_num)
    }
  }

  if (!length(rows)) return(NULL)
  bind_rows(lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- page_fluid(
  theme = bs_theme(bootswatch = "flatly"),

  tags$head(tags$style(HTML("
    .gene-card {
      background: #f0fff4;
      border-left: 4px solid #1a7a4a;
      padding: 10px 16px;
      border-radius: 4px;
      margin-bottom: 14px;
      font-size: 0.95em;
    }
    h4 { margin-bottom: 0.3rem; }
  "))),

  titlePanel("CRISPR/Cas9 Knockout gRNA Designer (hg38 / SpCas9 NGG)"),

  layout_sidebar(
    sidebar = sidebar(
      width = 310,

      h5("Gene Target"),
      textInput("gene", NULL, placeholder = "Gene symbol, e.g. TP53, KRAS"),
      actionButton("load_gene", "Load Gene & Exons", class = "btn-outline-primary w-100"),

      uiOutput("exon_selector_ui"),

      hr(),
      h5("Guide Parameters"),
      numericInput("glen", "Guide Length (nt)", 20, 17, 24),

      hr(),
      h5("Quality Filters"),
      sliderInput("gc_range",   "GC Content (%)",       0, 100, c(30, 80), step = 5),
      checkboxInput("rm_polyt", "Remove poly-T (>=4T)", TRUE),

      hr(),
      h5("Scoring Filters"),
      helpText("On-target score: higher is better. Off-target risk: lower is safer."),
      sliderInput("min_on",  "Min On-target Score",  0, 100, 40, step = 5),
      sliderInput("max_off", "Max Off-target Risk",  0, 100, 70, step = 5),

      hr(),
      actionButton("go", "Design gRNAs", class = "btn-primary w-100"),
      br(), br(),
      downloadButton("dl_csv", "Download CSV", class = "btn-outline-secondary w-100")
    ),

    uiOutput("gene_card"),

    tabsetPanel(
      tabPanel("Guide Table",
        br(),
        DTOutput("tbl")
      ),
      tabPanel("Score Distribution",
        br(),
        plotOutput("score_plot", height = "480px")
      ),
      tabPanel("Exon Map",
        br(),
        plotOutput("exon_plot", height = "320px")
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  rv_guides <- reactiveVal(NULL)
  rv_gene   <- reactiveVal(NULL)
  rv_exons  <- reactiveVal(NULL)

  # ── Step 1: load gene → populate exon selector ──
  observeEvent(input$load_gene, {
    req(nzchar(trimws(input$gene)))
    rv_gene(NULL); rv_exons(NULL); rv_guides(NULL)

    withProgress(message = "Loading gene info...", value = 0, {
      tryCatch({
        sym <- toupper(trimws(input$gene))

        setProgress(0.2, detail = "Entrez ID lookup")
        eid <- unname(mapIds(org.Hs.eg.db, keys = sym, column = "ENTREZID",
                             keytype = "SYMBOL", multiVals = "first"))
        validate(need(!is.na(eid), paste0("Gene not found: '", sym,
          "'. Use an official HGNC symbol (e.g. TP53, KRAS).")))

        setProgress(0.5, detail = "Fetching exon coordinates")
        exons <- get_gene_exons(eid)
        validate(need(!is.null(exons) && nrow(exons) > 0,
          paste0("No annotated exons for ", sym, " in hg38 knownGene.")))

        rv_exons(exons)
        rv_gene(list(
          symbol      = sym,
          chr         = exons$chr[1],
          gene_strand = exons$strand[1],
          n_exons     = nrow(exons)
        ))

        choices <- setNames(
          as.character(exons$exon_num),
          paste0("Exon ", exons$exon_num,
                 "  (", format(exons$start, big.mark = ",", scientific = FALSE),
                 "–", format(exons$end, big.mark = ",", scientific = FALSE),
                 ", ", exons$width, " bp)")
        )
        updateCheckboxGroupInput(session, "exon_select",
          choices  = choices,
          selected = as.character(exons$exon_num[seq_len(min(3, nrow(exons)))])
        )

      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 12)
      })
    })
  })

  # ── Exon selector UI ──
  output$exon_selector_ui <- renderUI({
    req(rv_gene())
    gi <- rv_gene()
    tagList(
      hr(),
      h5("Select Exon(s) to Target"),
      helpText(paste0(gi$n_exons, " merged exons found. ",
        "For KO, early exons shared across isoforms are preferred.")),
      checkboxGroupInput("exon_select", NULL, choices = character(0))
    )
  })

  # ── Step 2: design gRNAs ──
  observeEvent(input$go, {
    req(rv_gene(), rv_exons(), length(input$exon_select) > 0)
    rv_guides(NULL)

    withProgress(message = "Designing gRNAs...", value = 0, {
      tryCatch({
        exons     <- rv_exons()
        sel_nums  <- as.integer(input$exon_select)
        sel_exons <- filter(exons, exon_num %in% sel_nums)

        all_res <- list()
        for (i in seq_len(nrow(sel_exons))) {
          ex  <- sel_exons[i, ]
          if (!ex$chr %in% seqnames(GENOME)) next
          setProgress(i / nrow(sel_exons) * 0.65,
            detail = paste0("Scanning exon ", ex$exon_num))
          seq_str <- as.character(
            getSeq(GENOME, GRanges(ex$chr, IRanges(ex$start, ex$end))))
          res <- design_ko_guides(seq_str, ex$chr, ex$start, ex$exon_num, input$glen)
          if (!is.null(res)) all_res[[i]] <- res
        }

        setProgress(0.70, detail = "Scoring guides")
        validate(need(length(all_res) > 0,
          "No gRNAs found in selected exons."))

        guides <- bind_rows(all_res) %>%
          distinct(guide_sequence, guide_strand, .keep_all = TRUE) %>%
          mutate(
            gc_content      = vapply(guide_sequence, gc_pct,          numeric(1)),
            poly_t_flag     = vapply(guide_sequence, has_poly_t,       logical(1)),
            on_target_score = vapply(guide_sequence, on_target_score,  numeric(1)),
            off_target_risk = vapply(guide_sequence, off_target_risk,  numeric(1))
          ) %>%
          filter(
            gc_content      >= input$gc_range[1],
            gc_content      <= input$gc_range[2],
            on_target_score >= input$min_on,
            off_target_risk <= input$max_off
          )

        if (input$rm_polyt) guides <- filter(guides, !poly_t_flag)
        guides <- arrange(guides, exon, desc(on_target_score))
        rv_guides(guides)

      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 12)
      })
    })
  })

  # ── Gene card ──
  output$gene_card <- renderUI({
    req(rv_gene())
    gi <- rv_gene()
    n  <- if (!is.null(rv_guides())) nrow(rv_guides()) else "—"
    div(class = "gene-card",
      tags$b(gi$symbol), " | ",
      tags$b("Chr: "), gi$chr, " | ",
      tags$b("Strand: "), gi$gene_strand, " | ",
      tags$b("Exons (merged): "), gi$n_exons, " | ",
      tags$b("gRNAs passing filters: "), n
    )
  })

  # ── Guide table ──
  output$tbl <- renderDT({
    req(rv_guides())
    df <- rv_guides() %>%
      transmute(
        `Exon`                   = exon,
        `Guide Sequence (5'→3')` = guide_sequence,
        `PAM`                    = pam,
        `Strand`                 = guide_strand,
        `Chr`                    = chromosome,
        `Position (hg38)`        = format(position, big.mark = ",", scientific = FALSE),
        `GC (%)`                 = gc_content,
        `On-target Score`        = on_target_score,
        `Off-target Risk`        = off_target_risk
      )

    datatable(
      df, rownames = FALSE, filter = "top", selection = "multiple",
      options = list(pageLength = 25, scrollX = TRUE, dom = "lftip",
        columnDefs = list(list(className = "dt-center", targets = 0:8)))
    ) %>%
      formatStyle(
        "On-target Score",
        backgroundColor = styleInterval(c(40, 70),
          c("#fee2e2", "#fef3c7", "#dcfce7"))
      ) %>%
      formatStyle(
        "Off-target Risk",
        backgroundColor = styleInterval(c(30, 60),
          c("#dcfce7", "#fef3c7", "#fee2e2"))
      ) %>%
      formatStyle(
        "GC (%)",
        backgroundColor = styleInterval(c(29.9, 39.9, 70.1, 80.1),
          c("#fee2e2", "#fee2e2", "white", "#fee2e2", "#fee2e2"))
      ) %>%
      formatStyle(
        "Strand",
        color = styleEqual(c("+", "-"), c("#2166ac", "#d73027"))
      )
  }, server = TRUE)

  # ── Score distribution plot ──
  output$score_plot <- renderPlot({
    req(rv_guides())
    df  <- rv_guides()
    pal <- RColorBrewer::brewer.pal(max(3, n_distinct(df$exon)), "Set2")

    p1 <- ggplot(df, aes(x = on_target_score, fill = factor(exon))) +
      geom_histogram(binwidth = 5, color = "white", alpha = 0.85, position = "stack") +
      scale_fill_brewer(palette = "Set2", name = "Exon") +
      labs(title = "On-target Efficiency Score",
           x = "Score (higher = better)", y = "Count") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), legend.position = "top")

    p2 <- ggplot(df, aes(x = off_target_risk, fill = factor(exon))) +
      geom_histogram(binwidth = 5, color = "white", alpha = 0.85, position = "stack") +
      scale_fill_brewer(palette = "Set2", name = "Exon") +
      labs(title = "Off-target Risk Score",
           x = "Score (lower = safer)", y = "Count") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), legend.position = "top")

    p3 <- ggplot(df, aes(x = on_target_score, y = off_target_risk,
                          color = factor(exon), shape = guide_strand)) +
      geom_point(alpha = 0.75, size = 3) +
      scale_color_brewer(palette = "Set2", name = "Exon") +
      scale_shape_manual(values = c("+" = 16, "-" = 17), name = "Guide strand") +
      geom_vline(xintercept = 70, linetype = "dashed", color = "#1a7a4a", linewidth = 0.5) +
      geom_hline(yintercept = 30, linetype = "dashed", color = "#1a7a4a", linewidth = 0.5) +
      annotate("text", x = 72, y = 95, label = "High efficiency\nLow risk →",
               hjust = 0, size = 3, color = "#1a7a4a") +
      labs(title = "Efficiency vs Off-target Risk",
           x = "On-target Score", y = "Off-target Risk") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), legend.position = "top")

    grid.arrange(p1, p2, p3, ncol = 2,
      layout_matrix = rbind(c(1, 2), c(3, 3)))
  })

  # ── Exon map ──
  output$exon_plot <- renderPlot({
    req(rv_exons(), rv_gene())
    exons <- rv_exons()
    gi    <- rv_gene()
    sel   <- as.integer(isolate(input$exon_select))
    exons$selected <- exons$exon_num %in% sel

    ggplot(exons,
      aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = selected)) +
      geom_rect(color = "black", linewidth = 0.6, alpha = 0.85) +
      geom_text(aes(x = (start + end) / 2, y = 0.5,
                    label = paste0("E", exon_num)),
                size = 3.5, fontface = "bold") +
      scale_fill_manual(
        values = c(`TRUE` = "#4ade80", `FALSE` = "#cbd5e1"),
        labels = c(`TRUE` = "Selected", `FALSE` = "Not selected"),
        name   = NULL
      ) +
      scale_x_continuous(
        labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
      labs(
        title    = paste0(gi$symbol, " — Merged Exon Map (",
                          gi$chr, ", ", gi$gene_strand, " strand)"),
        subtitle = "Green = selected for gRNA design",
        x = "Genomic Position (hg38)", y = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        axis.text.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.line.y   = element_blank(),
        legend.position = "top",
        plot.title    = element_text(face = "bold")
      )
  })

  # ── Download ──
  output$dl_csv <- downloadHandler(
    filename = function() {
      paste0(toupper(input$gene), "_Cas9KO_gRNAs_", Sys.Date(), ".csv")
    },
    content = function(f) {
      write.csv(
        rv_guides() %>% select(exon, guide_sequence, pam, guide_strand,
          chromosome, position, gc_content, on_target_score, off_target_risk),
        f, row.names = FALSE
      )
    }
  )
}

shinyApp(ui, server)
