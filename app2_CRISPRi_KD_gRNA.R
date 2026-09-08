suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(dplyr)
  library(ggplot2)
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

# Find all NGG-PAM gRNAs in a sequence window.
# Position reported is the 5' end of the guide RNA in genome coordinates.
# Distance to TSS uses gene-direction sign convention:
#   negative = upstream of TSS (promoter), positive = downstream (gene body).
design_guides <- function(seq_str, chr, win_start, tss, gene_strand, glen) {
  slen <- nchar(seq_str)
  rows <- list()

  add <- function(...) rows[[length(rows) + 1]] <<- list(...)

  dist_to_tss <- function(pos) {
    if (gene_strand == "+") pos - tss else tss - pos
  }

  # ── + strand guides: look for [ACGT]GG (NGG PAM) ──
  pm <- gregexpr("[ACGT]GG", seq_str, perl = TRUE)[[1]]
  if (pm[1] > 0) {
    for (p in pm) {
      ge <- p - 1
      gs <- ge - glen + 1
      if (gs < 1) next
      sq <- substr(seq_str, gs, ge)
      if (grepl("[^ACGT]", sq)) next

      g_pos <- win_start + gs - 1   # 5' end of guide on + strand

      add(guide_sequence  = sq,
          pam             = substr(seq_str, p, p + 2),
          guide_strand    = "+",
          chromosome      = chr,
          position        = g_pos,
          distance_to_TSS = dist_to_tss(g_pos))
    }
  }

  # ── - strand guides: look for CC[ACGT] on + strand (= NGG on - strand) ──
  pm <- gregexpr("CC[ACGT]", seq_str, perl = TRUE)[[1]]
  if (pm[1] > 0) {
    for (p in pm) {
      gs <- p + 3
      ge <- gs + glen - 1
      if (ge > slen) next
      sq_p <- substr(seq_str, gs, ge)
      if (grepl("[^ACGT]", sq_p)) next

      sq  <- as.character(reverseComplement(DNAString(sq_p)))
      pam <- as.character(reverseComplement(DNAString(substr(seq_str, p, p + 2))))

      # 5' end of the guide on the - strand = highest genomic coordinate
      g_pos <- win_start + ge - 1

      add(guide_sequence  = sq,
          pam             = pam,
          guide_strand    = "-",
          chromosome      = chr,
          position        = g_pos,
          distance_to_TSS = dist_to_tss(g_pos))
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
      background: #f0f7ff;
      border-left: 4px solid #2166ac;
      padding: 10px 16px;
      border-radius: 4px;
      margin-bottom: 14px;
      font-size: 0.95em;
    }
    .badge-up   { background:#2166ac; color:white; padding:2px 8px; border-radius:10px; }
    .badge-down { background:#d73027; color:white; padding:2px 8px; border-radius:10px; }
    h4 { margin-bottom: 0.3rem; }
  "))),

  titlePanel("CRISPRi gRNA Designer (hg38 / SpCas9 NGG)"),

  layout_sidebar(
    sidebar = sidebar(
      width = 290,

      h5("Gene Target"),
      textInput("gene", NULL, placeholder = "Gene symbol, e.g. MYC, TP53"),

      hr(),
      h5("Search Window Around TSS"),
      fluidRow(
        column(6, numericInput("upstream",   "Upstream (bp)",   300, 0, 2000, 50)),
        column(6, numericInput("downstream", "Downstream (bp)", 300, 0, 2000, 50))
      ),

      hr(),
      h5("Guide Parameters"),
      numericInput("glen", "Guide Length (nt)", 20, 17, 24),

      hr(),
      h5("CRISPRi Window Filter"),
      helpText("CRISPRi is most effective in the -50 to +300 bp window."),
      checkboxInput("use_win", "Restrict to window", TRUE),
      conditionalPanel("input.use_win",
        fluidRow(
          column(6, numericInput("win_lo", "Min (bp)", -50)),
          column(6, numericInput("win_hi", "Max (bp)",  300))
        )
      ),

      hr(),
      h5("Quality Filters"),
      sliderInput("gc_range", "GC Content (%)", 0, 100, c(30, 80), step = 5),
      checkboxInput("rm_polyt", "Remove poly-T runs (>=4T)", TRUE),

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
      tabPanel("Position Plot",
        br(),
        plotOutput("pos_plot", height = "420px")
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  rv_guides <- reactiveVal(NULL)
  rv_gene   <- reactiveVal(NULL)

  observeEvent(input$go, {
    req(nzchar(trimws(input$gene)))
    rv_guides(NULL)
    rv_gene(NULL)

    withProgress(message = "Designing gRNAs...", value = 0, {
      tryCatch({
        sym <- toupper(trimws(input$gene))

        setProgress(0.1, detail = "Fetching gene annotation")
        eid <- mapIds(org.Hs.eg.db, sym, "ENTREZID", "SYMBOL", multiVals = "first")
        validate(need(!is.na(eid), paste0("Gene symbol not found: '", sym,
          "'. Please use an official HGNC symbol (e.g. MYC, TP53).")))

        all_txs <- transcriptsBy(TXDB, by = "gene")
        validate(need(as.character(eid) %in% names(all_txs),
          paste0("No annotated transcripts for ", sym, " in hg38 knownGene.")))

        txs    <- all_txs[[as.character(eid)]]

        # promoters(upstream=0, downstream=1) yields a 1-bp range at the TSS.
        # For both strands, start == end == TSS coordinate.
        tss_gr <- suppressWarnings(promoters(txs, upstream = 0, downstream = 1))
        tss_df <- as.data.frame(tss_gr) %>%
          mutate(tss_coord = start) %>%
          distinct(seqnames, strand, tss_coord)

        rv_gene(list(
          symbol      = sym,
          chr         = as.character(tss_df$seqnames[1]),
          gene_strand = as.character(tss_df$strand[1]),
          tss         = tss_df$tss_coord[1],
          n_tss       = nrow(tss_df)
        ))

        setProgress(0.3, detail = "Fetching genome sequence & scanning PAMs")
        all_res <- list()

        for (i in seq_len(nrow(tss_df))) {
          chr   <- as.character(tss_df$seqnames[i])
          gstr  <- as.character(tss_df$strand[i])
          tss_i <- tss_df$tss_coord[i]

          if (!chr %in% seqnames(GENOME)) next

          ws <- max(1, tss_i - input$upstream)
          we <- min(seqlengths(GENOME)[chr], tss_i + input$downstream)

          seq_str <- as.character(getSeq(GENOME, GRanges(chr, IRanges(ws, we))))
          res     <- design_guides(seq_str, chr, ws, tss_i, gstr, input$glen)

          if (!is.null(res)) {
            res$tss_position <- tss_i
            all_res[[i]]     <- res
          }
        }

        setProgress(0.85, detail = "Applying quality filters")
        validate(need(length(all_res) > 0,
          "No gRNAs found — try expanding the upstream/downstream search window."))

        guides <- bind_rows(all_res) %>%
          distinct(guide_sequence, guide_strand, .keep_all = TRUE) %>%
          mutate(
            gc_content  = vapply(guide_sequence, gc_pct,    numeric(1)),
            poly_t_flag = vapply(guide_sequence, has_poly_t, logical(1))
          ) %>%
          filter(gc_content >= input$gc_range[1],
                 gc_content <= input$gc_range[2])

        if (input$rm_polyt) guides <- filter(guides, !poly_t_flag)
        if (input$use_win)  guides <- filter(guides,
                                              distance_to_TSS >= input$win_lo,
                                              distance_to_TSS <= input$win_hi)

        guides <- arrange(guides, distance_to_TSS)
        rv_guides(guides)

      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 12)
      })
    })
  })

  # ── Gene info card ──
  output$gene_card <- renderUI({
    req(rv_gene())
    gi <- rv_gene()
    n  <- if (!is.null(rv_guides())) nrow(rv_guides()) else "calculating..."

    div(class = "gene-card",
      tags$b(gi$symbol), " | ",
      tags$b("Chr: "), gi$chr, " | ",
      tags$b("Gene strand: "), gi$gene_strand, " | ",
      tags$b("TSS: "), format(gi$tss, big.mark = ",", scientific = FALSE), " (hg38) | ",
      tags$b("Unique TSS: "), gi$n_tss, " | ",
      tags$b("gRNAs passing filters: "), n
    )
  })

  # ── Guide table ──
  output$tbl <- renderDT({
    req(rv_guides())

    df <- rv_guides() %>%
      transmute(
        `Guide Sequence (5'→3')`   = guide_sequence,
        `PAM`                      = pam,
        `Guide Strand`             = guide_strand,
        `Chromosome`               = chromosome,
        `Position (hg38, 5' end)`  = format(position, big.mark = ",", scientific = FALSE),
        `Distance to TSS (bp)`     = distance_to_TSS,
        `GC Content (%)`           = gc_content
      )

    datatable(
      df,
      rownames  = FALSE,
      filter    = "top",
      selection = "multiple",
      options   = list(
        pageLength = 25,
        scrollX    = TRUE,
        dom        = "lftip",
        columnDefs = list(list(className = "dt-center", targets = 1:6))
      )
    ) %>%
      formatStyle(
        "Distance to TSS (bp)",
        color      = styleInterval(c(-0.5), c("#2166ac", "#d73027")),
        fontWeight = "bold"
      ) %>%
      formatStyle(
        "GC Content (%)",
        backgroundColor = styleInterval(
          c(29.9, 39.9, 70.1, 80.1),
          c("#fde0d0", "#fde0d0", "white", "#fde0d0", "#fde0d0")
        )
      ) %>%
      formatStyle(
        "Guide Strand",
        color = styleEqual(c("+", "-"), c("#2166ac", "#d73027"))
      )
  }, server = TRUE)

  # ── Position plot ──
  output$pos_plot <- renderPlot({
    req(rv_guides(), rv_gene())
    gi  <- rv_gene()
    df  <- rv_guides()

    ggplot(df, aes(x = distance_to_TSS, y = 0, color = guide_strand)) +
      # Optimal CRISPRi window shading
      annotate("rect",
        xmin = -50, xmax = 300,
        ymin = -Inf, ymax = Inf,
        fill = "#c6dbef", alpha = 0.35
      ) +
      annotate("text",
        x = 125, y = Inf, vjust = -0.4,
        label = "Optimal CRISPRi window", size = 3.5, color = "#2166ac"
      ) +
      # TSS line
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.9) +
      annotate("text",
        x = 0, y = Inf, vjust = -0.4, hjust = -0.1,
        label = "TSS", fontface = "bold", color = "black", size = 4
      ) +
      # Guide points
      geom_jitter(size = 3, alpha = 0.75, height = 0.25, width = 0) +
      scale_color_manual(
        values = c("+" = "#2166ac", "-" = "#d73027"),
        name   = "Guide strand"
      ) +
      scale_x_continuous(
        name   = "Distance to TSS (bp)",
        labels = function(x) paste0(ifelse(x > 0, "+", ""), x)
      ) +
      labs(
        title    = paste0(gi$symbol, " — gRNA positions relative to TSS"),
        subtitle = paste0("Chr: ", gi$chr, "  |  Gene strand: ", gi$gene_strand,
                          "  |  TSS: ", format(gi$tss, big.mark = ",")),
        y        = NULL
      ) +
      theme_classic(base_size = 13) +
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
      paste0(toupper(input$gene), "_CRISPRi_gRNAs_", Sys.Date(), ".csv")
    },
    content = function(f) {
      out <- rv_guides() %>%
        select(guide_sequence, pam, guide_strand, chromosome, position,
               distance_to_TSS, gc_content, tss_position)
      write.csv(out, f, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
