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

# ── Scoring ───────────────────────────────────────────────────────────────────
# Criteria synthesised from Reynolds 2004 (Nat Biotech) and Schwarz/Khvorova 2003
# (thermodynamic asymmetry rule).  Score is 0–10.
#
# The sense strand = same sequence as the mRNA target (5'→3').
# The antisense strand is the functional (guide) strand loaded into RISC.
# Thermodynamic asymmetry: 5' end of antisense (= 3' end of sense) should be
# A/U-rich so that RISC preferentially loads the antisense strand.

score_sirna <- function(seq) {
  s  <- toupper(seq)
  n  <- nchar(s)
  ch <- strsplit(s, "")[[1]]

  pts <- 0L

  # [1] GC content  (Reynolds criterion 1: optimal 36–52%)
  gc  <- round(100 * sum(ch %in% c("G", "C")) / n, 1)
  gc_pts <- if (gc >= 36 && gc <= 52) 3L else if (gc >= 30 && gc <= 60) 1L else 0L
  pts <- pts + gc_pts

  # [2] A/U-rich at 3' half of sense (positions 13–19 for 19-mer)
  #     = 5' end of antisense (guide strand) → low stability → favours RISC loading
  end_pos  <- max(1L, n - 6L):n
  au_end   <- sum(ch[end_pos] %in% c("A", "T"))
  asym_pts <- if (au_end >= 5L) 3L else if (au_end >= 3L) 2L else if (au_end >= 1L) 1L else 0L
  pts <- pts + asym_pts

  # [3] A/U at position n (very 5' of antisense) – extra bonus
  end1_pts <- if (ch[n] %in% c("A", "T")) 1L else 0L
  pts <- pts + end1_pts

  # [4] G/C at position 1 (5' of sense) – stabilises sense 5' end, discourages
  #     sense-strand RISC loading
  gc1_pts <- if (ch[1] %in% c("G", "C")) 1L else 0L
  pts <- pts + gc1_pts

  # [5] A/U at central position (near RISC cleavage site)
  mid     <- ceiling(n / 2L)
  mid_pts <- if (ch[mid] %in% c("A", "T")) 1L else 0L
  pts <- pts + mid_pts

  # [6] No homopolymer run ≥ 4 nt  (TTTT = Pol-III termination; GGGG = G-quadruplex)
  nopoly_pts <- if (!grepl("(.)\\1{3,}", s, perl = TRUE)) 1L else 0L
  pts <- pts + nopoly_pts

  # Theoretical max: 3 + 3 + 1 + 1 + 1 + 1 = 10
  list(total     = pts,
       gc        = gc,
       gc_pts    = gc_pts,
       asym_pts  = asym_pts,
       end1_pts  = end1_pts,
       gc1_pts   = gc1_pts,
       mid_pts   = mid_pts,
       nopoly    = nopoly_pts == 1L)
}

# ── shRNA hairpin builder ──────────────────────────────────────────────────────
# Returns the DNA sequence to clone into a U6/H1 Pol-III vector.
# Structure: [sense] – [loop] – [antisense] – TTTTTT (terminator)
# A 'G' is prepended if the sense does not start with G (improves U6 transcription).

make_shrna <- function(sense, loop = "TTCAAGAGA") {
  antisense  <- as.character(reverseComplement(DNAString(sense)))
  prefix     <- if (!startsWith(sense, "G")) "G" else ""
  full       <- paste0(prefix, sense, loop, antisense, "TTTTTT")
  list(sense = sense, antisense = antisense, loop = loop,
       full = full, g_prepended = nchar(prefix) > 0L)
}

# ── Transcript scanner ────────────────────────────────────────────────────────
scan_transcript <- function(tx_seq, tx_name, target_len = 19L) {
  s  <- toupper(as.character(tx_seq))
  tx_len <- nchar(s)
  if (tx_len < target_len) return(NULL)

  rows <- vector("list", tx_len - target_len + 1L)
  k    <- 0L

  for (i in seq_len(tx_len - target_len + 1L)) {
    sq <- substr(s, i, i + target_len - 1L)
    if (grepl("[^ACGT]", sq)) next
    sc <- score_sirna(sq)
    k  <- k + 1L
    rows[[k]] <- data.frame(
      target_sequence = sq,
      tx_position     = i,
      transcript_id   = tx_name,
      tx_length       = tx_len,
      gc_content      = sc$gc,
      score           = sc$total,
      poly_t          = grepl("TTTT", sq, ignore.case = TRUE),
      stringsAsFactors = FALSE
    )
  }

  if (k == 0L) return(NULL)
  bind_rows(rows[seq_len(k)])
}

# ── UI ────────────────────────────────────────────────────────────────────────

LOOP_CHOICES <- c(
  "TTCAAGAGA  (Brummelkamp 2002, most common)" = "TTCAAGAGA",
  "AAGTTCTCT  (reverse complement of above)"   = "AAGTTCTCT",
  "CTCGAG     (short / XhoI site)"             = "CTCGAG",
  "CCACCG     (minimal stem-loop)"             = "CCACCG"
)

ui <- page_fluid(
  theme = bs_theme(bootswatch = "flatly"),

  tags$head(tags$style(HTML("
    .gene-card {
      background:#f0fff4; border-left:4px solid #276749;
      padding:10px 16px; border-radius:4px; margin-bottom:14px; font-size:0.95em;
    }
    .seq-box  { font-family:monospace; font-size:0.85em; word-break:break-all;
                background:#f8f8f8; border:1px solid #ddd;
                padding:4px 8px; border-radius:3px; }
    .score-badge { font-weight:bold; padding:2px 8px; border-radius:10px; }
  "))),

  titlePanel("shRNA Designer (hg38 / U6–H1 vector, SpCas9-independent)"),

  layout_sidebar(
    sidebar = sidebar(
      width = 290,

      h5("Gene Target"),
      textInput("gene", NULL, placeholder = "Gene symbol, e.g. MYC, GAPDH"),

      hr(),
      h5("shRNA Parameters"),
      numericInput("target_len", "Target / Stem Length (nt)", 19L, 17L, 21L),
      selectInput("loop", "Hairpin Loop Sequence", choices = LOOP_CHOICES),

      hr(),
      h5("Transcript Selection"),
      checkboxInput("canonical_only", "Longest transcript only", FALSE),
      helpText("By default all transcripts are scanned. Guides are prioritised by",
               "how many transcripts they target."),

      hr(),
      h5("Filters"),
      sliderInput("min_score", "Min Efficacy Score (0–10)", 0, 10, 5, step = 1),
      sliderInput("gc_range",  "GC Content (%)",            0, 100, c(30, 65), step = 5),
      checkboxInput("rm_polyt", "Remove guides with poly-T (≥4 T's)", TRUE),

      hr(),
      numericInput("n_top", "Max results to show", 50L, 5L, 500L, 5L),

      hr(),
      actionButton("go", "Design shRNAs", class = "btn-success w-100"),
      br(), br(),
      downloadButton("dl_csv", "Download CSV", class = "btn-outline-secondary w-100")
    ),

    uiOutput("gene_card"),

    tabsetPanel(
      id = "tabs",
      tabPanel("Results Table",
        br(),
        p(tags$b("Ranking:"), " guides sorted first by number of transcripts targeted",
          "(broader knockdown), then by efficacy score."),
        DTOutput("tbl")
      ),
      tabPanel("Score Distribution",
        br(),
        plotOutput("score_plot", height = "380px")
      ),
      tabPanel("shRNA Sequences",
        br(),
        helpText("Select rows in the Results Table to display full hairpin sequences here."),
        uiOutput("shrna_cards")
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  rv_results <- reactiveVal(NULL)
  rv_gene    <- reactiveVal(NULL)

  observeEvent(input$go, {
    req(nzchar(trimws(input$gene)))
    rv_results(NULL); rv_gene(NULL)

    withProgress(message = "Designing shRNAs...", value = 0, {
      tryCatch({
        sym <- toupper(trimws(input$gene))

        # ── 1. Gene lookup ──
        setProgress(0.05, detail = "Fetching gene annotation")
        eid <- mapIds(org.Hs.eg.db, sym, "ENTREZID", "SYMBOL", multiVals = "first")
        validate(need(!is.na(eid), paste0("Gene symbol not found: '", sym,
          "'. Use an official HGNC symbol (e.g. MYC, TP53, GAPDH).")))

        all_txs <- transcriptsBy(TXDB, by = "gene")
        validate(need(as.character(eid) %in% names(all_txs),
          paste0("No annotated transcripts for '", sym, "' in hg38 knownGene.")))

        txs   <- all_txs[[as.character(eid)]]
        tx_df <- as.data.frame(txs) %>% arrange(desc(width))

        if (input$canonical_only) {
          tx_df <- slice_head(tx_df, n = 1L)
          txs   <- txs[txs$tx_name %in% tx_df$tx_name]
        }

        rv_gene(list(
          symbol = sym,
          chr    = as.character(tx_df$seqnames[1]),
          strand = as.character(tx_df$strand[1]),
          n_tx   = nrow(tx_df)
        ))

        # ── 2. Extract spliced mRNA sequences ──
        setProgress(0.15, detail = "Extracting transcript sequences")
        exons_by_tx <- exonsBy(TXDB, by = "tx", use.names = TRUE)
        tx_names    <- tx_df$tx_name
        avail       <- intersect(tx_names, names(exons_by_tx))
        validate(need(length(avail) > 0,
          "No exon annotations found for this gene's transcripts."))

        # Cap at 20 transcripts to keep runtime reasonable
        if (length(avail) > 20L) avail <- avail[seq_len(20L)]

        setProgress(0.25, detail = "Extracting spliced sequences from hg38")
        # Build spliced mRNA sequences by concatenating strand-corrected exon sequences.
        # exonsBy(...) orders exons by exon_rank (5'→3' transcript order), so
        # getSeq() on each exon (which handles strand automatically) can be directly
        # pasted to produce the mature mRNA sequence — no extractTxSeqs() needed.
        tx_seqs <- DNAStringSet(setNames(
          vapply(avail, function(tx) {
            ex <- exons_by_tx[[tx]]
            paste(as.character(getSeq(GENOME, ex)), collapse = "")
          }, character(1L)),
          avail
        ))

        # ── 3. Scan all transcripts ──
        setProgress(0.35, detail = "Scanning for siRNA candidates")
        all_hits <- list()
        for (i in seq_along(tx_seqs)) {
          res <- scan_transcript(tx_seqs[[i]], names(tx_seqs)[i], input$target_len)
          if (!is.null(res)) all_hits[[i]] <- res
          setProgress(0.35 + 0.40 * i / length(tx_seqs))
        }

        validate(need(length(all_hits) > 0, "No siRNA candidates found."))

        # ── 4. Aggregate, rank, filter ──
        setProgress(0.80, detail = "Ranking and filtering")
        raw <- bind_rows(all_hits)

        # For each unique target sequence: count transcripts hit and find best position
        tx_cov <- raw %>%
          distinct(target_sequence, transcript_id) %>%
          group_by(target_sequence) %>%
          summarise(
            transcripts_hit = n(),
            transcript_ids  = paste(sort(unique(transcript_id)), collapse = "; "),
            .groups = "drop"
          )

        results <- raw %>%
          group_by(target_sequence) %>%
          # Keep the row with the highest score; break ties by earliest position
          arrange(desc(score), tx_position) %>%
          slice_head(n = 1L) %>%
          ungroup() %>%
          left_join(tx_cov, by = "target_sequence") %>%
          filter(score       >= input$min_score,
                 gc_content  >= input$gc_range[1],
                 gc_content  <= input$gc_range[2]) %>%
          { if (input$rm_polyt) filter(., !poly_t) else . } %>%
          arrange(desc(transcripts_hit), desc(score)) %>%
          slice_head(n = input$n_top) %>%
          mutate(
            antisense = vapply(target_sequence,
              function(s) as.character(reverseComplement(DNAString(s))), character(1)),
            shrna_full = vapply(target_sequence,
              function(s) make_shrna(s, input$loop)$full, character(1)),
            g_prepended = !startsWith(target_sequence, "G")
          )

        rv_results(results)

      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 12)
      })
    })
  })

  # ── Gene info card ──
  output$gene_card <- renderUI({
    req(rv_gene())
    gi <- rv_gene()
    n  <- if (!is.null(rv_results())) nrow(rv_results()) else "calculating..."

    div(class = "gene-card",
      tags$b(gi$symbol), " | ",
      tags$b("Chr: "), gi$chr, " | ",
      tags$b("Strand: "), gi$strand, " | ",
      tags$b("Transcripts scanned: "), gi$n_tx, " | ",
      tags$b("shRNAs passing filters: "), n
    )
  })

  # ── Results table ──
  output$tbl <- renderDT({
    req(rv_results())
    gi <- rv_gene()

    df <- rv_results() %>%
      transmute(
        `Target Sequence (5'→3')`  = target_sequence,
        `Antisense (5'→3')`        = antisense,
        `Transcripts Hit`          = paste0(transcripts_hit, " / ", gi$n_tx),
        `mRNA Position`            = tx_position,
        `GC Content (%)`           = gc_content,
        `Efficacy Score (/10)`     = score,
        `Poly-T`                   = ifelse(poly_t, "Yes", ""),
        `U6 G-prefix needed`       = ifelse(g_prepended, "Yes", "")
      )

    datatable(df,
              rownames  = FALSE,
              filter    = "top",
              selection = "multiple",
              options   = list(
                pageLength = 20L,
                scrollX    = TRUE,
                dom        = "lftip",
                columnDefs = list(list(className = "dt-center", targets = 2:7))
              )) %>%
      formatStyle("Efficacy Score (/10)",
        background = styleInterval(c(4.9, 6.9, 8.9),
          c("#fde0d0", "#fff3c8", "#d4edda", "#b8daff")),
        fontWeight = "bold") %>%
      formatStyle("GC Content (%)",
        backgroundColor = styleInterval(
          c(29.9, 35.9, 52.1, 65.1),
          c("#fde0d0", "#fff3c8", "white", "#fff3c8", "#fde0d0"))) %>%
      formatStyle("Poly-T", color = styleEqual("Yes", "#d73027")) %>%
      formatStyle("Transcripts Hit",
        fontWeight = styleEqual(
          paste0(seq_len(gi$n_tx), " / ", gi$n_tx),
          rep("bold", gi$n_tx)))
  }, server = TRUE)

  # ── Score distribution plot ──
  output$score_plot <- renderPlot({
    req(rv_results())
    df <- rv_results()

    ggplot(df, aes(x = score)) +
      geom_histogram(binwidth = 1, fill = "#276749", color = "white", alpha = 0.8,
                     boundary = 0) +
      geom_vline(xintercept = input$min_score - 0.5,
                 linetype = "dashed", color = "#d73027", linewidth = 0.9) +
      annotate("text",
        x = input$min_score + 0.1, y = Inf, hjust = 0, vjust = 1.5,
        label = paste0("Min score = ", input$min_score),
        color = "#d73027", size = 4) +
      scale_x_continuous("Predicted Efficacy Score (0–10)",
                         breaks = 0:10, limits = c(-0.5, 10.5)) +
      labs(title    = paste0(rv_gene()$symbol, " — shRNA Efficacy Score Distribution"),
           subtitle = "Score criteria: GC content, thermodynamic asymmetry, cleavage-site A/U, no homopolymers",
           y        = "Number of candidate guides") +
      theme_classic(base_size = 13) +
      theme(plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(size = 10, color = "grey40"))
  })

  # ── shRNA sequence cards ──
  output$shrna_cards <- renderUI({
    req(rv_results())
    sel <- input$tbl_rows_selected
    if (!length(sel)) {
      return(div(class = "alert alert-info",
        "Select one or more rows in the Results Table to view full shRNA sequences."))
    }

    df <- rv_results()[sel, , drop = FALSE]

    tagList(lapply(seq_len(nrow(df)), function(i) {
      row  <- df[i, ]
      sh   <- make_shrna(row$target_sequence, input$loop)
      loop <- sh$loop

      div(
        style = "margin-bottom:18px; border:1px solid #ccc; border-radius:6px; padding:16px;",

        tags$h5(paste0("Candidate #", sel[i],
          "  |  Score: ", row$score, "/10",
          "  |  Transcripts: ", row$transcripts_hit, "/", rv_gene()$n_tx,
          "  |  GC: ", row$gc_content, "%")),

        tags$table(
          style = "width:100%; border-collapse:collapse; font-size:0.9em;",

          tags$tr(
            tags$td(tags$b("Sense (target):"), style = "width:200px; padding:4px 0;"),
            tags$td(div(class = "seq-box", style = "color:#2166ac", sh$sense))
          ),
          tags$tr(
            tags$td(tags$b("Antisense (guide):"), style = "padding:4px 0;"),
            tags$td(div(class = "seq-box", style = "color:#d73027", sh$antisense))
          ),
          tags$tr(
            tags$td(tags$b("Loop:"), style = "padding:4px 0;"),
            tags$td(div(class = "seq-box", style = "color:#555", loop))
          ),
          tags$tr(
            tags$td(tags$b("Full shRNA (clone):"), style = "padding:4px 0;"),
            tags$td(
              div(class = "seq-box",
                if (sh$g_prepended)
                  tags$span(style = "color:#888; font-style:italic;", "G"),
                tags$span(style = "color:#2166ac", sh$sense),
                tags$span(style = "color:#666; font-weight:bold;", loop),
                tags$span(style = "color:#d73027", sh$antisense),
                tags$span(style = "color:#aaa", "TTTTTT")
              )
            )
          ),
          tags$tr(
            tags$td(tags$b("mRNA position:"), style = "padding:4px 0;"),
            tags$td(paste0(row$tx_position, "–",
                           row$tx_position + nchar(row$target_sequence) - 1L,
                           "  (transcript: ", row$transcript_id, ")"))
          )
        ),

        if (sh$g_prepended)
          div(class = "alert alert-warning",
              style = "padding:6px 10px; margin-top:8px; font-size:0.85em;",
              tags$b("Note:"), " a G was prepended (shown in grey) because the target does not",
              " start with G. This improves U6 promoter transcription efficiency.")
      )
    }))
  })

  # ── Download ──
  output$dl_csv <- downloadHandler(
    filename = function() paste0(toupper(input$gene), "_shRNA_", Sys.Date(), ".csv"),
    content  = function(f) {
      out <- rv_results() %>%
        select(target_sequence, antisense, shrna_full,
               transcript_ids, transcripts_hit, transcript_id, tx_position, tx_length,
               gc_content, score, poly_t, g_prepended)
      write.csv(out, f, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
