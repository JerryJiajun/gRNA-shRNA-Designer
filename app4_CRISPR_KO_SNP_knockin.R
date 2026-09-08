suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(dplyr)
  library(ggplot2)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
})

GENOME <- BSgenome.Hsapiens.UCSC.hg38

# ── Helpers ───────────────────────────────────────────────────────────────────

gc_pct <- function(s) round(100 * nchar(gsub("[^GCgc]", "", s)) / nchar(s), 1)

# Find all NGG-PAM guides within a window around a SNP.
# cut_site convention (both strands cut at same locus for SpCas9):
#   + strand guide, PAM at P:  cut_site = P - 3
#   - strand guide, CCN at Q:  cut_site = Q + 5  (≈ 3 nt into protospacer from PAM side)
find_guides <- function(seq_str, chr, win_start, snp_pos, ref_len, glen = 20L) {
  slen    <- nchar(seq_str)
  snp_end <- snp_pos + ref_len - 1L
  rows    <- list()

  add <- function(...) rows[[length(rows) + 1L]] <<- list(...)

  # + strand
  pm <- gregexpr("[ACGT]GG", seq_str, perl = TRUE)[[1]]
  if (pm[1] > 0L) for (p in pm) {
    ge <- p - 1L; gs <- ge - glen + 1L
    if (gs < 1L) next
    sq <- substr(seq_str, gs, ge)
    if (grepl("[^ACGT]", sq)) next

    g_start   <- win_start + gs - 1L
    g_end     <- win_start + ge - 1L
    pam_start <- win_start + p  - 1L
    cut_site  <- pam_start - 3L

    add(guide_sequence  = sq,
        pam             = substr(seq_str, p, p + 2L),
        guide_strand    = "+",
        chromosome      = chr,
        guide_start     = g_start,
        guide_end       = g_end,
        pam_start       = pam_start,
        pam_mut_pos     = pam_start + 2L,   # change 3rd G → A to break NGG
        pam_mut_from    = "G",
        pam_mut_to      = "A",
        cut_site        = cut_site,
        dist_to_snp     = cut_site - snp_pos,
        snp_in_guide    = snp_pos <= g_end   && snp_end >= g_start,
        snp_in_pam      = snp_pos <= pam_start + 2L && snp_end >= pam_start,
        snp_in_seed     = snp_pos <= g_end   && snp_end >= (g_end - 11L),
        disrupts_recut  = (snp_pos <= g_end && snp_end >= g_start) ||
                          (snp_pos <= pam_start + 2L && snp_end >= pam_start))
  }

  # - strand (CCN on + strand = NGG on - strand)
  pm <- gregexpr("CC[ACGT]", seq_str, perl = TRUE)[[1]]
  if (pm[1] > 0L) for (p in pm) {
    gs <- p + 3L; ge <- gs + glen - 1L
    if (ge > slen) next
    sq_p <- substr(seq_str, gs, ge)
    if (grepl("[^ACGT]", sq_p)) next
    sq  <- as.character(reverseComplement(DNAString(sq_p)))
    pam <- as.character(reverseComplement(DNAString(substr(seq_str, p, p + 2L))))

    g_start   <- win_start + gs - 1L
    g_end     <- win_start + ge - 1L
    pam_start <- win_start + p  - 1L
    cut_site  <- pam_start + 5L

    add(guide_sequence  = sq,
        pam             = pam,
        guide_strand    = "-",
        chromosome      = chr,
        guide_start     = g_start,
        guide_end       = g_end,
        pam_start       = pam_start,
        pam_mut_pos     = pam_start + 1L,   # change 2nd C → T  (CCN→CTN; - strand NGG→NAG)
        pam_mut_from    = "C",
        pam_mut_to      = "T",
        cut_site        = cut_site,
        dist_to_snp     = cut_site - snp_pos,
        snp_in_guide    = snp_pos <= g_end   && snp_end >= g_start,
        snp_in_pam      = snp_pos <= pam_start + 2L && snp_end >= pam_start,
        snp_in_seed     = snp_pos >= g_start && snp_end <= (g_start + 11L),
        disrupts_recut  = (snp_pos <= g_end && snp_end >= g_start) ||
                          (snp_pos <= pam_start + 2L && snp_end >= pam_start))
  }

  if (!length(rows)) return(NULL)
  bind_rows(lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}

# Build ssODN repair template.
# Structure: [left homology arm] + [ALT allele] + [right homology arm]
# Richardson 2016: use the strand whose 3' end is closer to the cut site.
# Optional: introduce a single-nt PAM-disrupting mutation in the arms.
design_ssODN <- function(chr, snp_pos, ref, alt, arm_len, guide,
                          mutate_pam = FALSE) {
  ref_len   <- nchar(ref)
  alt_len   <- nchar(alt)
  chr_len   <- seqlengths(GENOME)[chr]

  left_start  <- max(1L,       snp_pos - arm_len)
  left_end    <- snp_pos - 1L
  right_start <- snp_pos + ref_len
  right_end   <- min(chr_len,  snp_pos + ref_len - 1L + arm_len)

  left_arm  <- if (left_end  >= left_start)
    toupper(as.character(getSeq(GENOME, GRanges(chr, IRanges(left_start,  left_end))))) else ""
  right_arm <- if (right_end >= right_start)
    toupper(as.character(getSeq(GENOME, GRanges(chr, IRanges(right_start, right_end))))) else ""

  sense    <- paste0(left_arm, toupper(alt), right_arm)
  pam_note <- ""

  # Optional PAM mutation — mutate one base in the arm to break the PAM
  if (mutate_pam && !is.null(guide)) {
    mp   <- guide$pam_mut_pos
    mfr  <- guide$pam_mut_from
    mto  <- guide$pam_mut_to

    if (mp >= left_start && mp <= right_end && mp != snp_pos) {
      tidx <- if (mp < snp_pos)
        mp - left_start + 1L
      else
        nchar(left_arm) + alt_len + (mp - right_start) + 1L

      if (tidx >= 1L && tidx <= nchar(sense) &&
          toupper(substr(sense, tidx, tidx)) == toupper(mfr)) {
        substr(sense, tidx, tidx) <- mto
        pam_note <- paste0("PAM mutation: template pos ", tidx,
                           " (genomic ", chr, ":", mp, ")  ", mfr, " → ", mto)
      }
    }
  }

  antisense <- as.character(reverseComplement(DNAString(sense)))
  L         <- nchar(sense)

  # Richardson 2016: 3' end of sense  = right_end (genomic)
  #                  3' end of antisense = left_start (genomic, lowest coord)
  dist_s  <- abs(guide$cut_site - right_end)
  dist_as <- abs(guide$cut_site - left_start)
  rec     <- if (dist_s <= dist_as) "sense" else "antisense"

  # SNP position in sense template (1-based)
  snp_idx_sense <- nchar(left_arm) + 1L
  # SNP position in antisense template (revcomp)
  snp_idx_anti  <- L - snp_idx_sense - alt_len + 2L

  list(sense           = sense,
       antisense       = antisense,
       recommended     = rec,
       snp_idx_sense   = snp_idx_sense,
       snp_idx_anti    = snp_idx_anti,
       alt_rc          = as.character(reverseComplement(DNAString(alt))),
       left_arm_len    = nchar(left_arm),
       right_arm_len   = nchar(right_arm),
       total_len       = L,
       dist_3p_sense   = dist_s,
       dist_3p_anti    = dist_as,
       pam_note        = pam_note)
}

# Render a monospace sequence with the SNP region highlighted in colour
seq_display <- function(seq, snp_idx, snp_seq) {
  pre  <- substr(seq, 1L, snp_idx - 1L)
  mid  <- substr(seq, snp_idx, snp_idx + nchar(snp_seq) - 1L)
  post <- substr(seq, snp_idx + nchar(snp_seq), nchar(seq))
  div(class = "seq-box",
    tags$span(style = "color:#888", pre),
    tags$span(style = "color:#c0392b; font-weight:bold; text-decoration:underline;", mid),
    tags$span(style = "color:#888", post)
  )
}

# ── UI ────────────────────────────────────────────────────────────────────────

CHRS <- paste0("chr", c(1:22, "X", "Y", "M"))

ui <- page_fluid(
  theme = bs_theme(bootswatch = "flatly"),

  tags$head(tags$style(HTML("
    .snp-card  { background:#fff8e1; border-left:4px solid #f59e0b;
                 padding:10px 16px; border-radius:4px; margin-bottom:14px; }
    .seq-box   { font-family:monospace; font-size:0.82em; word-break:break-all;
                 background:#f8f8f8; border:1px solid #ddd;
                 padding:6px 10px; border-radius:3px; line-height:1.6; }
    .rec-badge { background:#276749; color:white; padding:2px 8px;
                 border-radius:8px; font-size:0.82em; }
    .alt-badge { background:#666; color:white; padding:2px 8px;
                 border-radius:8px; font-size:0.82em; }
  "))),

  titlePanel("CRISPR SNP Knock-in Designer (hg38 / SpCas9 NGG)"),

  layout_sidebar(
    sidebar = sidebar(
      width = 295,

      h5("SNP Coordinates (hg38, 1-based)"),
      selectInput("chr", "Chromosome", CHRS, "chr17"),
      numericInput("pos",  "Position",          7674220, 1),
      textInput("ref",     "Reference Allele",  "G",  placeholder = "e.g. A"),
      textInput("alt",     "Alternative Allele","A",  placeholder = "e.g. T"),
      helpText("Ref allele is verified against the hg38 genome on submission."),

      hr(),
      h5("Guide RNA Search"),
      sliderInput("guide_win", "Window around SNP (bp)", 10, 100, 50, step = 5),
      numericInput("glen",     "Guide Length (nt)", 20L, 17L, 24L),

      hr(),
      h5("ssODN Repair Template"),
      numericInput("arm_len",  "Homology Arm Length (bp)", 60L, 20L, 150L),
      checkboxInput("mut_pam", "Add PAM-disrupting mutation to template", TRUE),
      checkboxInput("both_strands", "Show both ssODN orientations", TRUE),
      helpText("Recommended strand follows Richardson et al. 2016",
               "(3’ end proximal to cut site improves HDR rate)."),

      hr(),
      actionButton("go", "Design", class = "btn-primary w-100"),
      br(), br(),
      downloadButton("dl_guides",    "Download Guides CSV",    class = "btn-outline-secondary w-100"),
      br(), br(),
      downloadButton("dl_templates", "Download Templates CSV", class = "btn-outline-secondary w-100")
    ),

    uiOutput("snp_card"),

    tabsetPanel(
      id = "tabs",
      tabPanel("Guide RNAs",
        br(),
        p(tags$b("Select one guide"), "in the table to see its repair template below."),
        DTOutput("guide_tbl"),
        br(),
        uiOutput("template_panel")
      ),
      tabPanel("Locus Schematic",
        br(),
        plotOutput("schematic", height = "460px")
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  rv_guides <- reactiveVal(NULL)
  rv_snp    <- reactiveVal(NULL)

  observeEvent(input$go, {
    rv_guides(NULL); rv_snp(NULL)

    tryCatch({
      chr <- input$chr
      pos <- as.integer(input$pos)
      ref <- toupper(trimws(input$ref))
      alt <- toupper(trimws(input$alt))

      validate(
        need(nchar(ref) > 0 && !grepl("[^ACGT]", ref), "Ref must be A/C/G/T only"),
        need(nchar(alt) > 0 && !grepl("[^ACGT]", alt), "Alt must be A/C/G/T only"),
        need(chr %in% seqnames(GENOME), paste(chr, "not in hg38"))
      )

      # Verify ref against genome
      genome_ref <- toupper(as.character(
        getSeq(GENOME, GRanges(chr, IRanges(pos, pos + nchar(ref) - 1L)))))
      validate(need(genome_ref == ref,
        paste0("Reference mismatch: genome has '", genome_ref, "' at ",
               chr, ":", pos, ", you provided '", ref, "'")))

      rv_snp(list(chr = chr, pos = pos, ref = ref, alt = alt))

      # Fetch sequence window large enough for guides
      margin    <- as.integer(input$guide_win) + as.integer(input$glen) + 5L
      win_start <- max(1L, pos - margin)
      win_end   <- min(seqlengths(GENOME)[chr], pos + nchar(ref) - 1L + margin)
      seq_str   <- toupper(as.character(
        getSeq(GENOME, GRanges(chr, IRanges(win_start, win_end)))))

      guides <- find_guides(seq_str, chr, win_start, pos, nchar(ref), as.integer(input$glen))

      if (!is.null(guides)) {
        guides <- guides %>%
          filter(abs(dist_to_snp) <= input$guide_win) %>%
          mutate(gc_content = vapply(guide_sequence, gc_pct, numeric(1))) %>%
          arrange(abs(dist_to_snp))
      }

      validate(need(!is.null(guides) && nrow(guides) > 0,
        "No guides found in this window. Try increasing the search window."))

      rv_guides(guides)

    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 12)
    })
  })

  # ── SNP card ──
  output$snp_card <- renderUI({
    req(rv_snp())
    sn <- rv_snp()
    n  <- if (!is.null(rv_guides())) nrow(rv_guides()) else 0L
    div(class = "snp-card",
      tags$b(sn$chr), ":", format(sn$pos, big.mark = ","), "  |  ",
      tags$b("REF: "), sn$ref, " → ", tags$b("ALT: "), sn$alt, "  |  ",
      tags$b("Guides found: "), n
    )
  })

  # ── Guide table ──
  output$guide_tbl <- renderDT({
    req(rv_guides())
    df <- rv_guides() %>%
      transmute(
        `Guide Sequence (5'→3')`  = guide_sequence,
        `PAM`                          = pam,
        `Strand`                       = guide_strand,
        `Chr`                          = chromosome,
        `Guide Start`                  = format(guide_start, big.mark=",", scientific=FALSE),
        `Cut Site`                     = format(cut_site,    big.mark=",", scientific=FALSE),
        `Dist to SNP (bp)`             = dist_to_snp,
        `GC (%)`                       = gc_content,
        `SNP in Guide`                 = ifelse(snp_in_guide, "Yes", ""),
        `SNP in PAM`                   = ifelse(snp_in_pam,   "Yes", ""),
        `Disrupts Re-cut`              = ifelse(disrupts_recut, "✓", "")
      )

    datatable(df, rownames = FALSE, selection = "single", filter = "top",
      options = list(pageLength = 15, scrollX = TRUE, dom = "lftip",
        columnDefs = list(list(className = "dt-center", targets = 1:10)))) %>%
      formatStyle("Dist to SNP (bp)",
        background = styleInterval(c(-20.5, 0.5, 20.5),
          c("#fff3c8","#d4edda","#d4edda","#fff3c8")), fontWeight = "bold") %>%
      formatStyle("Disrupts Re-cut",
        color = styleEqual("✓", "#276749"), fontWeight = "bold") %>%
      formatStyle("Strand",
        color = styleEqual(c("+","-"), c("#2166ac","#d73027")))
  }, server = TRUE)

  # ── Repair template panel ──
  output$template_panel <- renderUI({
    req(rv_guides(), rv_snp())
    sel <- input$guide_tbl_rows_selected
    if (!length(sel))
      return(div(class = "alert alert-info mt-2",
        "Select a guide above to generate its ssODN repair template."))

    guide <- as.list(rv_guides()[sel[1L], ])
    sn    <- rv_snp()

    od <- design_ssODN(sn$chr, sn$pos, sn$ref, sn$alt,
                        input$arm_len, guide, input$mut_pam)

    rec <- od$recommended

    sense_box <- div(style = "margin-bottom:10px;",
      tags$span(class = if (rec=="sense") "rec-badge" else "alt-badge",
        if (rec=="sense") "Sense — RECOMMENDED (Richardson 2016)" else "Sense"),
      tags$br(),
      seq_display(od$sense, od$snp_idx_sense, sn$alt)
    )

    anti_box <- div(style = "margin-bottom:10px;",
      tags$span(class = if (rec=="antisense") "rec-badge" else "alt-badge",
        if (rec=="antisense") "Antisense — RECOMMENDED (Richardson 2016)" else "Antisense"),
      tags$br(),
      seq_display(od$antisense, od$snp_idx_anti, od$alt_rc)
    )

    tagList(
      hr(),
      h5(paste0("Repair Template  —  guide on ", guide$guide_strand,
                " strand, cut at ", format(guide$cut_site, big.mark=","))),

      div(style = "font-size:0.9em; margin-bottom:8px;",
        tags$b("Left arm: "), od$left_arm_len, " bp  |  ",
        tags$b("ALT allele: "), sn$alt, "  |  ",
        tags$b("Right arm: "), od$right_arm_len, " bp  |  ",
        tags$b("Total: "), od$total_len, " nt"
      ),

      div(style = "font-size:0.9em; margin-bottom:8px; color:#555;",
        "3’ end distance to cut site: ",
        tags$b(paste0("sense = ", od$dist_3p_sense, " bp")), ",  ",
        tags$b(paste0("antisense = ", od$dist_3p_anti, " bp"))
      ),

      if (nchar(od$pam_note) > 0)
        div(class = "alert alert-warning",
            style = "padding:6px 10px; font-size:0.85em; margin-bottom:8px;",
            tags$b("PAM mutation applied: "), od$pam_note),

      if (!guide$disrupts_recut && !input$mut_pam)
        div(class = "alert alert-warning",
            style = "padding:6px 10px; font-size:0.85em; margin-bottom:8px;",
            tags$b("Note: "), "The SNP does not fall within the guide or PAM. ",
            "Consider enabling 'Add PAM-disrupting mutation' to prevent Cas9 ",
            "from re-cutting the edited allele."),

      sense_box,
      if (input$both_strands) anti_box,

      div(style = "font-size:0.82em; color:#666; margin-top:4px;",
        "Red underline = ALT allele (or its reverse complement in the antisense strand). ",
        "Grey = homology arms.")
    )
  })

  # ── Locus schematic ──
  output$schematic <- renderPlot({
    req(rv_guides(), rv_snp())
    sn     <- rv_snp()
    guides <- rv_guides()

    pad   <- input$guide_win + as.integer(input$glen) + 10L
    xmin  <- sn$pos - pad
    xmax  <- sn$pos + nchar(sn$ref) - 1L + pad

    # Assign y levels for guides (avoid crowding)
    n_g   <- nrow(guides)
    y_lvl <- seq(0.7, 0.7 + (n_g - 1) * 0.55, length.out = n_g)

    p <- ggplot() +
      # Genome backbone
      annotate("segment", x = xmin, xend = xmax, y = 0, yend = 0,
               color = "black", linewidth = 1.8) +
      # SNP marker
      geom_vline(xintercept = sn$pos, color = "#c0392b",
                 linetype = "dashed", linewidth = 1.1) +
      annotate("label",
        x = sn$pos, y = max(y_lvl) + 0.85,
        label = paste0(sn$ref, "→", sn$alt, "\n",
                       sn$chr, ":", format(sn$pos, big.mark=",")),
        color = "#c0392b", fill = "white", size = 3.5, fontface = "bold") +
      # ssODN bar
      annotate("segment",
        x    = sn$pos - input$arm_len,
        xend = sn$pos + nchar(sn$ref) - 1L + input$arm_len,
        y = -0.55, yend = -0.55, color = "#276749", linewidth = 5, alpha = 0.45) +
      annotate("text",
        x = sn$pos, y = -0.82,
        label = paste0("ssODN  (", input$arm_len, " bp + ALT + ", input$arm_len, " bp)"),
        size = 3.4, color = "#276749")

    gcolors <- c("+" = "#2166ac", "-" = "#d73027")

    for (i in seq_len(n_g)) {
      g   <- guides[i, ]
      col <- gcolors[g$guide_strand]
      y   <- y_lvl[i]

      # Guide body
      p <- p +
        annotate("segment",
          x = g$guide_start, xend = g$guide_end,
          y = y, yend = y, color = col, linewidth = 4, alpha = 0.65) +
        # PAM tick
        annotate("segment",
          x    = g$pam_start, xend = g$pam_start + 2L,
          y = y, yend = y, color = col, linewidth = 4, alpha = 1) +
        # Cut site marker
        annotate("segment",
          x = g$cut_site, xend = g$cut_site,
          y = y - 0.2, yend = y + 0.2,
          color = "black", linewidth = 1.2) +
        # Label
        annotate("text",
          x = (g$guide_start + g$guide_end) / 2,
          y = y + 0.28,
          label = paste0(
            g$guide_strand, " | d=", g$dist_to_snp, "bp",
            if (g$disrupts_recut) "  ✓recut" else ""),
          size = 3, color = col)
    }

    p +
      scale_x_continuous(
        name   = paste0(sn$chr, " coordinate (hg38)"),
        labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
      labs(
        title    = paste0("Knock-in locus: ", sn$chr, ":",
                           format(sn$pos, big.mark=","), "  ", sn$ref, "→", sn$alt),
        subtitle = paste0(n_g, " guide(s)  |  blue = +strand, red = −strand  |",
                          "  | = cut site  |  ✓recut = SNP disrupts re-cutting"),
        y        = NULL) +
      ylim(-1.1, max(y_lvl) + 1.3) +
      theme_classic(base_size = 13) +
      theme(axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y  = element_blank(),
            plot.title   = element_text(face = "bold"),
            plot.subtitle = element_text(size = 10, color = "grey40"))
  })

  # ── Downloads ──
  output$dl_guides <- downloadHandler(
    filename = function()
      paste0("guides_", input$chr, "_", input$pos, "_", Sys.Date(), ".csv"),
    content = function(f) {
      req(rv_guides())
      write.csv(rv_guides(), f, row.names = FALSE)
    }
  )

  output$dl_templates <- downloadHandler(
    filename = function()
      paste0("templates_", input$chr, "_", input$pos, "_", Sys.Date(), ".csv"),
    content = function(f) {
      req(rv_guides(), rv_snp())
      sn <- rv_snp()
      rows <- lapply(seq_len(nrow(rv_guides())), function(i) {
        g  <- as.list(rv_guides()[i, ])
        od <- design_ssODN(sn$chr, sn$pos, sn$ref, sn$alt,
                            input$arm_len, g, input$mut_pam)
        data.frame(
          guide_sequence  = g$guide_sequence,
          guide_strand    = g$guide_strand,
          cut_site        = g$cut_site,
          dist_to_snp     = g$dist_to_snp,
          disrupts_recut  = g$disrupts_recut,
          ssODN_sense     = od$sense,
          ssODN_antisense = od$antisense,
          recommended     = od$recommended,
          total_len_nt    = od$total_len,
          pam_note        = od$pam_note,
          stringsAsFactors = FALSE
        )
      })
      write.csv(bind_rows(rows), f, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
