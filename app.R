library(shiny)
library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)

# ---------------------------------------------------------------------------
# Load data (once, at startup)
# ---------------------------------------------------------------------------
term_data   <- readRDS("data/term.rds")
lang_meta   <- readRDS("data/lang_meta.rds")
dict_data   <- readRDS("data/dict.rds")
chip_layout <- readRDS("data/chip_layout.rds")
cielab      <- readRDS("data/cielab.rds")
spkr_meta   <- readRDS("data/spkr_meta.rds")
foci_data   <- readRDS("data/foci.rds")

wcs_lang_csv <- read.csv("wcs_languages.csv", stringsAsFactors = FALSE)
colnames(wcs_lang_csv)[1:4] <- c("lang_id", "language_csv", "family", "country_csv")

lang_choices <- setNames(lang_meta$lang_id, lang_meta$language)

# ---------------------------------------------------------------------------
# Helper: convert L*, a*, b* vectors to sRGB hex strings
# ---------------------------------------------------------------------------
lab_to_hex <- function(L, a, b) {
  mat     <- cbind(L, a, b)
  rgb_raw <- grDevices::convertColor(mat, from = "Lab", to = "sRGB")
  rgb_mat <- matrix(pmax(0, pmin(1, as.numeric(rgb_raw))), ncol = 3)
  grDevices::rgb(rgb_mat[, 1], rgb_mat[, 2], rgb_mat[, 3])
}

# ---------------------------------------------------------------------------
# Helper: black or white text depending on background luminance
# ---------------------------------------------------------------------------
contrast_text <- function(hex_colors) {
  rgb_mat <- col2rgb(hex_colors) / 255
  lum     <- 0.299 * rgb_mat[1, ] + 0.587 * rgb_mat[2, ] + 0.114 * rgb_mat[3, ]
  ifelse(lum > 0.45, "black", "white")
}

# ---------------------------------------------------------------------------
# Helper: find the chip closest to the CIELAB centroid for each term
# ---------------------------------------------------------------------------
centroid_chips <- function(naming_rows, cielab) {
  naming_rows |>
    left_join(cielab, by = "chip_id") |>
    group_by(term) |>
    summarise(
      cL = mean(L_star, na.rm = TRUE),
      ca = mean(a_star, na.rm = TRUE),
      cb = mean(b_star, na.rm = TRUE),
      .groups = "drop"
    ) |>
    left_join(
      naming_rows |> left_join(cielab, by = "chip_id"),
      by = "term"
    ) |>
    mutate(dist2 = (L_star - cL)^2 + (a_star - ca)^2 + (b_star - cb)^2) |>
    group_by(term) |>
    slice_min(dist2, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(term, chip_id)
}

# ---------------------------------------------------------------------------
# Helper: add focal-choice markers to an existing ggplot
# focal_chip_ids: integer vector of chip_ids to mark
# ---------------------------------------------------------------------------
add_focal_overlay <- function(p, focal_chip_ids, size = 1.5, stroke = 0.7) {
  if (length(focal_chip_ids) == 0) return(p)
  focal_pos <- chip_layout |> filter(chip_id %in% focal_chip_ids)
  p + geom_point(
    data        = focal_pos,
    aes(x = col_num, y = row_num),
    inherit.aes = FALSE,
    shape = 21, size = size, stroke = stroke, fill = "white", color = "black"
  )
}

# ---------------------------------------------------------------------------
# Helper: build a mode-map ggplot
#   naming_rows  : data frame with chip_id, term
#   labels_df    : data frame with chip_id, term (centroid chips for annotation)
#   palette      : named character vector term -> hex color
#   title        : plot title string or NULL
# ---------------------------------------------------------------------------
make_mode_map <- function(naming_rows, labels_df, palette, title = NULL,
                          axis_size = 7, title_size = 8, label_size = 2.2) {
  text_col_map        <- contrast_text(palette)
  names(text_col_map) <- names(palette)

  plot_df <- chip_layout |>
    left_join(naming_rows |> select(chip_id, term), by = "chip_id") |>
    left_join(labels_df  |> select(chip_id, label = term), by = "chip_id") |>
    mutate(text_col = text_col_map[term])

  ggplot(plot_df, aes(x = col_num, y = row_num, fill = term)) +
    geom_tile(color = "grey80", linewidth = 0.08) +
    geom_text(
      data        = ~ filter(.x, !is.na(label)),
      aes(label = label, color = text_col),
      size = label_size, fontface = "bold",
      inherit.aes = TRUE
    ) +
    scale_fill_manual(values = palette, na.value = "#DDDDDD", na.translate = FALSE) +
    scale_color_identity() +
    scale_y_reverse(breaks = 1:10, labels = LETTERS[1:10], expand = c(0.04, 0.04)) +
    scale_x_continuous(breaks = c(0, seq(10, 40, 10)), expand = c(0.015, 0.015)) +
    coord_fixed() +
    labs(title = title) +
    theme_minimal() +
    theme(
      panel.grid      = element_blank(),
      axis.title      = element_blank(),
      axis.text       = element_text(size = axis_size),
      legend.position = "none",
      plot.title      = element_text(size = title_size, face = "bold")
    )
}

# ---------------------------------------------------------------------------
# Helper: build an agreement-threshold map
#   chip_agree : data frame with chip_id, term, frac (modal term + its fraction)
#   threshold  : minimum fraction to color a chip
# ---------------------------------------------------------------------------
make_agree_map <- function(chip_agree, palette, threshold, title = NULL,
                           axis_size = 7, title_size = 8) {
  plot_df <- chip_layout |>
    left_join(chip_agree, by = "chip_id") |>
    mutate(fill_term = if_else(!is.na(frac) & frac >= threshold, term, NA_character_))

  ggplot(plot_df, aes(x = col_num, y = row_num, fill = fill_term)) +
    geom_tile(color = "grey80", linewidth = 0.08) +
    scale_fill_manual(values = palette, na.value = "#DDDDDD", na.translate = FALSE) +
    scale_y_reverse(breaks = 1:10, labels = LETTERS[1:10], expand = c(0.04, 0.04)) +
    scale_x_continuous(breaks = c(0, seq(10, 40, 10)), expand = c(0.015, 0.015)) +
    coord_fixed() +
    labs(title = title) +
    theme_minimal() +
    theme(
      panel.grid      = element_blank(),
      axis.title      = element_blank(),
      axis.text       = element_text(size = axis_size),
      legend.position = "none",
      plot.title      = element_text(size = title_size, face = "bold")
    )
}

# ---------------------------------------------------------------------------
# Helper: build the average-individual map for a selected term
#   chip_prob : data frame with chip_id, prob  (fraction of speakers naming chip with term)
#   term_color: hex color for the term
# ---------------------------------------------------------------------------
make_avg_ind_map <- function(chip_prob, term_color, title = NULL,
                             axis_size = 7, title_size = 8) {
  plot_df <- chip_layout |>
    left_join(chip_prob, by = "chip_id") |>
    mutate(prob = replace_na(prob, 0))

  ggplot(plot_df, aes(x = col_num, y = row_num, fill = prob)) +
    geom_tile(color = "grey80", linewidth = 0.08) +
    scale_fill_gradient(low = "#DDDDDD", high = term_color, limits = c(0, 1)) +
    scale_y_reverse(breaks = 1:10, labels = LETTERS[1:10], expand = c(0.04, 0.04)) +
    scale_x_continuous(breaks = c(0, seq(10, 40, 10)), expand = c(0.015, 0.015)) +
    coord_fixed() +
    labs(title = title) +
    theme_minimal() +
    theme(
      panel.grid      = element_blank(),
      axis.title      = element_blank(),
      axis.text       = element_text(size = axis_size),
      legend.position = "none",
      plot.title      = element_text(size = title_size, face = "bold")
    )
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .control-label       { font-size: 16px; font-weight: bold; }
    .selectize-input     { font-size: 15px; }
    .selectize-dropdown  { font-size: 15px; }
  ")),
  # When embedded in an iframe, report our content height to the parent page so
  # it can size the iframe to fit (no fixed-height whitespace or inner scrollbar).
  # Fire on Shiny's render events (a ResizeObserver on body misses the plot
  # images, which overflow the body box without enlarging it).
  tags$script(HTML("
    function wcsSendHeight() {
      if (window.parent === window) return;
      // body.scrollHeight is the content height; documentElement.scrollHeight
      // would include the iframe viewport and cause a resize feedback loop.
      var h = document.body.scrollHeight;
      window.parent.postMessage({ wcsHeight: h }, '*');
    }
    window.addEventListener('load', wcsSendHeight);
    window.addEventListener('resize', wcsSendHeight);
    $(document).on('shiny:value shiny:recalculated shiny:idle', function() {
      setTimeout(wcsSendHeight, 150);
    });
  "))),
  titlePanel("World Color Survey"),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      selectizeInput("lang_id", "Language:",
        choices  = lang_choices,
        selected = 1
      ),
      selectizeInput("term", "Term:",
        choices  = "All",
        selected = "All"
      )
    ),
    mainPanel(
      width = 10,
      uiOutput("lang_info_bar"),
      plotOutput("top_maps", height = "auto"),
      hr(),
      tags$h4("Individual speakers", style = "margin-bottom: 2px;"),
      plotOutput("indiv_maps", height = "auto")
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # Naming data for the selected language
  lang_terms <- reactive({
    term_data |> filter(lang_id == input$lang_id)
  })

  n_spkrs <- reactive({
    n_distinct(lang_terms()$spkr_id)
  })

  # Term abbreviations for the selected language, for term menu
  observeEvent(input$lang_id, {
    term_df <- dict_data |>
      filter(lang_id == input$lang_id) |>
      arrange(abbreviation) |>
      mutate(label = paste0(translation, " (", abbreviation, ")"))
    choices <- c("All" = "All", setNames(term_df$abbreviation, term_df$label))
    updateSelectizeInput(session, "term", choices = choices, selected = "All")
  })

  # Color palette: CIELAB centroid of each term across all speakers
  palette <- reactive({
    lang_terms() |>
      left_join(cielab, by = "chip_id") |>
      group_by(term) |>
      summarise(
        L = mean(L_star, na.rm = TRUE),
        a = mean(a_star, na.rm = TRUE),
        b = mean(b_star, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(hex = lab_to_hex(L, a, b)) |>
      select(term, hex) |>
      tibble::deframe()
  })

  # Population mode term per chip
  pop_naming <- reactive({
    lang_terms() |>
      count(chip_id, term) |>
      slice_max(n, n = 1, with_ties = FALSE, by = chip_id) |>
      select(chip_id, term)
  })

  # Modal term + agreement fraction per chip (for agreement maps)
  chip_agree <- reactive({
    lang_terms() |>
      count(chip_id, term) |>
      slice_max(n, n = 1, with_ties = FALSE, by = chip_id) |>
      mutate(frac = n / n_spkrs())
  })

  # Centroid-closest chip per term, for population map labels
  pop_labels <- reactive({
    centroid_chips(pop_naming(), cielab)
  })

  # Focal data for the selected language
  lang_foci <- reactive({
    foci_data |> filter(lang_id == input$lang_id)
  })

  # ---------------------------------------------------------------------------
  # Language metadata bar
  # ---------------------------------------------------------------------------
  output$lang_info_bar <- renderUI({
    meta <- lang_meta |> filter(lang_id == input$lang_id)
    if (nrow(meta) == 0) return(NULL)

    csv_row <- wcs_lang_csv |> filter(lang_id == as.integer(input$lang_id))
    family  <- if (nrow(csv_row) > 0) csv_row$family[1]      else "—"
    country <- if (nrow(csv_row) > 0) csv_row$country_csv[1] else "—"

    fw <- na.omit(c(meta$c1, meta$c2, meta$c3))
    fw <- fw[nchar(fw) > 0]
    fw_str <- if (length(fw) > 0) paste(fw, collapse = " · ") else "—"

    wellPanel(
      style = "padding: 6px 14px; background: #f0f4f8; border: 1px solid #d0d8e4; margin-bottom: 10px;",
      fluidRow(
        column(3, tags$b("Language: "), meta$language),
        column(2, tags$b("Family: "),   family),
        column(2, tags$b("Country: "),  country),
        column(5, tags$b("Fieldworkers: "), fw_str)
      )
    )
  })

  # ---------------------------------------------------------------------------
  # Top maps: 4 agreement maps (All) or term mode + average individual (term selected)
  # ---------------------------------------------------------------------------
  output$top_maps <- renderPlot({
    sel <- input$term
    pal <- palette()
    ca  <- chip_agree()

    if (sel == "All") {
      # Mode map
      p_mode <- make_mode_map(
        naming_rows = pop_naming(),
        labels_df   = pop_labels(),
        palette     = pal,
        title       = "Mode map",
        axis_size   = 6, title_size = 15, label_size = 2.8
      )
      # Agreement maps
      p30  <- make_agree_map(ca, pal, 0.30, title = "30% agreement",
                             axis_size = 6, title_size = 15)
      p70  <- make_agree_map(ca, pal, 0.70, title = "70% agreement",
                             axis_size = 6, title_size = 15)
      p100 <- make_agree_map(ca, pal, 1.00, title = "100% agreement",
                             axis_size = 6, title_size = 15)
      (p_mode | p30) / (p70 | p100)

    } else {
      # Term mode map: mode map with other terms blanked
      term_naming <- pop_naming() |>
        mutate(term = if_else(term == sel, term, NA_character_))
      term_labels <- pop_labels() |> filter(term == sel)

      p_term <- make_mode_map(
        naming_rows = term_naming,
        labels_df   = term_labels[0, ],
        palette     = pal,
        title       = paste0("Mode map — ", sel),
        axis_size   = 6, title_size = 15, label_size = 2.8
      )

      # Average individual map: fraction of speakers who named each chip with sel
      chip_prob <- lang_terms() |>
        filter(term == sel) |>
        count(chip_id) |>
        mutate(prob = n / n_spkrs())

      term_color <- if (sel %in% names(pal)) pal[[sel]] else "#888888"
      p_avg <- make_avg_ind_map(
        chip_prob  = chip_prob,
        term_color = term_color,
        title      = paste0("Avg. individual — ", sel),
        axis_size  = 6, title_size = 15
      )

      # Focal overlay: all speakers' focal choices for this term
      focal_ids <- lang_foci() |>
        filter(term == sel) |>
        pull(chip_id) |>
        unique()

      p_term <- add_focal_overlay(p_term, focal_ids)
      p_avg  <- add_focal_overlay(p_avg,  focal_ids)

      p_term | p_avg
    }
  }, height = function() {
    # Size the plot to the maps it holds so they don't float in whitespace.
    # Each map is coord_fixed with a 10-row x 41-col grid drawn in a 2-column
    # layout: 2 rows of maps for "All", 1 row when a single term is selected.
    w <- session$clientData$output_top_maps_width
    if (is.null(w) || w < 10) {
      return(if (!is.null(input$term) && input$term != "All") 210L else 380L)
    }
    n_rows   <- if (!is.null(input$term) && input$term != "All") 1L else 2L
    panel_h  <- ((w - 60L) / 2) * (10 / 41)  # per-map panel height, minus axis/spacing
    as.integer(n_rows * (panel_h + 40L) + 15L)  # +40 per row for title & x-axis labels
  })

  # ---------------------------------------------------------------------------
  # Individual speaker maps
  # ---------------------------------------------------------------------------

  # All speaker–chip combinations with per-speaker CIELAB fill colors
  indiv_all <- reactive({
    lt    <- lang_terms()
    spkrs <- sort(unique(lt$spkr_id))

    spkr_pal <- lt |>
      left_join(cielab, by = "chip_id") |>
      group_by(spkr_id, term) |>
      summarise(
        L = mean(L_star, na.rm = TRUE),
        a = mean(a_star, na.rm = TRUE),
        b = mean(b_star, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(fill_hex = lab_to_hex(L, a, b))

    expand_grid(spkr_id = spkrs, chip_id = chip_layout$chip_id) |>
      left_join(chip_layout, by = "chip_id") |>
      left_join(lt |> select(spkr_id, chip_id, term), by = c("spkr_id", "chip_id")) |>
      left_join(spkr_pal |> select(spkr_id, term, fill_hex), by = c("spkr_id", "term")) |>
      mutate(fill_hex = replace_na(fill_hex, "#DDDDDD"))
  })

  # Centroid-closest chip per (speaker, term) for individual annotations
  indiv_labels <- reactive({
    lt <- lang_terms()
    lt |>
      left_join(cielab, by = "chip_id") |>
      group_by(spkr_id, term) |>
      summarise(
        cL = mean(L_star, na.rm = TRUE),
        ca = mean(a_star, na.rm = TRUE),
        cb = mean(b_star, na.rm = TRUE),
        .groups = "drop"
      ) |>
      left_join(lt |> left_join(cielab, by = "chip_id"),
                by = c("spkr_id", "term")) |>
      mutate(dist2 = (L_star - cL)^2 + (a_star - ca)^2 + (b_star - cb)^2) |>
      group_by(spkr_id, term) |>
      slice_min(dist2, n = 1, with_ties = FALSE) |>
      ungroup() |>
      select(spkr_id, term, chip_id)
  })

  # Speaker facet labels: "Speaker N · age sex"
  spkr_label_map <- reactive({
    meta <- spkr_meta |> filter(lang_id == input$lang_id)
    setNames(
      paste0("Speaker ", meta$spkr_id, " · ", meta$age, " ", meta$sex),
      as.character(meta$spkr_id)
    )
  })

  output$indiv_maps <- renderPlot({
    sel <- input$term
    ia  <- indiv_all()
    il  <- indiv_labels()

    # Apply term filter
    if (sel != "All") {
      ia <- ia |> mutate(
        fill_hex = if_else(!is.na(term) & term != sel, "#DDDDDD", fill_hex),
        term     = if_else(term == sel, term, NA_character_)
      )
      il <- il |> filter(term == sel)
    }

    plot_df <- ia |>
      left_join(il |> select(spkr_id, chip_id, label = term),
                by = c("spkr_id", "chip_id")) |>
      mutate(
        label    = if (sel != "All") NA_character_ else label,
        text_col = contrast_text(fill_hex)
      )

    # Focal data for individual panels
    foc <- lang_foci()
    if (sel != "All") foc <- foc |> filter(term == sel)
    foc_pos <- foc |>
      left_join(chip_layout |> select(chip_id, col_num, row_num), by = "chip_id") |>
      mutate(spkr_id = as.integer(spkr_id))

    n_spkr_val <- n_spkrs()
    n_cols     <- min(n_spkr_val, 5L)
    lbl_map    <- spkr_label_map()

    p <- ggplot(plot_df, aes(x = col_num, y = row_num, fill = fill_hex)) +
      geom_tile(color = "grey80", linewidth = 0.04) +
      # geom_text(
      #   data        = ~ filter(.x, !is.na(label)),
      #   aes(label = label, color = text_col),
      #   size = 3.2, fontface = "bold",
      #   inherit.aes = TRUE
      # ) +
      scale_fill_identity() +
      scale_color_identity() +
      scale_y_reverse() +
      scale_x_continuous() +
      facet_wrap(
        ~spkr_id, ncol = n_cols,
        labeller = labeller(spkr_id = lbl_map)
      ) +
      coord_fixed() +
      theme_minimal() +
      theme(
        panel.grid      = element_blank(),
        axis.title      = element_blank(),
        axis.text       = element_blank(),
        strip.text      = element_text(size = 12),
        legend.position = "none"
      )

    if (nrow(foc_pos) > 0) {
      p <- p + geom_point(
        data        = foc_pos,
        aes(x = col_num, y = row_num),
        inherit.aes = FALSE,
        shape = 21, size = 0.9, stroke = 0.5, fill = "white", color = "black"
      )
    }

    p
  },
  height = function() {
    w <- session$clientData$output_indiv_maps_width
    if (is.null(w) || w < 10) return(300L)
    n_spkr_val <- n_spkrs()
    n_cols     <- min(n_spkr_val, 5L)
    n_rows     <- ceiling(n_spkr_val / n_cols)
    panel_h    <- ((w - 30L) / n_cols) * (10 / 41)
    as.integer(n_rows * (panel_h + 24L) + 15L)
  })
}

shinyApp(ui, server)
