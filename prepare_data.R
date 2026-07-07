# prepare_data.R
# Run once from wcs/ before launching.
# Reads raw WCS data files and writes compact RDS files into data/.

library(tidyverse)

raw <- "./rawdata/wcs"

message("Reading term data...")
term_data <- read_tsv(
  file.path(raw, "term.txt"),
  col_names = c("lang_id", "spkr_id", "chip_id", "term"),
  col_types  = "iiic",
  na         = character()   # "NA" is a valid WCS term abbreviation
)
saveRDS(term_data, "data/term.rds")
message(sprintf("  %d rows saved to data/term.rds", nrow(term_data)))

message("Reading language metadata...")
macrom_lut <- c(
  "87" = "á",  # á
  "8A" = "ä",  # ä
  "8B" = "ã",  # ã
  "96" = "ñ",  # ñ
  "97" = "ó",  # ó
  "9C" = "ú"   # ú
)
decode_macrom <- function(x) {
  for (code in names(macrom_lut))
    x <- gsub(paste0("\\{\\\\x", code, "\\}"), macrom_lut[[code]], x, ignore.case = TRUE)
  x
}
lang_meta <- read_tsv(
  file.path(raw, "lang.txt"),
  col_names = c("lang_id", "language", "country", "c1", "c2", "c3", "filename", "status"),
  col_types  = "iccccccc"
) |>
  select(lang_id, language, country, c1, c2, c3) |>
  mutate(
    across(where(is.character), decode_macrom),
    across(c(c1, c2, c3), ~ na_if(.x, "*"))
  )
saveRDS(lang_meta, "data/lang_meta.rds")
message(sprintf("  %d languages saved to data/lang_meta.rds", nrow(lang_meta)))

message("Reading dictionary...")
dict_data <- read_tsv(
  file.path(raw, "dict.txt"),
  col_names = c("lang_id", "term_num", "translation", "abbreviation"),
  col_types  = "iicc",
  comment    = "#",
  na         = character()
)
saveRDS(dict_data, "data/dict.rds")
message(sprintf("  %d terms saved to data/dict.rds", nrow(dict_data)))

message("Reading chip layout...")
chip_layout <- read_tsv(
  file.path(raw, "chip.txt"),
  col_names = c("chip_id", "row_letter", "col_num", "chip_name"),
  col_types  = "icic"
) |>
  mutate(row_num = match(row_letter, LETTERS))
saveRDS(chip_layout, "data/chip_layout.rds")
message(sprintf("  %d chips saved to data/chip_layout.rds", nrow(chip_layout)))

message("Reading CIELAB chip colors...")
cielab <- read_tsv(
  "rawdata/cnum-vhcm-lab-new.txt",
  col_names = c("chip_id", "V", "H", "C_val", "MunH", "MunV", "L_star", "a_star", "b_star"),
  col_types  = "iiiiccddd",
  comment    = "#"
) |> select(chip_id, L_star, a_star, b_star)
saveRDS(cielab, "data/cielab.rds")
message(sprintf("  %d chips saved to data/cielab.rds", nrow(cielab)))

message("Reading speaker metadata...")
spkr_meta <- read_tsv(
  file.path(raw, "spkr-lsas.txt"),
  col_names = c("lang_id", "spkr_id", "age", "sex"),
  col_types  = "iiic"
)
saveRDS(spkr_meta, "data/spkr_meta.rds")
message(sprintf("  %d speakers saved to data/spkr_meta.rds", nrow(spkr_meta)))

message("Reading focal choices...")
foci_path  <- file.path(raw, "foci-exp.txt")
foci_raw   <- readChar(foci_path, file.info(foci_path)$size, useBytes = TRUE)
foci_lines <- strsplit(foci_raw, "\r", fixed = TRUE)[[1]]
foci_df    <- read_tsv(
  I(paste(foci_lines, collapse = "\n")),
  col_names = c("lang_id", "spkr_id", "focal_index", "term", "chip_code"),
  na        = c("*", "?"),
  col_types = "iiicc"
) |>
  filter(!is.na(term), !is.na(chip_code))

# Convert chip_code (e.g. "D9") to chip_id via chip_layout$chip_name
foci <- foci_df |>
  left_join(chip_layout |> select(chip_id, chip_name), by = c("chip_code" = "chip_name")) |>
  filter(!is.na(chip_id)) |>
  select(lang_id, spkr_id, term, chip_id) |>
  distinct()
saveRDS(foci, "data/foci.rds")
message(sprintf("  %d focal choices saved to data/foci.rds", nrow(foci)))

message("Done.")
