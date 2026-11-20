# =============================================================================
# final_project.R
# Tract-level TCE exposure from EPA Superfund (SEMS) and AirToxScreen data, Pennsylvania
#
# Outputs (in data/processed/):
#   tract_superfund_tce.csv          distance to nearest TCE Superfund site, sites within 2 km
#   tract_airtox_tce.csv             population-weighted AirToxScreen TCE concentration
#   tract_tce_superfund_airtox.csv   both combined, one row per 2020 PA tract
#
# Data sources:
#   Superfund: pulled directly from the EPA Envirofacts API (no download needed)
#   AirToxScreen: downloaded from a URL you paste into CONFIG, or a file in data/raw/
#
# Run with the working directory set to the project folder
# (C:/academics/bmin5030/final). Requires a Census API key for block populations.
# =============================================================================

setwd("C:/academics/bmin5030/final")

library(tidyverse)
library(sf)
library(tigris)
library(tidycensus)
library(readxl)
library(janitor)
options(tigris_use_cache = TRUE)


# ---- CONFIG -----------------------------------------------------------------

paths <- list(
  # AirToxScreen 2020: paste the download link for the PA/TCE subset or the Region 3
  # TCE file here, or leave NULL and put the file in data/raw/ yourself
  airtox.url = "https://gaftp.epa.gov/rtrmodeling_public/AirToxScreen/2020/Pollutant%20Summaries/Region3/TRICHLOROETHYLENE_EXPCONC_R3.zip",
  airtox     = "data/raw/airtoxscreen_2020_pa_tce.csv",
  out.dir    = "data/processed",
  fig.dir    = "output/figures"
)

skip <- list(airtox = 0)   # rows above the header row in the AirToxScreen file

# AirToxScreen column names: leave NULL to auto-detect, or set the cleaned name
cols <- list(
  airtox.id        = NULL,   # block (15-digit) or tract (11-digit) FIPS
  airtox.pollutant = NULL,
  airtox.value     = NULL    # total concentration across all source groups
)

sems.states  <- c("PA", "NJ", "NY", "OH", "WV", "MD", "DE")  # PA plus neighbors
npl.statuses <- NULL    # NULL keeps all SEMS sites; or e.g. c("Currently on the Final NPL")
buffer.m     <- 2000    # "nearby" distance for site counts
border.m     <- 10000   # keep out-of-state sites within 10 km of PA
crs.analysis <- 5070    # CONUS Albers, meters

tce.pattern <- regex("trichloroeth(yl)?ene|\\btce\\b", ignore_case = TRUE)  # excludes trichloroethane

ef.base     <- "https://data.epa.gov/efservice/"
ef.max.rows <- 100000L   # integer (L) so it pastes as "100000", not "1e+05"

dir.create(paths$out.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$fig.dir, recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw",    recursive = TRUE, showWarnings = FALSE)


# ---- Helpers ----------------------------------------------------------------

# Read an Envirofacts API query as all-character columns; stop on API errors
ef_get <- function(url) {
  df <- read_csv(url, col_types = cols(.default = "c"), show_col_types = FALSE)
  if (ncol(df) == 1 && str_detect(names(df)[1], "error"))
    stop("Envirofacts API error for\n  ", url, "\n  ", names(df)[1])
  if (nrow(df) >= ef.max.rows)
    warning("Query returned ", nrow(df), " rows, the maximum requested; results may be ",
            "truncated. Increase ef.max.rows.\n  ", url)
  clean_names(df)
}

read_any <- function(path, skip = 0) {
  if (!file.exists(path)) stop("File not found: ", path, "\nCheck `paths` in CONFIG.")
  ext <- tolower(tools::file_ext(path))
  df <- switch(ext,
    csv  = read_csv(path, skip = skip, show_col_types = FALSE, guess_max = 100000),
    txt  = read_delim(path, skip = skip, show_col_types = FALSE, guess_max = 100000),
    xlsx = , xls = read_excel(path, skip = skip, guess_max = 100000),
    stop("Unsupported format .", ext, ". Save it as .csv or .xlsx and try again."))
  df <- clean_names(df)
  if (mean(str_detect(names(df), "^x\\d+$")) > 0.5)
    warning("Most column names look like x1, x2, ... The header is probably lower down: ",
            "increase `skip` for this file in CONFIG.")
  df
}

# Find a column by override or by the first matching pattern
find_col <- function(df, override, patterns, what, required = TRUE) {
  nm <- names(df)
  if (!is.null(override)) {
    if (!override %in% nm) stop("Column '", override, "' (", what, ") not found.")
    return(override)
  }
  for (p in patterns) {
    hit <- nm[str_detect(nm, p)]
    if (length(hit) > 0) {
      message(what, ": using column '", hit[1], "'")
      return(hit[1])
    }
  }
  if (!required) {
    message(what, ": no column found, skipping")
    return(NULL)
  }
  stop("Could not find the ", what, " column. Columns are:\n  ",
       paste(nm, collapse = ", "), "\nSet it in `cols` in CONFIG.")
}

norm_id <- function(x) {
  x <- if (is.numeric(x)) sprintf("%.0f", x) else as.character(x)
  toupper(str_trim(x))
}


# ---- Shared geography -------------------------------------------------------

pa.tracts.tiger <- tracts(state = "PA", year = 2020, cb = FALSE) |>
  st_transform(crs.analysis)

# Simplified (cartographic) boundaries for maps only: they exclude Lake Erie and
# other large water areas. Distances always use pa.tracts.tiger.
pa.tracts.map <- tracts(state = "PA", year = 2020, cb = TRUE) |>
  st_transform(crs.analysis)

pa.boundary <- states(year = 2020, cb = TRUE) |>
  filter(STUSPS == "PA") |>
  st_transform(crs.analysis)


# =============================================================================
# PART 1: SUPERFUND (EPA Envirofacts API, SEMS tables)
# =============================================================================

# ---- 1a. TCE contaminant records, all SEMS sites nationwide -----------------
# The API's server-side join between these tables fails, so the two tables are
# pulled separately and joined here on site_id = fk_site_id.

sems.contam.tce <- ef_get(paste0(
  ef.base, "sems.envirofacts_contaminants/preferred_contaminant_name/contains/trichloroeth/1:",
  ef.max.rows, "/CSV")) |>
  filter(str_detect(preferred_contaminant_name, tce.pattern))

# Check exactly which names matched (no trichloroethane) and which media
sems.contam.tce |> count(preferred_contaminant_name, sort = TRUE) |> print()
sems.contam.tce |> count(media_name, sort = TRUE) |> print()
message(n_distinct(sems.contam.tce$fk_site_id), " SEMS sites nationwide with TCE recorded")


# ---- 1b. Sites in PA and neighboring states ---------------------------------

sems.sites <- map(sems.states, \(st) ef_get(paste0(
  ef.base, "sems.envirofacts_site/fk_ref_state_code/equals/", st, "/1:", ef.max.rows, "/CSV"))) |>
  list_rbind()

sems.sites |> count(fk_ref_state_code) |> print()

tce.site.info <- sems.contam.tce |>
  group_by(site_id = fk_site_id) |>
  summarize(media       = paste(sort(unique(media_name)), collapse = "; "),
            groundwater = any(str_detect(media_name, regex("ground", ignore_case = TRUE)),
                              na.rm = TRUE),
            n.records   = n(),
            .groups     = "drop")

tce.sites.tbl <- sems.sites |>
  distinct(site_id, .keep_all = TRUE) |>
  inner_join(tce.site.info, by = "site_id") |>
  mutate(lat = as.numeric(primary_latitude_decimal_val),
         lon = as.numeric(primary_longitude_decimal_val))

# NPL status of TCE sites: decide whether to keep all SEMS sites or only NPL sites
tce.sites.tbl |> count(npl_status_name, sort = TRUE) |> print()
if (!is.null(npl.statuses)) {
  tce.sites.tbl <- tce.sites.tbl |> filter(npl_status_name %in% npl.statuses)
}

no.coords <- tce.sites.tbl |> filter(is.na(lat) | is.na(lon))
message(nrow(no.coords), " TCE sites have no coordinates and are dropped")
if (nrow(no.coords) > 0) print(select(no.coords, site_id, epa_id, name, fk_ref_state_code))

tce.sites <- tce.sites.tbl |>
  filter(!is.na(lat), !is.na(lon)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs.analysis)

# Keep PA sites and out-of-state sites within border.m of PA
tce.sites <- tce.sites[lengths(st_is_within_distance(tce.sites, pa.boundary, border.m)) > 0, ]
message(nrow(tce.sites), " TCE Superfund sites in or near Pennsylvania (",
        sum(tce.sites$fk_ref_state_code == "PA"), " in PA)")


# ---- 1c. Tract-level Superfund metrics --------------------------------------
# Distances are measured from each tract's 2020 population-weighted center (Census
# Bureau "centers of population"), so a large rural tract containing a site is not
# scored 0 km when most residents live far from it. Site locations are points
# (each site's recorded coordinate), not full site footprints.

cenpop.url  <- "https://www2.census.gov/geo/docs/reference/cenpop2020/tract/CenPop2020_Mean_TR.txt"
cenpop.file <- "data/raw/CenPop2020_Mean_TR.txt"
if (!file.exists(cenpop.file)) download.file(cenpop.url, cenpop.file, mode = "wb")

tract.centers <- read_csv(cenpop.file, col_types = cols(.default = "c")) |>
  clean_names() |>
  filter(statefp == "42") |>
  transmute(GEOID      = paste0(statefp, countyfp, tractce),
            center.pop = as.numeric(population),
            lat        = as.numeric(latitude),
            lon        = as.numeric(longitude))

message(sum(pa.tracts.tiger$GEOID %in% tract.centers$GEOID), " of ",
        nrow(pa.tracts.tiger), " PA tracts matched to a population center")

# Points in the same order as pa.tracts.tiger; tracts with no population (or no
# center) fall back to a point inside the tract polygon
pa.tract.points <- pa.tracts.tiger |>
  select(GEOID) |>
  left_join(tract.centers, by = "GEOID") |>
  mutate(use.center = !is.na(lat) & !is.na(lon) & coalesce(center.pop, 0) > 0)

center.pts <- pa.tract.points |>
  st_drop_geometry() |>
  mutate(lat = if_else(use.center, lat, 0), lon = if_else(use.center, lon, 0)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs.analysis)

fallback.pts <- st_point_on_surface(st_geometry(pa.tracts.tiger))
pt.geom <- st_geometry(center.pts)
pt.geom[!pa.tract.points$use.center] <- fallback.pts[!pa.tract.points$use.center]
pa.tract.points <- st_sf(GEOID = pa.tract.points$GEOID,
                         use.center = pa.tract.points$use.center,
                         geometry = pt.geom)
message(sum(!pa.tract.points$use.center),
        " tracts use a point inside the tract instead (no population center)")

dist_to_nearest_km <- function(from, to) {
  as.numeric(st_distance(from, to[st_nearest_feature(from, to), ], by_element = TRUE)) / 1000
}

tract.superfund <- tibble(
  GEOID                = pa.tracts.tiger$GEOID,
  # Main measure: from the tract's population center
  npl.tce.dist.km      = dist_to_nearest_km(pa.tract.points, tce.sites),
  npl.tce.n.2km        = lengths(st_is_within_distance(pa.tract.points, tce.sites, buffer.m)),
  # Previous measure (from the tract boundary), kept for comparison / sensitivity
  npl.tce.dist.edge.km = dist_to_nearest_km(pa.tracts.tiger, tce.sites)
)

# Groundwater-only version (vapor intrusion pathway)
if (any(tce.sites$groundwater)) {
  gw.sites <- tce.sites |> filter(groundwater)
  tract.superfund$npl.tce.gw.dist.km <- dist_to_nearest_km(pa.tract.points, gw.sites)
}

# How much did switching to population centers change things?
tract.superfund |>
  summarize(n.edge.zero   = sum(npl.tce.dist.edge.km == 0),
            n.center.lt2  = sum(npl.tce.dist.km < 2),
            median.change = median(npl.tce.dist.km - npl.tce.dist.edge.km),
            spearman      = cor(npl.tce.dist.km, npl.tce.dist.edge.km, method = "spearman")) |>
  print()

summary(tract.superfund)
write_csv(tract.superfund, file.path(paths$out.dir, "tract_superfund_tce.csv"))
write_csv(st_drop_geometry(tce.sites), file.path(paths$out.dir, "superfund_tce_sites.csv"))

map.superfund <- pa.tracts.map |>
  left_join(tract.superfund, by = "GEOID") |>
  mutate(dist.band = cut(npl.tce.dist.km, breaks = c(-Inf, 1, 2, 5, 10, Inf),
                         labels = c("< 1 km", "1-2 km", "2-5 km", "5-10 km", "> 10 km")))

p1 <- ggplot() +
  geom_sf(data = map.superfund, aes(fill = dist.band), lwd = 0) +
  geom_sf(data = tce.sites, color = "black", size = 0.8) +
  scale_fill_viridis_d(name = "Distance from tract\npopulation center to\nTCE Superfund site", direction = -1) +
  coord_sf(xlim = st_bbox(pa.boundary)[c(1, 3)], ylim = st_bbox(pa.boundary)[c(2, 4)]) +
  theme_void() +
  ggtitle("Distance from each tract to the nearest Superfund site with TCE")
p1
ggsave(file.path(paths$fig.dir, "superfund_tce_distance.png"), p1, width = 9, height = 6, dpi = 300)


# =============================================================================
# PART 2: AIRTOXSCREEN
# =============================================================================

# ---- 2a. Read and filter to TCE ---------------------------------------------

# Download from paths$airtox.url if given (unzipping if needed)
if (!is.null(paths$airtox.url)) {
  dl.file <- file.path("data/raw", basename(sub("\\?.*$", "", paths$airtox.url)))
  if (!file.exists(dl.file)) {
    options(timeout = max(1200, getOption("timeout")))   # default 60 s is too short
    tmp.file <- paste0(dl.file, ".part")
    ok <- tryCatch(download.file(paths$airtox.url, tmp.file, mode = "wb") == 0,
                   error = function(e) { message(conditionMessage(e)); FALSE })
    if (!ok) {
      unlink(tmp.file)
      stop("Download failed. Rerun to try again, or download the file in a browser ",
           "and save it as ", dl.file)
    }
    file.rename(tmp.file, dl.file)   # only a complete download gets the real name
  }
  if (tolower(tools::file_ext(dl.file)) == "zip") {
    unzipped <- unzip(dl.file, exdir = "data/raw/airtox")
    data.files <- unzipped[str_detect(tolower(unzipped), "\\.(csv|xlsx|xls|txt)$")]
    if (length(data.files) == 0)
      stop("The zip contains no csv/xlsx file: ", paste(basename(unzipped), collapse = ", "))
    paths$airtox <- data.files[1]
  } else {
    paths$airtox <- dl.file
  }
  message("Using AirToxScreen file: ", paths$airtox)
}

airtox <- read_any(paths$airtox, skip = skip$airtox)
glimpse(airtox)

col.id    <- find_col(airtox, cols$airtox.id,
                      c("^block$", "block_fips", "block", "^fips$", "fips", "geoid", "tract"),
                      "Block/tract ID")
col.poll  <- find_col(airtox, cols$airtox.pollutant, c("pollutant", "hap", "chemical"),
                      "Pollutant", required = FALSE)
col.value <- find_col(airtox, cols$airtox.value,
                      c("^total_conc", "^total$", "total", "exp.*conc", "amb.*conc", "conc"),
                      "Concentration (should be the TOTAL across source groups)")

if (!is.null(col.poll)) {
  airtox |> filter(str_detect(.data[[col.poll]], tce.pattern)) |>
    count(.data[[col.poll]]) |> print()
  airtox <- airtox |> filter(str_detect(.data[[col.poll]], tce.pattern))
}

airtox.tce <- airtox |>
  transmute(id = norm_id(.data[[col.id]]),
            value = suppressWarnings(as.numeric(.data[[col.value]]))) |>
  # Drop summary rows (e.g., "Entire State") and anything outside PA
  filter(str_detect(id, "^\\d+$"), !is.na(value)) |>
  mutate(id = str_pad(id, width = if_else(nchar(id) > 11, 15, 11), side = "left", pad = "0")) |>
  filter(str_starts(id, "42"))

id.lengths <- table(nchar(airtox.tce$id))
print(id.lengths)
summary(airtox.tce$value)


# ---- 2b. Aggregate to tracts ------------------------------------------------

if (all(nchar(airtox.tce$id) == 11)) {
  # Already tract-level (older AirToxScreen years)
  tract.airtox <- airtox.tce |> transmute(GEOID = id, airtox.tce = value)

} else {
  # Block-level (2020+): population-weighted mean within each tract
  pop.file <- file.path(paths$out.dir, "pa_block_pop_2020.rds")
  if (file.exists(pop.file)) {
    block.pop <- readRDS(pop.file)
  } else {
    pa.counties <- fips_codes |> filter(state == "PA") |> pull(county_code)
    block.pop <- map(pa.counties, \(cty)
      get_decennial(geography = "block", variables = "P1_001N", year = 2020,
                    sumfile = "pl", state = "PA", county = cty)) |>
      list_rbind() |>
      select(id = GEOID, pop = value)
    saveRDS(block.pop, pop.file)
  }

  airtox.blocks <- airtox.tce |>
    filter(nchar(id) == 15) |>
    left_join(block.pop, by = "id")
  message(round(100 * mean(!is.na(airtox.blocks$pop)), 1),
          "% of AirToxScreen blocks matched to 2020 block population")

  tract.airtox <- airtox.blocks |>
    mutate(GEOID = str_sub(id, 1, 11), pop = replace_na(pop, 0)) |>
    group_by(GEOID) |>
    summarize(n.blocks   = n(),
              airtox.pop = sum(pop),
              airtox.tce = if (sum(pop) > 0) weighted.mean(value, pop) else mean(value),
              .groups = "drop")
}

message(sum(!pa.tracts.tiger$GEOID %in% tract.airtox$GEOID),
        " PA tracts have no AirToxScreen value (often unpopulated tracts)")
summary(tract.airtox$airtox.tce)
write_csv(tract.airtox, file.path(paths$out.dir, "tract_airtox_tce.csv"))

p2 <- ggplot() +
  geom_sf(data = left_join(pa.tracts.map, tract.airtox, by = "GEOID"),
          aes(fill = airtox.tce), lwd = 0) +
  scale_fill_viridis_c(name = "TCE (µg/m³)", trans = "log10", na.value = "grey85") +
  theme_void() +
  ggtitle("AirToxScreen 2020 modeled TCE concentration, PA tracts")
p2
ggsave(file.path(paths$fig.dir, "airtox_tce.png"), p2, width = 9, height = 6, dpi = 300)


# =============================================================================
# PART 3: COMBINE
# =============================================================================

tract.combined <- tibble(GEOID = pa.tracts.tiger$GEOID) |>
  left_join(tract.superfund, by = "GEOID") |>
  left_join(select(tract.airtox, GEOID, airtox.tce), by = "GEOID")

# Do the two exposure sources agree? (low correlation is expected: different pathways)
cor(tract.combined$npl.tce.dist.km, tract.combined$airtox.tce,
    method = "spearman", use = "complete.obs")

write_csv(tract.combined, file.path(paths$out.dir, "tract_tce_superfund_airtox.csv"))
message("Done. Wrote ", nrow(tract.combined), " tracts to ",
        file.path(paths$out.dir, "tract_tce_superfund_airtox.csv"))
