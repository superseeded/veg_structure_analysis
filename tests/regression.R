# Run from the repository root: Rscript tests/regression.R
# Load the actual function definitions without attaching optional web/export
# packages. Numerical/spatial tests need sf; the optional plot check needs ggplot2.
library(sf)
utils <- new.env(parent = globalenv())
for (expr in parse("R/utils.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
    eval(expr, utils)
  }
}
attach(utils, name = "vegetation_regression_functions")
checks <- 0L
check <- function(name, expression) {
  force(expression)
  checks <<- checks + 1L
  cat("PASS:", name, "\n")
}
equal <- function(actual, expected, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = tolerance, check.attributes = FALSE)))
}
errors <- function(expression) stopifnot(inherits(try(force(expression), silent = TRUE), "try-error"))

check("zero-safe evenness and represented-category convention", {
  equal(shannon_evenness(c(10, 10, 0)), 1)
  equal(shannon_evenness(c(10, 0, 0)), 0)
  equal(shannon_evenness(c(0, 0)), 0)
  equal(shannon_evenness(numeric()), 0)
  errors(shannon_evenness(c(1, NA)))
  errors(shannon_evenness(c(1, -1)))
})
check("entropy retains precision until final score rounding", {
  proportions <- c(1, 2, 9) / 12
  equal(shannon_evenness(c(1, 2, 9)), -sum(proportions * log(proportions)) / log(3))
})
check("explicit native labels and one-class inputs", {
  equal(native_proportion(c("Native", " non-native ")), .5)
  equal(native_proportion(rep("native", 4)), 1)
  equal(native_proportion(rep("introduced", 4)), 0)
  equal(native_proportion(c("exotic", "native", "naturalised")), 1/3)
  errors(native_proportion(c("native", "unknown")))
  errors(native_proportion(c("native", NA)))
})
check("mounding radii are dimensionless and scale-invariant", {
  for (height in c(.25, 1, 4, 10)) {
    equal(mounding(0, height), 1)
    equal(mounding(height * .6, height), .8)
    equal(mounding(height, height), 0)
  }
  equal(mounding(c(0, .6, 1), 1), c(1, .8, 0))
  equal(mounding(c(0, 0), c(0, 4)), c(0, 1))
  equal(mounding(0, 0), 0)
  equal(mounding(2, 1), 0)
})
check("upright foliage is finite at its base", {
  equal(upright(c(0, 2, 4), 4), c(1, 1, 1))
  equal(upright(0, 4), 1)
  equal(upright(4, 4), 1)
  equal(upright(0, 0), 0)
  equal(upright(5, 4), 0)
})

boundary <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0)
))), crs = 3857))
grid <- create_grid(boundary, cellarea = 10)
grid$prop_ol <- .2
year_two <- grid
year_two$prop_ol <- .8
spatial <- list(year_01 = list(cut_0100 = grid), year_02 = list(cut_0100 = year_two))
attr(spatial, "gridspacing") <- sqrt(2 * 10 / (3 * sqrt(3))) * sqrt(3)
plants <- st_as_sf(data.frame(
  species = paste("Species", letters[1:12]), density = rep(c("open", "semi-closed", "closed"), 4),
  texture = rep(c("fine", "medium", "coarse"), 4), max_height = rep(c(.5, 1, 2, 4), 3),
  endemism = rep(c("native", "exotic"), 6), phenology = "jan feb mar apr may jun",
  x = seq_len(12) / 2, y = 5
), coords = c("x", "y"), crs = 3857)
results <- suppressWarnings(analyse_spatial_data(plants, spatial))
check("equal category counts receive full balance scores", {
  equal(results[[1]]$density_score, 100)
  equal(results[[1]]$size_score, 100)
  equal(results[[1]]$texture_score, 100)
  equal(results[[1]]$phenology_score, 50)
})
check("native rings always contain 100 ticks", {
  equal(nrow(results[[1]]$endemism_data), 100)
  equal(results[[1]]$endemism_score, 50)
  introduced <- plants
  introduced$endemism <- "introduced"
  zero <- suppressWarnings(analyse_spatial_data(introduced, spatial))[[1]]
  equal(zero$endemism_score, 0)
  equal(nrow(zero$endemism_data), 100)
  equal(sum(zero$endemism_data$values), 0)
})
check("period scores use all years, not the last year", {
  equal(period_metric_score(results, "coverage"), 50)
  equal(period_metric_score(results, "connectivity"), 50)
  equal(period_metric_score(results[1], "coverage"), 20)
  equal(period_metric_score(results[1], "connectivity"), 0)
})
check("overall score is independent of global score variables", {
  before <- scorecard_scores(results)
  assign("unrelated_score", -10000, envir = .GlobalEnv)
  assign("coverage_score", 999, envir = .GlobalEnv)
  equal(scorecard_scores(results), before)
  stopifnot(length(before) == 8, all(is.finite(before)))
  rm(unrelated_score, coverage_score, envir = .GlobalEnv)
})
check("both connectivity entry points honour the same threshold", {
  for (threshold in c(.1, .3, .9)) {
    analysed <- suppressWarnings(analyse_spatial_data(plants, spatial, threshold = threshold))
    separate <- suppressWarnings(estimate_connectivity(spatial, threshold = threshold))
    equal(vapply(analysed, function(year) year$connectivity_data[[1]], numeric(1)), separate[[1]])
  }
  errors(analyse_spatial_data(plants, spatial, threshold = -1))
  errors(estimate_connectivity(spatial, threshold = NA_real_))
})
check("an isolated cell produces finite zero neighbour connectivity", {
  single <- grid[1, ]
  single$id <- 1
  single$prop_ol <- 1
  isolated <- list(year_01 = list(cut_0100 = single))
  attr(isolated, "gridspacing") <- attr(spatial, "gridspacing")
  equal(suppressWarnings(analyse_spatial_data(plants, isolated))[[1]]$connectivity_score, 0)
  equal(suppressWarnings(estimate_connectivity(isolated))[[1]], 0)
})
check("polygon attributes are validated separately from point attributes", {
  required <- c("coverage", "density", "endemism", "form", "ini_height", "max_height",
                "max_width", "phenology", "ref_height", "spacing", "species", "texture", "year_max")
  valid <- st_sf(as.data.frame(setNames(rep(list(1), length(required)), required)),
                 geometry = st_sfc(st_point(c(1, 1)), crs = 3857))
  check_spatial_input(valid)
  errors(check_spatial_input(valid, boundary))
})
check("stored study scorecards average their complete time series", {
  for (site in c(torquay = "park_tq", averley = "park_av", booyeembara = "park_b")) {
    folder <- c(park_tq = "torquay", park_av = "averley", park_b = "booyeembara")[[site]]
    saved <- new.env()
    load(file.path("output", folder, paste0(site, "_results")), envir = saved)
    years <- saved[[ls(saved)[1]]]
    equal(period_metric_score(years, "coverage"),
          100 * mean(vapply(years, function(y) y$coverage_data, numeric(1))))
    equal(period_metric_score(years, "connectivity"),
          100 * mean(vapply(years, function(y) mean(y$connectivity_data), numeric(1))))
    stopifnot(length(scorecard_scores(years)) == 8)
  }
})
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  check("coverage chart shows the period mean and does not create globals", {
    plot <- plot_circ_bar(results, "coverage", colours = "#c9837d")
    labels <- unlist(lapply(ggplot_build(plot)$data, function(layer) layer$label))
    stopifnot("50/100\nCOVERAGE\nSCORE" %in% labels)
    stopifnot(!exists("coverage_score", envir = .GlobalEnv, inherits = FALSE))
  })
} else cat("SKIP: coverage plot rendering (ggplot2 unavailable)\n")
cat(checks, "regression checks passed.\n")
