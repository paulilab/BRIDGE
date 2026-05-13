# BRIDGE unified color palette
# All plot colors should be drawn from these constants so the appearance
# stays consistent across volcano, PCA, datapoints, heatmap and integration
# plots.

# ---------------------------------------------------------------------------
# Named constants ----
# ---------------------------------------------------------------------------

#' BRIDGE color palette
#'
#' A named list of color constants used across all BRIDGE plot types.
#'
#' \describe{
#'   \item{accent}{Primary brand blue used for single-color geometry (scatter
#'     points, trend lines, spinners).}
#'   \item{up}{Color for up-regulated features in volcano plots.}
#'   \item{down}{Color for down-regulated features in volcano plots.}
#'   \item{ns}{Color for non-significant features.}
#'   \item{highlight}{Color for user-highlighted genes.}
#'   \item{discrete}{Ordered vector of colors for categorical mappings (groups,
#'     genes). Interpolated automatically when more levels are needed.}
#'   \item{heatmap_low / heatmap_mid / heatmap_high}{Endpoints for the
#'     diverging expression heatmap scale.}
#'   \item{trend_low / trend_mid / trend_high}{Endpoints for the fuzzy-cluster
#'     membership score color ramp (integration trend lines).}
#' }
#'
#' @export
BRIDGE_COLORS <- list(
    # Single-geometry accent (replaces the old #2b8cbe)
    accent       = "#4C72B0",

    # Volcano / differential expression directions
    up           = "#DD4949",   # coral-red
    down         = "#4C72B0",   # indigo-blue
    ns           = "#CCCCCC",   # light gray
    highlight    = "#2C3E50",   # near-black

    # Discrete categorical palette (8 base colors, interpolated for larger sets)
    # Inspired by the Seaborn "muted" set; verified colorblind-distinguishable.
    discrete = c(
        "#4C72B0",  # indigo blue
        "#DD4949",  # coral red
        "#55A868",  # sage green
        "#C44E52",  # brick
        "#8172B2",  # lavender
        "#937860",  # sand
        "#DA8BC3",  # pink
        "#8C8C8C"   # medium gray
    ),

    # Diverging scale for expression heatmaps (blue → white → red)
    heatmap_low  = "#4C72B0",
    heatmap_mid  = "#FFFFFF",
    heatmap_high = "#DD4949",

    # Membership score ramp for fuzzy-cluster trend lines
    trend_low    = "#55A868",   # green  → low membership
    trend_mid    = "#F0C929",   # yellow → mid
    trend_high   = "#DD4949"    # red    → high
)

# ---------------------------------------------------------------------------
# Helper functions ----
# ---------------------------------------------------------------------------

#' Generate a discrete color vector of length \emph{n}
#'
#' Returns exactly \code{n} colors drawn from \code{BRIDGE_COLORS$discrete}.
#' When \code{n} exceeds the 8 base colors the palette is interpolated via
#' \code{\link[grDevices]{colorRampPalette}}.
#'
#' @param n Integer. Number of colors required.
#' @return Character vector of hex color codes.
#' @export
bridge_discrete_pal <- function(n) {
    base <- BRIDGE_COLORS$discrete
    if (n <= length(base)) return(base[seq_len(n)])
    grDevices::colorRampPalette(base)(n)
}

#' Diverging heatmap color ramp (circlize)
#'
#' Constructs a \code{circlize::colorRamp2} object for use as the \code{col}
#' argument in \code{ComplexHeatmap::Heatmap}.
#'
#' @param col_limit Numeric. Upper limit of the scale. When \code{col_min} is
#'   \code{NULL} the scale is symmetric: \code{[-col_limit, 0, col_limit]}.
#'   Defaults to \code{2}.
#' @param n_steps Integer. Number of interpolation steps (must be odd so the
#'   midpoint is the mid colour).  Defaults to \code{11}.
#' @param col_min Numeric. Optional lower limit. When supplied the scale runs
#'   \code{[col_min, mid, col_limit]} instead of being symmetric.
#' @return A \code{colorRamp2} function.
#' @export
bridge_heatmap_col <- function(col_limit = 2, n_steps = 11, col_min = NULL) {
    stopifnot(n_steps %% 2 == 1)
    lower <- if (is.null(col_min)) -col_limit else col_min
    breaks <- seq(lower, col_limit, length.out = n_steps)
    mid    <- ceiling(n_steps / 2)
    colors <- c(
        grDevices::colorRampPalette(
            c(BRIDGE_COLORS$heatmap_low, BRIDGE_COLORS$heatmap_mid)
        )(mid),
        grDevices::colorRampPalette(
            c(BRIDGE_COLORS$heatmap_mid, BRIDGE_COLORS$heatmap_high)
        )(mid)[-1]
    )
    circlize::colorRamp2(breaks, colors)
}

#' Membership-score trend-line color ramp (circlize)
#'
#' Used to color fuzzy-cluster trend lines by their membership score.
#'
#' @return A \code{colorRamp2} function mapping \code{[0.01, 0.5, 0.9]} to
#'   the low / mid / high trend colors.
#' @export
bridge_trend_col <- function() {
    circlize::colorRamp2(
        c(0.01, 0.5, 0.9),
        c(BRIDGE_COLORS$trend_low, BRIDGE_COLORS$trend_mid, BRIDGE_COLORS$trend_high)
    )
}
