#' Reconstruct an SSA Trend
#'
#' Reconstruct selected singular spectrum analysis components. Leading and
#' trailing missing values (including the undefined first record anomaly) are
#' excluded from fitting and restored in the result. Interior missing values
#' are rejected rather than interpolated. All-missing series return all NA.
#' @param x Numeric time series.
#' @param L Integer SSA window length, at least 2 and smaller than the number
#'   of observations remaining after removing missing edges.
#' @param groups Named list of eigentriple indices with a component named trend.
#' @return Numeric vector of the same length as x.
#' @export
#' @examples
#' rb_ssa_trend(c(NA, seq_len(20)), L = 5)
rb_ssa_trend <- function(x, L = 15, groups = list(trend = 1)) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop("x must be a numeric vector.")
  }
  out <- rep(NA_real_, length(x))
  valid <- which(!is.na(x))
  if (!length(valid)) return(out)
  index <- seq.int(min(valid), max(valid))
  values <- x[index]
  if (any(!is.finite(values))) {
    stop("SSA requires finite observations without interior missing values.")
  }
  if (length(L) != 1L || !is.finite(L) || L != as.integer(L) ||
      L < 2L || L >= length(values)) {
    stop("L must be an integer >= 2 and smaller than the usable series length.")
  }
  if (!is.list(groups) || is.null(groups$trend)) {
    stop("groups must be a list with a component named trend.")
  }
  s <- Rssa::ssa(values, L = L, kind = "1d-ssa")
  r <- Rssa::reconstruct(s, groups = groups)
  out[index] <- as.numeric(r$trend)
  out
}

#' Reconstruct SSA Trends for Matrix Rows
#'
#' Apply [rb_ssa_trend()] independently to each row while preserving matrix shape.
#' @param data_matrix Numeric matrix with series in rows and time in columns.
#' @inheritParams rb_ssa_trend
#' @return Numeric matrix with the same dimensions and dimnames as data_matrix.
#' @export
rb_ssa_trend_matrix <- function(data_matrix, L = 15, groups = list(trend = 1)) {
  if (!is.matrix(data_matrix) || !is.numeric(data_matrix)) {
    stop("data_matrix must be a numeric matrix.")
  }
  out <- matrix(NA_real_, nrow(data_matrix), ncol(data_matrix),
                dimnames = dimnames(data_matrix))
  for (i in seq_len(nrow(data_matrix))) {
    out[i, ] <- rb_ssa_trend(data_matrix[i, ], L = L, groups = groups)
  }
  out
}

# ---- Compatibility aliases ----------------------------------------------
# Older scripts can keep their names and named arguments. These aliases
# refer to the same implementation; use the rb_* API for new analyses.

#' @rdname rb_ssa_trend
#' @details `ssa_trend()` is a compatibility alias. Use `rb_ssa_trend()`
#'   in new code; arguments and returned values are unchanged.
#' @export
ssa_trend <- rb_ssa_trend

#' @rdname rb_ssa_trend_matrix
#' @details `ssa_trend_matrix()` is a compatibility alias. Use `rb_ssa_trend_matrix()`
#'   in new code; arguments and returned values are unchanged.
#' @export
ssa_trend_matrix <- rb_ssa_trend_matrix
