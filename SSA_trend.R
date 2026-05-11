
library(Rssa)


ssa_trend <- function(x, L = 15, groups = list(trend = 1)) {
  # x: time series vector
  # L: SSA window length
  # groups: SSA eigentriples used for trend reconstruction
  
  s <- ssa(x, L = L, kind = "1d-ssa")
  r <- reconstruct(s, groups = groups)
  
  return(as.numeric(r$trend))
}

ssa_trend_matrix <- function(data_matrix, L = 15, groups = list(trend = 1)) {
    data_matrix <- as.numeric(data_matrix)
    p <- dim(data_matrix)[1]  # number of rows (stations)
    n <- dim(data_matrix)[2]  # number of columns (years)
    # ssa_trend <- matrix(NA_real_, ncol = n, nrow = p)

    ssa_trend <- apply(Ranom_iid, 1, function(x) {
        ssa_trend(x, L = L, groups = groups)
    })
    return(ssa_trend)
}







