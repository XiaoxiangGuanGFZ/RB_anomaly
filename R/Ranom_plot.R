


#' Plot a Record-Breaking Anomaly Series
#' 
#' Draw anomaly bars and a horizontal zero reference.
#' 
#' @param Ranom Numeric vector of record-breaking anomalies in percent.
#' @param year Numeric observation years, one per input column, in temporal order.
#' @param p_print Logical; print the plot when TRUE. The plot is always returned.
#' @return A ggplot object.
#' @export
plot_Ranom <- function(Ranom, year, p_print = TRUE) {

  df = data.frame(
    time = year,
    value = Ranom
  )
  p <- ggplot(df, aes(x = time)) +

    # observed RB anomaly bars
    geom_col(
      aes(y = value),
      fill = "grey70",
      color = "grey45",
      linewidth = 0.2,
      width = 0.85
    ) +

    # zero reference line
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      color = "black"
    ) +

    labs(
      x = "year",
      y = "Record-breaking anomaly (%)"
    ) +

    theme_minimal()

  if (p_print) {
    print(p)
  }

  return(p)
}




#' Plot Anomalies, Trends, and Confidence Intervals
#' 
#' Draw observed anomalies, the supplied trend, and 90 and 95 percent interval ribbons. Accepts both LOESS and SSA results with the documented columns.
#' 
#' @param df_Ranom Result from [rb_anomaly_ci()], its data component, or a
#'   legacy data frame with time, value, loess_trend and interval columns.
#'   Legacy 90-percent-only LOESS results are also supported. New results plot
#'   statistic_value against matching bands; legacy results keep their trend.
#' @param p_print Logical; print the plot when TRUE. The plot is always returned.
#' @return A ggplot object.
#' @export
plot_Ranom_CI <- function(df_Ranom, p_print = TRUE) {

    # columns in the input data frame df_Ranom:
    # time, value, loess_trend,
    # ci_lower_90, ci_upper_90,
    # ci_lower_95, ci_upper_95

    if (is.list(df_Ranom) && !is.data.frame(df_Ranom) &&
        !is.null(df_Ranom$data)) df_Ranom <- df_Ranom$data
    df <- df_Ranom
    if ("statistic_value" %in% names(df)) {
      df$value <- df$Ranom
      df$loess_trend <- df$statistic_value
    }
    if (all(c("ci_lower", "ci_upper") %in% names(df))) {
      df$ci_lower_90 <- df$ci_lower
      df$ci_upper_90 <- df$ci_upper
      df$ci_lower_95 <- NA_real_
      df$ci_upper_95 <- NA_real_
    }
    p <- ggplot(df, aes(x = time)) +

        # observed RB anomaly bars
        geom_col(
          aes(y = value),
          fill = "grey70",
          color = "grey45",
          linewidth = 0.2,
          width = 0.85
        ) +

        # 95% reference interval: wider / outer band
        geom_ribbon(
            aes(ymin = ci_lower_95, ymax = ci_upper_95),
            fill = "#4C78A8",
            alpha = 0.35
        ) +

        # 90% reference interval: narrower / inner band
        geom_ribbon(
            aes(ymin = ci_lower_90, ymax = ci_upper_90),
            fill = "#9ECAE1",
            alpha = 0.45
        ) +

        # zero reference line
        geom_hline(
            yintercept = 0,
            linewidth = 0.4,
            color = "black"
        ) +


        # smoothed / long-term trend
        geom_line(
            aes(y = loess_trend),
            color = "black",
            linewidth = 1.1
        ) +

        labs(
            x = "year",
            y = "Record-breaking anomaly (%)"
        ) +

        theme_minimal()

    if (p_print) {
      print(p)
    }

    return(p)
}



#' Format Trend Estimates
#'
#' A formattable formatter for HTML tables.
#' @param x Numeric trend estimates to format.
#' @param data Optional data frame used when evaluating formatting expressions.
#' @param text Optional display text corresponding to x.
#' @format A formatter function.
#' @export
trend_formatter <- formattable::formatter(
  "span",
  style = x ~ formattable::style(
    color = ifelse(x > 0, "red", "green"), # Green for positive, red for negative
    font.weight = "bold"
  ),
  x ~ sprintf("%+.3f", x) # Format to include + or - sign with one decimal
)

#' Format P-Values
#'
#' A formattable formatter for HTML tables.
#' @param x Numeric p-values to format.
#' @param data Optional data frame used when evaluating formatting expressions.
#' @param text Optional display text corresponding to x.
#' @format A formatter function.
#' @export
pvalue_formatter <- formattable::formatter(
  "span",
  style = x ~ formattable::style(
    background.color = ifelse(x < 0.05, "gray", "transparent"), # Gray background for p < 0.05
    color = ifelse(x < 0.05, "white", "black") # White text if background is gray for contrast
  ),
  x ~ sprintf("%.3f", x) # Format p-values to three decimal places
)

