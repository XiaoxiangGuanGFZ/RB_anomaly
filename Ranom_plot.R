
library(ggplot)


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




plot_Ranom_CI <- function(df_Ranom, p_print = TRUE) {

    # columns in the input data frame df_Ranom:
    # time, value, loess_trend,
    # ci_lower_90, ci_upper_90,
    # ci_lower_95, ci_upper_95

    df = df_Ranom
    p <- ggplot(df, aes(x = time)) +

        # observed RB anomaly bars
        geom_col(
          aes(y = value),
          fill = "grey70",
          color = "grey45",
          linewidth = 0.2,
          width = 0.85
        ) +

        # 95% confidence interval: darker / inner band
        geom_ribbon(
            aes(ymin = ci_lower_95, ymax = ci_upper_95),
            fill = "#4C78A8",
            alpha = 0.35
        ) +

        # 90% confidence interval: lighter / outer band
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


