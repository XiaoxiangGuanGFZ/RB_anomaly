
# ------- set the work space and load the required libraries --------
# required packages:
# tidyverse, trend, formattable, ggplot2
require(tidyverse)
require(ggplot2)
require(trend)
require(formattable)
# library(patchwork)

ws <- "D:/HW2024_SD/example/"

# load the mega flood dataset,
# megaflood data source: https://github.com/tuwhydro/megafloods
# Two data frames in this RData file: df_AMS_All, and df_gauge
# df_AMS_All:
#  - 2D, contains the annual maximum specific discharges in m3/s/km2
#  - each column stands for one gauge
#  - first column: year
# df_gauge: meta data of this megaflood dataset
#
library(RBanomaly)
data("mega_floods", package = "RBanomaly")
dim(df_AMS_All)
dim(df_gauge)

# Algorithms are loaded through library(RBanomaly) above.


df_AMS = t(df_AMS_All) ## transpose, you can check the dimensions
# dim(df_AMS)
year =  df_AMS[1, ]    ## the temporal coverage of the megaflood dataset

df <- df_AMS[-1, ]     ## pure AMS matrix, the year column is removed
N = dim(df)[2]         ## N: length of the observation (time series)
index_von <- which(year == 1925)
index_bis <- which(year == 2021)
Ranom_data = df[, index_von:index_bis]  ## input ready to be apply
dim(Ranom_data)
Ranom_year = year[index_von:index_bis]  ## input ready to be apply
length(Ranom_year)


# filter the data based on 3 constraints, reduce the data size to accelerate the calculation
# here, a number of gauges were removed,
id_filter = df_gauge$id[df_gauge$area <= 1000 & df_gauge$regions1 == 1 & df_gauge$dams == "no"]
n_gauges = length(id_filter)
Ranom_data_filter = Ranom_data[id_filter, ]
dim(Ranom_data_filter)
n_stations = length(id_filter)

# implement the model to calculate the flood anomaly
# runs = 5000 for robust results, smaller value to test the model
df_out = LOESS_Trend_CI(Ranom_data_filter, year = Ranom_year, runs = 1000)

write.table(
  df_out,
  file = paste0(ws, "RB_flood_anomaly.csv"),
  col.names = T, row.names = F, quote = F, append = F, sep = ","
)

# -------- visualization ----------
# the columns in the output dataframe
# `value`: the flood anomaly series
# `loess_trend`: loess trend of the observed flood anomaly
# `ci_upper`, and `ci_lower`: the confidence interval of flood anomaly under a stationary environment.

# -- interpretation
# if the `loess_trend` is outside the confidence bound, we say that the observed flood anomaly is significantly
# different from that in a stationary environment, and we could say climate change increases or decreases the flood anomaly.
# if the `loess_trend` is within the confidence bound, we say it is rather stable, no significant trend is detected on flood anomaly,
# meaning that it is not quite different from that under a stationary climate condition

ggplot() +
  geom_bar(
    data = df_out, aes(x = time, y = value),
    stat = "identity", fill = "gray"
  ) +
  geom_ribbon(
    data = df_out,
    aes(x = time, ymin = ci_upper, ymax = ci_lower),
    alpha = 0.3,
    fill = 'lightblue'
  ) +
  geom_line(
    data = df_out,
    aes(x = time, y = loess_trend), color = "blue"
  ) +
  labs(y = "annual record-breaking anomaly (%)", x = "year",
       title = paste0("anomaly of record-breaking floods"),
       subtitle = paste0("number of stations: ", n_stations)) +
  theme_minimal()


df_Ranom$rb_ssa_trend = rb_ssa_trend(df_Ranom$Ranom, L = 20)

ggplot(data = df_Ranom) +
  geom_col(aes(x = year, y = Ranom)) +
  geom_line(aes(x = year, y=rb_ssa_trend), color = "navy", size = 1.2) +
  labs(y = "annual record-breaking anomaly (%)", x = "year") +
  theme_minimal()

# ---- test the SSA_trend

df <- read.csv(
  "D:/guan/projects/RB_anomaly_flood/data/Ranom_summer_small.csv"
)

L_test <- c(5, 10, 15, 20, 25, 30, 35)
L_test <- c(5, 10, 15)

for (i in seq_along(L_test)) {
  L = L_test[i]
  df$trend <- rb_ssa_trend(df$value, L = L)
  p <- ggplot(data = df) +
    geom_col(aes(x = time, y = value),
             fill = "grey70",
             color = "grey45",
             linewidth = 0.2,
             width = 0.85) +
    geom_line(aes(x = time, y = trend), color = "black", size = 1.2) +
    labs(y = "annual record-breaking anomaly (%)", x = "year") +
    theme_minimal()
  print(p)
  # ggsave(filename = paste0(ws, "results/CAMELS_DE_flood/figures/test_SSA_trend/L_", L, ".pdf"),
  #        plot = p, device = "pdf", dpi = 300, unit = "cm", width = 12, height = 10)

}
lm(data = df, value~time)
library(Rssa)

x <- df$value
L <- 15

s <- ssa(x, L = L, kind = "1d-ssa")

plot(s)                   # eigenvalues
plot(s, type = "vectors") # eigenvectors
plot(s, type = "paired")  # paired components
wcor <- wcor(s, groups = 1:10)
plot(wcor)


r1 <- reconstruct(s, groups = list(trend = 1))
r2 <- reconstruct(s, groups = list(trend = 1:2))
r3 <- reconstruct(s, groups = list(trend = 1:3))
r4 <- reconstruct(s, groups = list(trend = 1:4))
r5 <- reconstruct(s, groups = list(trend = 1:5))

df$trend_1   <- as.numeric(r1$trend)
df$trend_12  <- as.numeric(r2$trend)
df$trend_123 <- as.numeric(r3$trend)
df$trend_124 <- as.numeric(r4$trend)
df$trend_125 <- as.numeric(r5$trend)


ggplot(df, aes(x = time)) +
  geom_col(aes(y = value),
           fill = "grey70",
           color = "grey45",
           linewidth = 0.2,
           width = 0.85) +
  geom_line(aes(y = trend_1), color = "black", linewidth = 1.0) +
  geom_line(aes(y = trend_12), color = "red", linewidth = 1.0) +
  geom_line(aes(y = trend_123), color = "blue", linewidth = 1.0) +
  geom_line(aes(y = trend_124), color = "green", linewidth = 1.0) +
  labs(y = "annual record-breaking anomaly (%)",
       x = "year") +
  theme_minimal()


list_df <- list()

list_df[[1]] <- Ranom_data_filter
year = Ranom_year

results <- Func_RB_region_list_CI_fast(
    list_df,
    year,
    seed = 1,
    runs = 2000,
    L = 20,
    groups = list(trend = 1)
)

rb_plot_anomaly_ci(results)
