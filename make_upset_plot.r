


### Make plots of self-reported diagnoses by class membership
### Upset plot and donut chart
### 31/02/21

library(ggplot2)
library(ComplexUpset)
library(janitor)
library(dplyr)
library(grid)
library(gridExtra)
library(png)
library(cowplot)

output_dir <- "~/brc_scratch/output/manic_lca"


######################## EXAMPLE DATA #######################

# movies <- as.data.frame(ggplot2movies::movies)
# genres <- colnames(movies)[18:24]

# for simplicity of examples, only use the complete data points
# movies[movies$mpaa == "", "mpaa"] <- NA
# movies <- na.omit(movies)

# upset1 <- upset(
#  movies,
#  genres,
#  base_annotations = list(
#    "Intersection size" = intersection_size(
#      counts = F,
#      aes = aes(fill = mpaa)
#    )
#  ),
#  width_ratio = 0.1
# )

# bitmap(file = paste0(output_dir, "/", "upset_tmp_2.png"),
#  width = 26, height = 11, type = "png16m", units = "cm", res = 300
# )

# print(upset1)

# dev.off()


###

working_df <- readRDS(paste0(output_dir, "/working_df_diagnoses_for_plots_tmp.rds"))

names(working_df)[names(working_df) == "predclass"] <- "Classes"

disorders_ls <- names(working_df[1:6])

working_df_indx <- working_df[, c(disorders_ls)]

### split data in to "no diagnoses" and ">1 diagnosis" sets

working_df_cc <- working_df[which(rowSums(working_df_indx[, -1]) > 0), ]

working_df_zeros <- working_df[which(rowSums(working_df_indx[, -1]) == 0), ]

### basic upset plot

upset_basic <- upset(working_df_cc, disorders_ls,
  base_annotations = list("Intersection size" = intersection_size(counts = F, mapping = aes(fill = Classes)) + scale_fill_manual(values = c(
    "1" = "#E41A1C", "2" = "#377EB8",
    "4" = "#984EA3", "5" = "#FF7F00", "3" = "#4DAF4A"
  ))),
  width_ratio = 0.2, name = "Group (diagnosis combination)"
)

# red "#E41A1C" blue "#377EB8" green "#4DAF4A" purple "#984EA3" orange "#FF7F00" - Set1

bitmap(
  file = paste0(output_dir, "/", "upset_plot_ukb_basic.png"),
  width = 32, height = 14, type = "png16m", units = "cm", res = 300
)

print(upset_basic)

dev.off()

### more complex upset plot

rating_scale <- scale_fill_manual(values = c(
  "1" = "#E41A1C", "2" = "#377EB8",
  "4" = "#984EA3", "5" = "#FF7F00", "3" = "#4DAF4A"
))

show_hide_scale <- scale_color_manual(values = c("show" = "black", "hide" = "transparent"), guide = FALSE)

upset3 <- upset(
  working_df_cc, disorders_ls, wrap=T ,
  set_sizes = (upset_set_size(position = "right")), name = "Group (diagnosis combination)",
  width_ratio = 0.15, min_size = 2, base_annotations = list("Intersection size" = intersection_size(
    text = list(
      vjust = -0.65,
      hjust = 0.3,
      angle = 0, color = "#E41A1C"
    )
  )),
  annotations = list(
    "Class\n (percent of intersection)" = list(
      aes = aes(x = intersection, fill = Classes),
      geom = list(
        geom_bar(stat = "count", position = "fill", na.rm = TRUE),
        geom_text(
          aes(
            label = !!aes_percentage(relative_to = "intersection"),
            color = "hide"
          ),
          stat = "count",
          position = position_fill(vjust = 0.5)
        ),
        scale_y_continuous(labels = scales::percent_format()), show_hide_scale,
        rating_scale
      )
    )
  ), guides = "over"
)

bitmap(
  file = paste0(output_dir, "/", "upset_plot_ukb_complex.png"),
  width = 32, height = 14, type = "png16m", units = "cm", res = 300
)

print(upset3)

dev.off()

### donut chart

count.data <- tabyl(working_df_zeros$Classes)

names(count.data) <- c("Classes", "N", "Percent")

count.data$lab <- paste0(count.data$N, "(", signif(count.data$Percent * 100, 2), "%)")

count.data <- count.data %>%
  arrange(desc(Classes)) %>%
  mutate(lab.ypos = cumsum(Percent) - 0.5 * Percent)

print(count.data)

donut1 <- ggplot(count.data, aes(x = 2, y = Percent, fill = Classes)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(x = 3.0, aes(y = lab.ypos, label = lab), size = 5, color = "grey10") +
  scale_fill_manual(values = c(
    "1" = "#E41A1C", "2" = "#377EB8",
    "4" = "#984EA3", "5" = "#FF7F00", "3" = "#4DAF4A"
  )) +
  theme_void() +
  xlim(0.1, 2.9) +
  ggtitle(paste0("No diagnoses\n (N = ", sum(count.data$N), ")")) +
  theme(plot.title = element_text(hjust = 0.5))


bitmap(
  file = paste0(output_dir, "/", "donut_chart_no_diagnoses_ukb.png"),
  width = 24, height = 20, type = "png16m", units = "cm", res = 300
)

print(donut1)

dev.off()

### combine plots in to panel plot

png1_dest <- paste0(output_dir, "/donut_chart_no_diagnoses_ukb.png")

png2_dest <- paste0(output_dir, "/upset_plot_ukb_complex.png")

img2 <- grid::rasterGrob(as.raster(readPNG(png2_dest)),
                         interpolate = FALSE)
img1 <- grid::rasterGrob(as.raster(readPNG(png1_dest)),
                         interpolate = FALSE)

bitmap(
  file = paste0(output_dir, "/", "combined_upset_donut_plot_ukb.png"),
  width = 30, height = 10, type = "png16m", units = "cm", res = 300
)

	print(plot_grid(img1, img2, labels = c("A", "B"), ncol = 2, greedy = F, label_size = 16, hjust = -1, scale = c(0.8, 1)))

dev.off()

#
#
