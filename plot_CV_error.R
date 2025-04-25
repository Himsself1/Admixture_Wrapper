# * Libraries

list_of_packages <- c(
  "ggplot2", "reshape2",
  "forcats", "ggthemes",
  "patchwork", "gridExtra",
  "grid", "argparse",
  "Cairo"
)

for (i in list_of_packages) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i)
  }
  library( i, character.only = T )
}

parser <- ArgumentParser()

parser$add_argument("-plot_folder",
  default = T, action = "store", type = "character",
  help = "Folder of output plots"
)

parser$add_argument("-input_file",
  default = T, action = "store", type = "character",
  help = "Input of CV errors"
)

parser$add_argument("-name",
  default = T, action = "store", type = "character",
  help = "Name of the run"
)


arguments <- parser$parse_args()
if (!dir.exists(arguments$plot_folder)) {
  dir.create(arguments$plot_folder, recursive = T)
  print("Created plot folder and parents")
}

to_plot <- read.csv(file = arguments$input_file, header = FALSE)

colnames(to_plot) <- c("K", "Error")

plot_ks <- ggplot(data = to_plot, aes(as.factor(K), Error, fill = "red")) +
  geom_point(color = "black", alpha = 0.8, shape = 21, size = 6) +
  geom_line() +
  labs(
    title = "Admixture CV Error",
    x = "K",
    y = "Error"
  ) +
  theme( 
    plot.title = element_text(size = 30, hjust = 0.5),
    strip.text = element_text(size = 17),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

plot_file_name_pdf <- file.path(
  arguments$plot_folder,
  paste0(c(arguments$name, "_cv_errors.pdf"), collapse = '')
)

Cairo(file = plot_file_name_pdf, type = "pdf", width = 1024, height = 1024, dpi = 80)
print(plot_ks)
dev.off()
