# * Description
# This script takes the output of admixture and returns the plots

# * Libraries

list_of_packages <- c(
  "ggplot2", "reshape2",
  "forcats", "ggthemes",
  "patchwork", "gridExtra",
  "grid", "ggh4x",
  "stringr", "ggforce",
  "argparse", "stringr",
  "Cairo"
)

for (i in list_of_packages) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i)
  }
  library( i, character.only = T )
}

# * Reading Command Line

parser <- ArgumentParser()

parser$add_argument("-input_folder",
  default = T, action = "store", type = "character",
  help = "Folder of Q files"
)

parser$add_argument("-prefix",
  default = T, action = "store", type = "character",
  help = "Prefix of Q files"
)

parser$add_argument("-plot_folder",
  default = T, action = "store", type = "character",
  help = "Folder of output plots"
)

parser$add_argument("-label_file",
  default = T, action = "store", type = "character",
  help = "File that contains the Individuals' labels that are to be plotted in csv format"
)

parser$add_argument("-meta",
  default = T, action = "store", type = "character",
  help = "File that contains Individuals' meta information (.fam file)"
)

arguments <- parser$parse_args()

# * Loop that iterates over all Q files

meta_file <- arguments$meta
input_folder <- arguments$input_folder
plot_folder <- arguments$plot_folder
label_file <- arguments$label_file
prefix <- arguments$prefix

print(paste0(c("Meta File: ", meta_file), collapse = ''))
print(paste0(c("Input Folder: ", input_folder), collapse = ''))
print(paste0(c("Plot Folder: ", plot_folder), collapse = ''))
print(paste0(c("Labels: ", label_file), collapse = ''))

if (!dir.exists(plot_folder)) {
  dir.create(plot_folder, recursive = T)
  print("Created plot folder and parents")
}

q.files <- grep("\\.Q", list.files(input_folder), value = TRUE)
files_to_plot <- grep(prefix, q.files, value = TRUE)

for (files_to_plot in q.files) {
  
  temp_kappa_match <- regmatches(q.file.name, regexpr(pattern = "[0-9]+\\.Q", q.file.name))
  kappa <- as.numeric(strsplit(temp_kappa_match, "\\.")[[1]][1])
  ## prefix <- regmatches(q.file.name, regexpr(pattern = "[0-9]+\\.Q", q.file.name), invert = TRUE)[[1]][1]

  admix_data <- read.table(paste0(c(input_folder, q.file.name), collapse = "/"), header = FALSE)
  colnames(admix_data) <- paste0(rep("comp_", kappa), 1:kappa)

  labels <- read.csv(label_file, header = TRUE)
  colnames(labels) <- c("ID", "Date", "Include_1", "Include_2", "Location", "Layer_1", "Layer_2", "Layer_3")
  meta_data <- read.table(meta_file) ## We need this to fix the order of samples in the labels file
  ## colnames(meta_data) <- c("id", "sex", "Population")

  ## The following lines text mine names individual names for the date estimates
  temp_str_1 <- str_extract( labels$Include_1, "\\d+\\w*(\\-|_)\\d+" )
  temp_letters <- str_extract( temp_str_1, "[A-Z]+" )
  old_estimate <- as.numeric(str_extract( temp_str_1, "^\\d+" ))
  young_estimate <- as.numeric(str_extract( temp_str_1, "\\d+$" ))
  temp_sign_left <- ifelse((old_estimate > young_estimate),-1,1 )
  temp_sign_right <- ifelse(is.na(temp_letters), -1,1 )
  date_midpoint <- (old_estimate*temp_sign_left + young_estimate*temp_sign_right)/2
  labels$midpoint <- date_midpoint
  labels$Layer_3 <- reorder(labels$Layer_3, labels$midpoint, FUN = mean)
  
  ## Labels file isn't in the same order as meta file.
  correct_order_of_labels <- unlist(lapply(meta_data[, 2], function(x) {
    which(labels[, 1] %in% x)
  }))

  ## Not all individuals might be represented in labels file.
  names_in_labels <- unlist(lapply(labels[, 1], function(x) {
    which(meta_data[, 2] %in% x)
  }))

  labels <- labels[sort(correct_order_of_labels, decreasing = FALSE), ]
  admix_data <- admix_data[names_in_labels, ]
  meta_data <- meta_data[names_in_labels, ]

  to_melt <- data.frame(
    ID = labels$ID,
    Date = labels$Date,
    Location = labels$Location,
    Layer_1 = labels$Layer_1,
    Layer_2 = labels$Layer_2,
    Layer_3 = labels$Layer_3,
    Midpoint = labels$midpoint
  )
  to_melt <- cbind(to_melt, admix_data)

  melted_data <- melt(to_melt, measure.vars = paste0(rep("comp_", kappa), 1:kappa))

  ## Write the function that will make the plot

  melted_to_plot <- melted_data %>%
    dplyr::group_by(Location) %>%
    dplyr::arrange(Midpoint)
  ## colnames(melted_to_plot)

  k2plot <-
    ggplot(melted_data, aes(value, ID, fill = as.factor(variable))) +
    geom_col(color = "white", width = 0.99, linewidth = 0.005) +
    facet_nested(
      Location + Layer_3 + Layer_2 ~ .,
      shrink = TRUE, solo_line = T,
      scales = "free", switch = "y", space = "free_y",
      independent = "x", drop = TRUE, render_empty = FALSE,
      strip = strip_nested(
        text_y = elem_list_text(angle = 0),
        by_layer_y = TRUE, size = "variable", clip = "off"
      )
    ) +
    ## theme_minimal() +
    labs(
      x = "Individuals",
      title = paste0(c("K=", kappa), collapse = ""),
      y = "Ancestry"
    ) +
    theme(
      title = element_text(size = 24),
      axis.title = element_blank(),
      panel.spacing.y = unit(0.05, "lines"),
      strip.text.y = element_text(
        size = 14,
        margin = margin( 5, 4, 5, 4, unit = "mm" )
      ),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      ## strip.text.y = element_text( angle = 180 ),
      ggh4x.facet.nestline = element_line(colour = "blue")
    ) +
    scale_y_discrete( expand = expansion( mult = c(0.01,0.01) ) ) +
    scale_x_continuous( expand = expansion( mult = c(0.01,0.01) ) ) +
    ## scale_x_discrete(label=function(x) abbreviate(x, minlength=3, strict = TRUE)) +
    scale_fill_gdocs(guide = "none")
  
  k2plot_pdf <-
    ggplot(melted_data, aes(value, ID, fill = as.factor(variable))) +
    geom_col(width = 0.99, color = "white", linewidth = 0.005) +
    facet_nested(
      Location + Layer_3 + Layer_2 ~ .,
      shrink = TRUE, solo_line = T,
      scales = "free", switch = "y", space = "free_y",
      independent = "x", drop = TRUE, render_empty = FALSE,
      strip = strip_nested(
        text_y = elem_list_text(angle = 0),
        by_layer_y = TRUE, size = "variable", clip = "off"
      )
    ) +
    ## theme_minimal() +
    labs(
      x = "Individuals",
      title = paste0(c("K=", kappa), collapse = ""),
      y = "Ancestry"
    ) +
    theme(
      title = element_text( size = 24 ),
      axis.title = element_blank(),
      panel.spacing.y = unit(0.4, "lines"),
      strip.text.y = element_text(
        size = 20,
        margin = margin( 5, 4, 5, 4, unit = "mm" )
      ),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      ## strip.text.y = element_text( angle = 180 ),
      ggh4x.facet.nestline = element_line(colour = "blue")
    ) +
    ## facetted_pos_scales( y = 1 ) +
    scale_y_discrete( expand = expansion( mult = c(0.01,0.01) ) ) +
    scale_x_continuous( expand = expansion( mult = c(0.01,0.01) ) ) +
    ## scale_x_discrete(label=function(x) abbreviate(x, minlength=3, strict = TRUE)) +
    scale_fill_gdocs(guide = "none")
  
  plot_file_name <- paste0(plot_folder, "/", prefix, kappa, ".png", collapse = "")
  print(plot_file_name)
  Cairo( file = plot_file_name, type = "png", height = 3072, width = 1240, dpi = 45, pointsize = 12 )
  print(k2plot)
  dev.off()

  plot_file_name_pdf <- paste0(plot_folder, "/", prefix, kappa, ".pdf", collapse = "")
  print(plot_file_name_pdf)
  Cairo(file = plot_file_name_pdf, type = "pdf", height = 3072, dpi = 45, pointsize = 12 )
  print(k2plot_pdf)
  dev.off()
  
}
