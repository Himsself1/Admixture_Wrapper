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
  "Cairo", "plotly",
  "htmlwidgets", "htmltools"
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

parser$add_argument("-projected_input_folder",
  default = '', action = "store", type = "character",
  help = "Folder of Q files of projected individuals"
)


## parser$add_argument("-prefix",
##   default = T, action = "store", type = "character",
##   help = "Prefix of Q files"
## )

parser$add_argument("-plot_folder",
  default = T, action = "store", type = "character",
  help = "Folder of output plots"
)

parser$add_argument("-label_file",
  default = T, action = "store", type = "character",
  help = "File that contains Individuals' meta information (.fam file)"
)

parser$add_argument("-projected_label_file",
  default = '', action = "store", type = "character",
  help = "File that contains Individuals' meta information (.fam file) for projected individuals."
)

parser$add_argument("-meta",
  default = T, action = "store", type = "character",
  help = "File that contains the Individuals' labels that are to be plotted in tsv format"
)

parser$add_argument("-name",
  default = T, action = "store", type = "character",
  help = "Name of the plots"
)


arguments <- parser$parse_args()

# * Loop that iterates over all Q files

meta_file <- arguments$meta
input_folder <- arguments$input_folder
projected_input_folder <- arguments$projected_input_folder
plot_folder <- arguments$plot_folder
label_file <- arguments$label_file
projected_label_file <- arguments$projected_label_file
## prefix <- arguments$prefix
name <- arguments$name

print(paste0(c("Meta File: ", meta_file), collapse = ''))
print(paste0(c("Input Folder: ", input_folder), collapse = ''))
print(paste0(c("Projected Input Folder: ", projected_input_folder), collapse = ''))
print(paste0(c("Plot Folder: ", plot_folder), collapse = ''))
print(paste0(c("Labels: ", label_file), collapse = ''))
print(paste0(c("Projected Labels: ", projected_label_file), collapse = ''))

if (!dir.exists(plot_folder)) {
  dir.create(plot_folder, recursive = T)
  print("Created plot folder and parents")
}

q.files <- grep("Q$", list.files(input_folder, full.names = T), value = TRUE)
if( projected_input_folder > "" ){
  q.files.projected <- grep("Q$", list.files(projected_input_folder, full.names = T), value = TRUE)
}

for (counter in 1:length(q.files)) {
  files_to_plot <- q.files[counter]
  ## Find out which K is analyzed.
  temp_kappa_match <- regmatches(files_to_plot, regexpr(pattern = "[0-9]+\\.Q$", files_to_plot))
  kappa <- as.numeric(strsplit(temp_kappa_match, "\\.")[[1]][1])
  ## prefix <- regmatches(q.file.name, regexpr(pattern = "[0-9]+\\.Q", q.file.name), invert = TRUE)[[1]][1]
  
  admix_data <- read.table(files_to_plot, header = FALSE)
  colnames(admix_data) <- paste0(rep("comp_", kappa), 1:kappa)
  labels <- read.table(label_file, sep = "\t", header = F)
  meta_data <- read.table(meta_file, header = T, sep = '\t') ## We need this to fix the order of samples in the labels file
  ## Matches order of individuals of .ind/.fam and metadata
  if( grep( "\\.fam$", label_file )){
    correct_order_of_labels <- unlist(lapply(labels[, 2], function(x) {
      which(meta_data[, 1] %in% x)
    }))
    names_in_labels <- unlist(lapply(meta_data[, 1], function(x) {
      which(labels[, 2] %in% x)
    }))
  } else if(grep( "\\.ind$", label_file )){
    correct_order_of_labels <- unlist(lapply(labels[, 1], function(x) {
      which(meta_data[, 1] %in% x)
    }))
    names_in_labels <- unlist(lapply(meta_data[, 1], function(x) {
      which(labels[, 1] %in% x)
    }))
  } else {
    sprintf( "Label file needs to be .fam or .ind" )
    stop("Unknown label file format")
  }

  if( projected_input_folder > "" ){
    files_to_plot_projected <- q.files.projected[counter]

    projected_admix_data <- read.table(files_to_plot_projected, header = FALSE)
    colnames(projected_admix_data) <- paste0(rep("comp_", kappa), 1:kappa)
    projected_labels <- read.table(projected_label_file, sep = "\t", header = F)
    projected_meta_data <- read.table(meta_file, header = T, sep = '\t') ## We need this to fix the order of samples in the labels file
  ## Matches order of individuals of .ind/.fam and metadata
    if( grep( "\\.fam$", projected_label_file )){
      projected_correct_order_of_labels <- unlist(lapply(projected_labels[, 2], function(x) {
        which(projected_meta_data[, 1] %in% x)
      }))
      projected_names_in_labels <- unlist(lapply(projected_meta_data[, 1], function(x) {
        which(projected_labels[, 2] %in% x)
      }))
    } else if(grep( "\\.ind$", projected_label_file )){
      projected_correct_order_of_labels <- unlist(lapply(projected_labels[, 1], function(x) {
        which(projected_meta_data[, 1] %in% x)
      }))
      projected_names_in_labels <- unlist(lapply(projected_meta_data[, 1], function(x) {
        which(projected_labels[, 1] %in% x)
      }))
    } else {
      sprintf( "Label file needs to be .fam or .ind" )
      stop("Unknown label file format")
    }
    
    projected_meta_data <- projected_meta_data[sort(projected_correct_order_of_labels, decreasing = FALSE), ]
    projected_meta_data <- cbind( projected_meta_data, rep("Projected", nrow(projected_meta_data)) )
    colnames(projected_meta_data) <- c("ID", "Layer_1", "Layer_2", "Projection")
    projected_admix_data <- admix_data[projected_names_in_labels, ]
    projected_labels <- labels[projected_names_in_labels, ]

    projected_to_melt <- cbind(projected_meta_data, projected_admix_data)
    projected_melted_data <- melt(projected_to_melt, measure.vars = paste0(rep("comp_", kappa), 1:kappa))
  }
  
  ## Not all individuals might be represented in labels file.

  meta_data <- meta_data[sort(correct_order_of_labels, decreasing = FALSE), ]
  meta_data <- cbind( meta_data, rep("Unprojected", nrow(meta_data)) )
  colnames(meta_data) <- c("ID", "Layer_1", "Layer_2", "Projection")
  admix_data <- admix_data[names_in_labels, ]
  labels <- labels[names_in_labels, ]

  to_melt <- cbind(meta_data, admix_data)

  melted_data <- melt(to_melt, measure.vars = paste0(rep("comp_", kappa), 1:kappa))

  ## Write the function that will make the plot

  ## melted_to_plot <- melted_data %>%
  ##   dplyr::group_by(Layer_1) ## %>%
  ## dplyr::arrange(desc(comp_1), id)
  ## colnames(melted_to_plot)

# *** Plot in png
  
  k2plot <- ggplot(melted_data, aes(value, ID, fill = as.factor(variable))) +
    geom_col(color = "white", width = 0.99, linewidth = 0.005) +
    facet_nested(
      Layer_1 + Layer_2 ~ .,
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
      ## axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      ## strip.text.y = element_text( angle = 180 ),
      ggh4x.facet.nestline = element_line(colour = "blue")
    ) +
    scale_y_discrete( expand = expansion( mult = c(0.01,0.01) ), position = "right" ) +
    scale_x_continuous( expand = expansion( mult = c(0.01,0.01) ) ) +
    ## scale_x_discrete(label=function(x) abbreviate(x, minlength=3, strict = TRUE)) +
    scale_fill_gdocs(guide = "none")

# *** Plot in pdf
  
  k2plot_pdf <-
    ggplot(melted_data, aes(value, ID, fill = as.factor(variable))) +
    geom_col(width = 0.99, color = "white", linewidth = 0.005) +
    facet_nested(
      Layer_1 + Layer_2 ~ .,
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
      axis.text.y = element_text(size = 20),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      ## strip.text.y = element_text( angle = 180 ),
      ggh4x.facet.nestline = element_line(colour = "blue")
    ) +
    ## facetted_pos_scales( y = 1 ) +
    scale_y_discrete( expand = expansion( mult = c(0.01,0.01) ), position = "right" ) +
    scale_x_continuous( expand = expansion( mult = c(0.01,0.01) ) ) +
    ## scale_x_discrete(label=function(x) abbreviate(x, minlength=3, strict = TRUE)) +
    scale_fill_gdocs(guide = "none")
  
  plot_file_name <- file.path(plot_folder, paste0(c(name, kappa, ".png"), collapse = "") )
  print(plot_file_name)
  Cairo( file = plot_file_name, type = "png", height = 1440, width = 1024, dpi = 45, pointsize = 12 )
  print(k2plot)
  dev.off()

  plot_file_name_pdf <- file.path(plot_folder, paste0(c(name, kappa, ".pdf"), collapse = "") )
  print(plot_file_name_pdf)
  Cairo(file = plot_file_name_pdf, type = "pdf", height = 1440, dpi = 45, pointsize = 12 )
  print(k2plot_pdf)
  dev.off()

# ** Plots with projection

  if( projected_input_folder > "" ){
    melted_data <- rbind(melted_data, projected_melted_data)
  
# *** Plot in png
  
    projected_k2plot <- ggplot(melted_data, aes(value, ID, fill = as.factor(variable))) +
      geom_col(color = "white", width = 0.99, linewidth = 0.005) +
      facet_nested(
        Projection + Layer_1 + Layer_2 ~ .,
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
        ## axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        ## strip.text.y = element_text( angle = 180 ),
        ggh4x.facet.nestline = element_line(colour = "blue")
      ) +
      scale_y_discrete( expand = expansion( mult = c(0.01,0.01) ), position = "right" ) +
      scale_x_continuous( expand = expansion( mult = c(0.01,0.01) ) ) +
      ## scale_x_discrete(label=function(x) abbreviate(x, minlength=3, strict = TRUE)) +
      scale_fill_gdocs(guide = "none")
    
# *** Plot in pdf
  
    projected_k2plot_pdf <-
      ggplot(melted_data, aes(value, ID, fill = as.factor(variable))) +
      geom_col(width = 0.99, color = "white", linewidth = 0.005) +
      facet_nested(
        Projection + Layer_1 + Layer_2 ~ .,
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
        axis.text.y = element_text(size = 20),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        ## strip.text.y = element_text( angle = 180 ),
        ggh4x.facet.nestline = element_line(colour = "blue")
      ) +
      ## facetted_pos_scales( y = 1 ) +
      scale_y_discrete( expand = expansion( mult = c(0.01,0.01) ), position = "right" ) +
      scale_x_continuous( expand = expansion( mult = c(0.01,0.01) ) ) +
      ## scale_x_discrete(label=function(x) abbreviate(x, minlength=3, strict = TRUE)) +
      scale_fill_gdocs(guide = "none")
  
    projected_plot_file_name <- file.path(plot_folder, paste0(c(name, kappa, ".projected.png"), collapse = "") )
    print(projected_plot_file_name)
    Cairo( file = projected_plot_file_name, type = "png", height = 1440, width = 1024, dpi = 45, pointsize = 12 )
    print(projected_k2plot)
    dev.off()
    
    projected_plot_file_name_pdf <- file.path(plot_folder, paste0(c(name, kappa, ".projected.pdf"), collapse = "") )
    print(projected_plot_file_name_pdf)
    Cairo(file = projected_plot_file_name_pdf, type = "pdf", height = 1440, dpi = 45, pointsize = 12 )
    print(projected_k2plot_pdf)
    dev.off()
    
    ## This code is supposed to print interactive html plots, but they look horrible
    ## plot_file_name_html <- file.path(plot_folder, paste0(c(name, kappa, ".html"), collapse = "") )
    ## print(plot_file_name_html)
    ## html_plot <- ggplotly(k2plot_pdf)
    ## ## saveWidget(k2plot, plot_file_name_html)
    ## saveWidget(html_plot, basename(plot_file_name_html))
  }

}
