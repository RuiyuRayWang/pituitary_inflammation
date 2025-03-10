plot_umap_refactored <- function (augur, sc, mode = c("default", "rank"), reduction = "umap", direction = 1,
          palette = "cividis", top_n = 0, limits = NULL, cell_type_col = "cell_type", size_sm = 6, size_lg = 7, lgd.pos = "right", lgd.name = "AUC")
{
  require(tidyverse)
  require(ggrepel)
  require(magrittr)
  require(viridis)
  
  mode = match.arg(mode)
  aucs = augur$AUC
  if (mode == "rank") {
    aucs %<>% mutate(rank = rank(auc), rank_pct = rank/n(), 
                     rank_pct = scales::rescale(rank_pct, c(0, 1)), fill = rank_pct)
    legend_name = "Rank (%)"
    label_fun = function(x) x * 100
    breaks = c(0, 1)
    color_labels = c(0, 100)
  }
  else {
    aucs %<>% mutate(fill = auc)
    legend_name = lgd.name
    if (!is.null(limits)) {
      breaks = limits
      color_labels = breaks
    }
    else {
      breaks = range(aucs$auc)
      color_labels = format(breaks, format = "f", digits = 2)
    }
  }
  if ("Seurat" %in% class(sc)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("install \"Seurat\" R package for Augur compatibility with ", 
           "input Seurat object", call. = FALSE)
    }
    meta = sc@meta.data %>% as.data.frame()
    red_coord = sc@reductions[[reduction]]@cell.embeddings
  }
  else if ("monocle3" %in% class(sc)) {
    if (!requireNamespace("monocle3", quietly = TRUE)) {
      stop("install \"monocle3\" R package for Augur compatibility with ", 
           "input monocle3 object", call. = FALSE)
    }
    meta = monocle3::pData(sc) %>% droplevels() %>% as.data.frame()
    reduction = toupper(reduction)
    red_coord = sc@int_colData@listData$reducedDims[[reduction]]
  }
  else if ("SingleCellExperiment" %in% class(sc)) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("install \"SingleCellExperiment\" R package for Augur ", 
           "compatibility with input SingleCellExperiment object", 
           call. = FALSE)
    }
    meta = SummarizedExperiment::colData(sc) %>% droplevels() %>% 
      as.data.frame()
    reduction = toupper(reduction)
    red_coord = sc@int_colData@listData$reducedDims[[reduction]]
    rownames(red_coord) = colnames(sc)
  }
  colnames(red_coord)[1:2] = c("coord_x", "coord_y")
  meta$auc = aucs$fill[match(meta[[cell_type_col]], aucs$cell_type)]
  labeled_types = aucs %>% arrange(desc(fill)) %>% head(top_n) %>% 
    pull(cell_type)
  plot_data = red_coord %>% as.data.frame() %>% rownames_to_column(var = "barcode") %>% 
    left_join(meta %>% rownames_to_column(var = "barcode") %>% 
                dplyr::select(barcode, all_of(cell_type_col), auc), 
              by = "barcode")
  colnames(plot_data)[colnames(plot_data) == cell_type_col] = "cell_type"
  labels = plot_data %>% filter(cell_type %in% labeled_types) %>% 
    group_by(cell_type) %>% summarise(coord_x = median(coord_x), coord_y = median(coord_y), auc = mean(auc)) %>% drop_na()
  # xlab = paste(ifelse(reduction == "umap", "UMAP", reduction), 1)
  # ylab = paste(ifelse(reduction == "umap", "UMAP", reduction), 2)
  ## Only tested for Seurat objects
  xlab = colnames(Embeddings(sc, reduction = "umap"))[1]
  ylab = colnames(Embeddings(sc, reduction = "umap"))[2]
  p = ggplot(plot_data, aes(x = coord_x, y = coord_y, color = auc, 
                            fill = cell_type)) + 
    geom_point(size = 0.4, stroke = 0, shape = 16) + labs(x = xlab, y = ylab) + 
    guides(fill = "none", color = guide_colorbar(nbin = 10, raster = FALSE, ticks = FALSE, title.position = "top", title.hjust = 0.5)) + 
    geom_text_repel(data = labels, aes(label = cell_type), color = "black", size = 2, segment.size = 0, box.padding = 0.5, min.segment.length = 0.33) + 
    guides(color = guide_colorbar(frame.colour = "black", ticks = FALSE)) + theme_bw() + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank(), axis.title.x = element_text(size = size_lg), 
          axis.title.y = element_text(size = size_lg), panel.grid = element_blank(), 
          strip.text = element_text(size = size_lg), strip.background = element_blank(), 
          axis.line.y = element_blank(), axis.line.x = element_blank(), 
          legend.position = lgd.pos, legend.justification = "bottom", 
          legend.text = element_text(size = size_sm), legend.title = element_text(size = size_sm, hjust = 0.5), 
          legend.key.size = unit(0.5, "lines"), legend.margin = margin(rep(0, 4)), legend.background = element_blank(), 
          plot.title = element_text(size = size_lg, hjust = 0.5), aspect.ratio = 1)
  p
  if (length(palette) == 1 && palette %in% c("viridis", "cividis", "plasma", "magma", "inferno")) {
    p = p + scale_color_viridis(option = palette, name = legend_name, 
                                labels = color_labels, limits = breaks, breaks = breaks, 
                                na.value = "white")
  }
  else if (length(palette) == 1 && palette %in% c("BrBG", "PiYG", 
                                                  "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", 
                                                  "Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", 
                                                  "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", 
                                                  "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")) {
    p = p + scale_color_distiller(palette = palette, name = legend_name, 
                                  labels = color_labels, limits = breaks, breaks = breaks, direction = direction,
                                  na.value = "white")
  }
  else {
    p = p + scale_color_gradientn(colours = palette, name = legend_name, 
                                  labels = color_labels, limits = breaks, breaks = breaks, 
                                  na.value = "white")
  }
  return(p)
}