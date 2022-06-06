#' Set a default value if an object is null
#' 
#' @param lhs An object to set if it's null
#' @param rhs The value to provide if x is null
#' 
#' @return rhs if lhs is null, else lhs
#' 
#' @author Hadley Wickham
#' @references https://adv-r.hadley.nz/functions.html#missing-arguments
#' 
`%||%` <- function(lhs, rhs) {
  if(!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

#' Set a default value if an object is NOT null
#' 
#' @param lhs An object to set if it's NOT null
#' @param rhs The value to provide if x is NOT null
#' 
#' @return lhs if lhs is null, else rhs
#' 
#' @author Hadley Wickham
#' @references https://adv-r.hadley.nz/functions.html#missing-arguments
#' 
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)){
    return(rhs)
  } else {
    return(lhs)
  }
}

#' @title Calculate Percent Cells Expressing a Given Feature
#' 
#' @description 
#' 
#' @details 
#' 
#' @param object A Seurat object.
#' @param features Features to calculate. By default, calculate all features.
#' @param group.by Whether to calculate percent expression on group of cells. By default, calculate on all cells ("all").
#' @param assays Assay used to calculate pct.expr.
#' 
#' @examples 
#' 
#' @return A Seurat object with
#' 
#' @author Ryan-Zhu
#' @references https://github.com/satijalab/seurat/issues/371#issuecomment-486384854
#' 
#' @import SeuratObject
#' @export
#' 
PercentExpr <- function(
  object,
  features = NULL,
  assays = "RNA"
){
  features <- features %||% rownames(object)
  # counts = object[[assays]]@counts
  counts = GetAssayData(object = object, assay = assays, slot = "counts")
  ncells = ncol(object)
  
  features_use <- features[features %in% row.names(object)]
  features_dump <- features[!features %in% row.names(object)]
  
  if(length(features_dump) > 0){
    warning(paste0("The following feature(s) are not found in the assay \'", assays,"\': ",paste(features_dump, collapse = ", "),".\nOmitting them in calculation of pct.expr."))
  }
  
  pct.expr <- rowSums(counts[features_use,]>0)/ncells
  pct.expr <- data.frame(pct.expr)
  object[[assays]][["pct.expr"]] <- pct.expr

  return(object)
}

#' @title Calculate Scaling Factor
#' 
#' @import tibble
CalculateScalingFactor <- function(
  object,
  group.by = 'ident',
  assays = "RNA"
) {
  df <- FetchData(object, vars = c("nFeature_RNA", group.by))
  df <- df %>% group_by(cell_type_brief) %>% summarise(median = median(nFeature_RNA))
  df <- mutate(df, scaling_fct = mean(median)/median)
  df <- column_to_rownames(.data = df, var = group.by)
  return(df)
}

#' @title
#' Calculate Marker Specificity
#' 
#' @description 
#' 
#' @details 
#' 
#' @param object
#' @param features
#' @param group.by
#' @param profile_use
#' @param assays
#'
#' @return A Seurat object with marker specificity scores attached to feature metadata.
#' @references Qiu, X. et al. Reversed graph embedding resolves complex single-cell trajectories. Nature Methods 14, 979-982 (2017).
#' 
#' @import dplyr
#' @import cummeRbund
#' @import reshape2
#' @import SeuratObject
#' @export
CalculateMarkerSpecificity <- function(
  object,
  features = NULL,
  group.by = 'ident',
  profile_use = "avgexpr",
  assays = "RNA"
) {
  features <- features %||% rownames(object)
  features_use <- features[features %in% row.names(object)]
  features_miss <- features[!features %in% row.names(object)]
  
  if(length(features_miss) > 0){
    warning(paste0("The following feature(s) are not found in the assay \'", assays,"\': ",paste(features_miss, collapse = ", "),".\nOmitting them in calculation of specificity."))
  }
  
  obj_sub <- subset(object, features = features_use)
  
  if(profile_use == "avgexpr"){
    feature_profile <- AverageExpression(object = object, features = features_use, assays = assays, group.by = group.by)
    feature_profile <- feature_profile[[assays]]
  } else if (profile_use == "proportion") {
    object.list <- SplitObject(object = obj_sub, split.by = group.by)
    object.list <- lapply(X = object.list, FUN = PercentExpr)
    
    df_scale_fct <- CalculateScalingFactor(object = object, group.by = group.by, assays = assays)
    
    feature_profile <- data.frame(row.names = features_use)
    for (i in seq_along(object.list)){
      cell_type = names(object.list[i])
      # object.list[[i]][[assays]][[paste0("norm.pct.expr", cell_type)]] <- object.list[[i]][[assays]][["pct.expr"]] * df_scale_fct[cell_type,"scaling_fct"]
      # object[[assays]][[paste0("pct.expr.", cell_type)]] <- object.list[[i]][[assays]][["pct.expr"]]
      # object[[assays]][[paste0("norm.pct.expr", cell_type)]] <- object.list[[i]][[assays]][["pct.expr"]] * df_scale_fct[cell_type,"scaling_fct"]
      feature_profile[[cell_type]] <- object.list[[i]][[assays]][["pct.expr"]] * df_scale_fct[cell_type,"scaling_fct"]
    }
  } else {
    stop("wrong argument: 'profile_use' should be one of 'avgexpr' or 'proportion'.")
  }
  
  marker_specificity <- lapply(1:ncol(feature_profile), function(cell_type_i){
    perfect_specificity <- rep(0.0, ncol(feature_profile))
    perfect_specificity[cell_type_i] <- 1.0
    apply(feature_profile, 1, function(x) { 
      if (sum(x) > 0) 1 - cummeRbund::JSdistVec(cummeRbund::makeprobsvec(x), perfect_specificity)
      else 0
    })
  })
  marker_specificity <- t(do.call(rbind, marker_specificity))
  colnames(marker_specificity) <- paste0("spec.", colnames(feature_profile))
  marker_specificity <- as.data.frame(marker_specificity)
  for (i in 1:ncol(marker_specificity)){
    object[[assays]][[colnames(marker_specificity)[i]]] <- dplyr::subset(x = marker_specificity, select = colnames(marker_specificity)[i])
  }
  
  return(object)
}