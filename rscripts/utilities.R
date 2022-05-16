#' @title Calculate Percent Cells Expressing a Given Feature
#' 
#' @description 
#' 
#' 
#' 
#' @details 
#' 
#' Inspired by a script written by Ryan-Zhu (https://github.com/satijalab/seurat/issues/371#issuecomment-486384854)
#' 
#' @param object A Seurat object.
#' @param features Features to calculate. By default, calculate all features.
#' @param assays Assay used to calculate pct.expr.
#' 
#' @examples 
#' 
#' @return A Seurat object with
#' @export
#' 
PercentExpr <- function(
  object,
  features = NULL,
  assays = "RNA"
){
  features <- features %||% rownames(object)
  counts = object[[assays]]@counts
  ncells = ncol(counts)
  
  features_use <- features[features %in% row.names(counts)]
  features_miss <- features[!features %in% row.names(counts)]
  
  if(length(features_miss) > 0){
    warning(paste0("The following feature(s) are not found in the assay \'", assays,"\': ",paste(features_miss, collapse = ", "),".\nOmitting them in calculation of pct.expr."))
  }
  
  if(length(features_use) > 0){
    pct.expr <- rowSums(counts[features_use,]>0)/ncells
    pct.expr <- data.frame(pct.expr)
    object[[assays]][["pct.expr"]] <- pct.expr
  }
  return(object)
}

#' @title  Calculate Marker Specificity
#' @description 
#' 
#' @details 
#' 
#' @import cummeRbund
#' @import reshape2
#' @export
CalculateMarkerSpecificity <- function(
  object,
  features = NULL,
  idents = 'ident',
  assays = "RNA",
  # min.pct.exp = NULL,
  return.Seurat = TRUE
) {
  features <- features %||% rownames(object)
  features_use <- features[features %in% row.names(object)]
  features_miss <- features[!features %in% row.names(object)]
  
  if(length(features_miss) > 0){
    warning(paste0("The following feature(s) are not found in the assay \'", assays,"\': ",paste(features_miss, collapse = ", "),".\nOmitting them in calculation of specificity."))
  }
  
  # min.pct.exp <- min.pct.exp %||% min(prop.table(table(object[[idents]])))
  avg_expr <- AverageExpression(object = object, features = features_use, assays = assays, group.by = idents)
  avg_expr <- avg_expr[[1]]
  
  marker_specificity <- lapply(1:ncol(avg_expr), function(cell_type_i){
    perfect_specificity <- rep(0.0, ncol(avg_expr))
    perfect_specificity[cell_type_i] <- 1.0
    apply(avg_expr, 1, function(x) { 
      if (sum(x) > 0) 1 - cummeRbund::JSdistVec(cummeRbund::makeprobsvec(x), perfect_specificity)
      else 0
    })
  })
  marker_specificity <- t(do.call(rbind, marker_specificity))
  colnames(marker_specificity) <- paste0("spec.", colnames(avg_expr))
  marker_specificity <- as.data.frame(marker_specificity)
  if (return.Seurat){
    for (i in 1:ncol(marker_specificity)){
      object[[assays]][[colnames(marker_specificity)[i]]] <- subset(x = marker_specificity, select = colnames(marker_specificity)[i])
    }
    return(object)
  } else {
    return(marker_specificity)
  }
}

# features = c("Cartpt","Nptx2","Lgals9","Nrtn","Serpina3n","Lox","Insl6","Fndc4","Lbp","Fgg","Isg15","Vgf","Il1r2")

## DEPRECATED
# # Credits to Ryan-Zhu (https://github.com/satijalab/seurat/issues/371)
# # updated 1/31/2020 to accommodate V3.1
# # updated 2/4/2020 to output "NA" for genes not detected in certain subgroups
# PrctCellExprGene <- function(
#   object, 
#   genes, 
#   group.by = "all") {
#   if(group.by == "all"){
#     prct = unlist(lapply(genes,calc_helper, object=object))
#     result = data.frame(Markers = genes, Cell_proportion = prct)
#     return(result)
#   }
#   
#   else{
#     list = SplitObject(object, group.by)
#     factors = names(list)
#     
#     results = lapply(list, PrctCellExpringGene, genes=genes)
#     for(i in 1:length(factors)){
#       results[[i]]$Feature = factors[i]
#     }
#     combined = do.call("rbind", results)
#     return(combined)
#   }
# }
