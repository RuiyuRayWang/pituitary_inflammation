#' Set the default identity of cells
#' 
#' Used to reorder factor level order in CellChat objects.
#' 
#' 
#' @param object CellChat object
#' @param ident.use the name of the variable in object.meta;
#' @param levels set the levels of factor
#' @references Suoqin Jin. https://github.com/sqjin/CellChat/issues/149
#' @return
#' @export
#'
#' @examples
setIdent <- function(object, ident.use = NULL, levels = NULL){
  object@idents <- as.factor(object@meta[[ident.use]])
  if (!is.null(levels)) {
    object@idents <- factor(object@idents, levels = levels)
  }
  if (length(object@net) > 0) {
    if (all(dimnames(object@net$prob)[[1]] %in% levels(object@idents) )) {
      message("Reorder cell groups! ")
      idx <- match(dimnames(object@net$prob)[[1]], levels(object@idents))
      object@net$prob <- object@net$prob[idx, idx, ]
      object@net$pval <- object@net$pval[idx, idx, ]
      cat("The cell group order after reordering is ", dimnames(object@net$prob)[[1]],'\n')
    } else {
      message("Rename cell groups but do not change the order! ")
      cat("The cell group order before renaming is ", dimnames(object@net$prob)[[1]],'\n')
      dimnames(object@net$prob) <- list(levels(object@idents), levels(object@idents), dimnames(object@net$prob)[[3]])
      dimnames(object@net$pval) <- dimnames(object@net$prob)
      cat("The cell group order after renaming is ", dimnames(object@net$prob)[[1]],'\n')
    }
    
  }
  return(object)
}