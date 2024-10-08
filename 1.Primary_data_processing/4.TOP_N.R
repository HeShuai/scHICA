options(width = 240)

TOP_N <- function(x, n, pct.1 = 0.1, first = "avg_log2FC", second = "p_val_adj", sig.padj = NULL, fc.threshold = c(-0.25, 0.25)){
  if(table(names(x) %in% c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))["TRUE"] == 7){
    reordered_daf <- data.frame()
    if (!is.null(pct.1) & (pct.1 > 0 & pct.1 < 1)) {
      x <- x[x$pct.1 >= pct.1, ]
    }
    if (!is.null(sig.padj) & is.numeric(sig.padj)) {
      x <- x[x$p_val_adj <= sig.padj, ]
    }
    if (length(fc.threshold) == 2) { 
      fc.threshold <- sort(fc.threshold)
      x <- x[(x$avg_log2FC < fc.threshold[1] | x$avg_log2FC > fc.threshold[2]), ]
    } else if (length(fc.threshold) == 1) {
      x <- x[x$avg_log2FC >= fc.threshold, ]
    }
    
    for (i in unique(x$cluster)) {
      if(n > dim(x[x$cluster == i, ])[1]){
        message("FBI warning, n < ", n)
        tmp <- dplyr::arrange(.data = x[x$cluster == i, ], desc(get(first)), desc(get(second)))[1:(dim(x[x$cluster == i, ])[1]), ]
        reordered_daf <- rbind(reordered_daf, tmp)
      } else {
        tmp <- dplyr::arrange(.data = x[x$cluster == i, ], desc(get(first)), desc(get(second)))[1:n, ]
        reordered_daf <- rbind(reordered_daf, tmp)
      }
    }
  }
  else {
    message("be careful !")
  }
  return(reordered_daf)
}


