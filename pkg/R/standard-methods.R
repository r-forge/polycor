coef.polycor <- function(object, correlation=TRUE, thresholds=TRUE, ...){
  result <- if (correlation) c(rho=object$rho) else numeric(0)
  if (thresholds){
    if (object$type == "polychoric"){
      row.cuts <- object$row.cuts
      if (!all(is.na(row.cuts))){
        names(row.cuts) <- paste0("row.threshold.", seq_along(row.cuts))
        result <- c(result, row.cuts)
      }
      col.cuts <- object$col.cuts
      if (!all(is.na(col.cuts))){
        names(col.cuts) <- paste0("col.threshold.", seq_along(col.cuts))
        result <- c(result, col.cuts)
      }
    } else {
      cuts <- object$cuts
      if (!all(is.na(cuts))){
        names(cuts) <- paste0("threshold.", seq_along(cuts))
        result <- c(result, cuts)
      }
    }
  }
  result
}

vcov.polycor <- function(object, correlation=TRUE, thresholds=TRUE, ...){
  if (!correlation && !thresholds) return(NULL)
  vc <- object$var
  if (is.null(vc) ||  all(is.na(vc))) return(NA)
  if (length(vc) > 1){
    rownames(vc) <- colnames(vc) <- names(coef(object))
  }
  if (correlation && (!thresholds)) {
    if (length(vc) == 1) {
      return(vc)
    } else {
      return(vc[1, 1])
    }
  }
  if ((!correlation) && thresholds){
    if (length(vc) == 1) {
      return(NA) 
    } else {
      return(vc[-1, -1])
    }
  }
  return(vc)
}

summary.polycor <- function(object, ...){
  print(object, ...)
}
