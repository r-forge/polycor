# last modified 2021-12-07 by J. Fox

"print.polycor" <-
function(x, digits = max(3, getOption("digits") - 3), ...){
  if (x$type == "polychoric"){
    if (!all(is.na(x$var))){
      se <- sqrt(diag(x$var))
      se.rho <- se[1]
    } else {
      se <- NA
      se.rho <- NA
    }
    est <- if (x$ML) "ML est." else "2-step est."
    if (!is.na(se.rho)){
      cat("\nPolychoric Correlation, ", est, " = ", signif(x$rho, digits),
        " (", signif(se.rho, digits), ")", sep="")
      if (x$df > 0)
          cat("\nTest of bivariate normality: Chisquare = ", 
          signif(x$chisq, digits), ", df = ", x$df, ", p = ", 
          signif(pchisq(x$chisq, x$df, lower.tail=FALSE), digits), "\n", sep="")
      else cat("\n")
    } else {
      cat("\nPolychoric Correlation, ", est, " = ", signif(x$rho, digits), "\n\n")
    }
    r <- length(x$row.cuts)
    c <- length(x$col.cuts)
    if (r == 0) return(invisible(x))
    if (!all(is.na(se))){
      row.cuts.se <- se[2:(r+1)]
      col.cuts.se <- se[(r+2):(r+c+1)]
    } else {
      row.cuts.se <- rep(NA, r)
      col.cuts.se <- rep(NA, c)
    }
    rowThresh <- signif(cbind(x$row.cuts, row.cuts.se), digits)
    if (r > 1) cat("\n  Row Thresholds") 
    else cat("\n  Row Threshold") 
    rownames(rowThresh) <- if (r > 1) 1:r else " "
    if (all(is.na(rowThresh[, 2]))) print(rowThresh[, 1, drop=FALSE]) 
    else {
      colnames(rowThresh) <- c("Threshold", "Std.Err.")
      cat("\n")
      print(rowThresh)
    }
    colThresh <- signif(cbind(x$col.cuts, col.cuts.se), digits)
    if (c > 1) cat("\n\n  Column Thresholds")
    else cat("\n\n  Column Threshold")
    rownames(colThresh) <- if (c > 1) 1:c else " "
    if (all(is.na(colThresh[, 2]))) print(colThresh[, 1, drop=FALSE]) 
    else {
      colnames(colThresh) <- c("Threshold", "Std.Err.")
      cat("\n")
      print(colThresh)
      }
    }
  else if (x$type == "polyserial"){
    if (!all(is.na(x$var))) {
      se <- sqrt(diag(x$var))
      se.rho <- se[1]
    } else {
      se <- NA
      se.rho <- NA
    }
    est <- if (x$ML) "ML est." else "2-step est."
    if (!all(is.na(se))){
      cat("\nPolyserial Correlation, ", est, " = ", signif(x$rho, digits),
        " (", signif(se.rho, digits), ")", sep="")
      cat("\nTest of bivariate normality: Chisquare = ", signif(x$chisq, digits),
        ", df = ", x$df, ", p = ", signif(pchisq(x$chisq, x$df, lower.tail=FALSE), digits),
        "\n\n", sep="")
    } else {
      cat("\nPolyserial Correlation, ", est, " = ", signif(x$rho, digits), "\n\n")
    }
    if (length(se) > 1) cuts.se <- se[-1] else cuts.se <- rep(NA, length(x$cuts))
    thresh <- signif(rbind(x$cuts, cuts.se), digits)
    colnames(thresh) <- 1:length(x$cuts)
    rownames(thresh) <- c("Threshold", "Std.Err.")
    if (all(is.na(thresh[2, ]))) thresh <- thresh[-2, , drop=FALSE]
    print(thresh)
    }
  else print(unclass(x))
  invisible(x)
  }
