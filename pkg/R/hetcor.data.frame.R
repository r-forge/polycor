# last modified 2022-01-10 by J. Fox

# the following function to be imported from admisc and then deleted here:

# `tryCatchWEM` <- function(expr, capture = FALSE) {
#     toreturn <- list()
#     
#     output <- withVisible(withCallingHandlers(
#         tryCatch(expr, error = function(e) {
#             toreturn$error <<- e$message
#             NULL
#         }),
#         warning = function(w) {
#             toreturn$warning <<- c(toreturn$warning, w$message)
#             invokeRestart("muffleWarning")
#         },
#         message = function(m) {
#             toreturn$message <<- paste(toreturn$message, m$message, sep = "")
#             invokeRestart("muffleMessage")
#         }
#     ))
#     
#     if (capture && output$visible && !is.null(output$value)) {
#         toreturn$output <- capture.output(output$value)
#         toreturn$value <- output$value
#     }
#     
#     if (length(toreturn) > 0) {
#         return(toreturn)
#     }
# }


hetcor.data.frame <- function(data, ML=FALSE, std.err=TRUE, use=c("complete.obs", "pairwise.complete.obs"),
                              bins=4, pd=TRUE, parallel=FALSE, ncores=detectCores(logical=FALSE), 
                              thresholds=FALSE, ...){
    
    se.r <- function(r, n){
        rho <- r*(1 + (1 - r^2)/(2*(n - 3))) # approx. unbiased estimator
        v <- (((1 - rho^2)^2)/(n + 6))*(1 + (14 + 11*rho^2)/(2*(n + 6)))
        sqrt(v)
    }
    
    computeCor <- function(pair){
        type <- ""
        se <- NA
        test <- NA
        i <- rows[pair]
        j <- cols[pair]
        x <- data[, i]
        y <- data[, j]
        n <- sum(complete.cases(x, y))
        if (n == 0) {
            test <- r <- se <- NA
            warning("no cases for pair ", j, ", ", i)
        }
        if (inherits(x, c("numeric", "integer")) && inherits(y, c("numeric", "integer"))) {
            r <- cor(x, y, use="complete.obs")
            type <- "Pearson"
            if (std.err) {
                se <- se.r(r, n)
                test <- pchisq(chisq(x, y, r, bins=bins), bins^2 - 2, lower.tail=FALSE)
            }
            Thresholds <- list(NULL)
        }
        else if (inherits(x, c("factor", "logical", "character")) && 
                 inherits(y, c("factor", "logical", "character"))) {
            type <- "Polychoric"
            result <- tryCatchWEM(polychor(x, y, ML=ML, std.err=std.err, thresholds=thresholds),
                                  capture=TRUE)
            error <- !is.null(result$error)
            if (!is.null(result$warning)){
                warning("polychoric correlation between variables ", vnames[j], " and ", vnames[i],
                        if (length(result$warning) == 1) " produced a warning:\n" else " produced warnings:\n",
                        paste(paste("  ", result$warning), collapse="\n"))
            }
            if (error){
                msg <- result$error
                warning("could not compute polychoric correlation between variables ", vnames[j], " and ", vnames[i], "\n",
                        "   Error message: ", msg)
                result <- NA
            }
            if (std.err && !error){
                result <- result$value
                if (!(length(result) == 1 && is.na(result))){
                    r <- result$rho
                    se <- sqrt(result$var[1,1])
                    test <- if (result$df > 0)
                        pchisq(result$chisq, result$df, lower.tail=FALSE)
                    else NA
                }
                else {
                    r <- if (is.list(result)) result$value else result
                    test <- se <- NA
                }
            }
            else {
                r <- if (is.list(result)) result$value else result
                test <- se <- NA
            }
            Thresholds <- if (thresholds) {
                list(row.cuts=as.vector(result$row.cuts), 
                     col.cuts=as.vector(result$col.cuts))
            } else {
                NULL
            }
        }
        else {
            if (inherits(x, c("factor", "logical", "character")) && 
                inherits(y, c("numeric", "integer")))
                result <- tryCatchWEM(polyserial(y, x, ML=ML, std.err=std.err, bins=bins, 
                                                 thresholds=thresholds),
                                      capture=TRUE)
            else if (inherits(x, c("numeric", "integer")) && 
                     inherits(y, c("factor", "logical", "character")))
                result <- tryCatchWEM(polyserial(x, y, ML=ML, std.err=std.err, bins=bins),
                                      capture=TRUE)
            else {
                stop("columns must be numeric, factors, logical, or character.")
            }
            type <- "Polyserial"
            error <- !(is.null(result$error))
            if (!is.null(result$warning)){
                warning("polyserial correlation between variables ", vnames[j], " and ", vnames[i],
                        if (length(result$warning) == 1) " produced a warning:\n" else " produced warnings:\n",
                        paste(paste( "  ", result$warning), collapse="\n"))
            }
            if (error){
                msg <- result$error
                warning("could not compute polyserial correlation between variables ", vnames[j], " and ", vnames[i], "\n",
                        "   Error message: ", msg)
                result <- NA
            }
            if (std.err && !error){
                result <- result$value
                if (!(length(result) == 1 && is.na(result))){
                    r <- result$rho
                    se <- sqrt(result$var[1,1])
                    test <- pchisq(result$chisq, result$df, lower.tail=FALSE)
                }
                else {
                    r <- if (is.list(result)) result$value else result
                    test <- se <- NA
                }
            }
            else {
                r <- if (is.list(result)) result$value else result
                se <- test <- NA
            }
            Thresholds <- if (thresholds) {
                list(cuts=as.vector(result$cuts))
            } else {
                NULL
            }
        }
        list(n=n, r=r, Type=type, SE=se, Test=test, Thresholds=Thresholds)
    }
    
    vnames <- names(data)
    if (any(sapply(data, function(x) inherits(x, "character")))){
        message("data contain one or more character variables",
                "\nthe values of which are ordered alphabetically")
    }
    use <- match.arg(use)
    if (use == "complete.obs") {
        data <- na.omit(data)
        n <- nrow(data)
    }
    p <- length(data)
    if (p < 2) stop("fewer than 2 variables.")
    R <- matrix(1, p, p)
    Type <- matrix("", p, p)
    SE <- matrix(0, p, p)
    N <- matrix(0, p, p)
    Test <- matrix(0, p, p)
    if (thresholds){
        Thresholds <- vector(p^2, mode="list")
        Thresholds <- matrix(Thresholds, p, p)
    }
    diag(N) <- if (use == "complete.obs") nrow(data)
    else sapply(data, function(x) sum(!is.na(x)))
    if (all(diag(N) == 0)) stop("no non-missing cases")
    npairs <- p*(p -1)/2
    rows <- matrix(1:p, p, p)
    cols <- t(rows)
    rows <- rows[lower.tri(rows)]
    cols <- cols[lower.tri(cols)]
    result <- if (parallel && ncores > 1){
        message("Note: using a cluster of ", ncores, " cores")
        cl <- parallel::makeCluster(ncores)
        on.exit(parallel::stopCluster(cl))
        parallel::clusterApply(cl, 1:npairs, computeCor)
    } else {
        lapply(1:npairs, computeCor)
    }
    for (pair in 1:npairs){
        i <- rows[pair]
        j <- cols[pair]
        res <- result[[pair]]
        N[i, j] <- N[j, i] <- res$n
        R[i, j] <- R[j, i] <- res$r
        Type[i, j] <- Type[j, i] <- res$Type
        SE[i, j] <- SE[j, i] <- res$SE
        Test[i, j] <- Test[j, i] <- res$Test
        if (thresholds) {
            Thresholds[[i, j]] <- res$Thresholds
            Thresholds[[j, i]] <- res$Type
        }
    }
    if (pd && !any(is.na(R)) && min(eigen(R, only.values=TRUE)$values) < 0){
        cor <- Matrix::nearPD(R, corr=TRUE)
        if (!cor$converged) warning("attempt to make correlation matrix positive-definite failed")
        else warning("the correlation matrix has been adjusted to make it positive-definite")
        R <- as.matrix(cor$mat)
    }
    rownames(R) <- colnames(R) <- names(data)
    result <- list(correlations=R, type=Type, NA.method=use, ML=ML)
    if (thresholds) result$thresholds <- Thresholds
    if (std.err) {
        rownames(SE) <- colnames(SE) <- names(data)
        rownames(N) <- colnames(N) <- names(N)
        rownames(Test) <- colnames(Test) <- names(data)
        result$std.errors <- SE
        result$n <- if (use == "complete.obs") n else N
        result$tests <- Test
    }
    if (0 < (nNA <- sum(is.na(R[lower.tri(R)])))){
        warning(nNA, if (nNA == 1) " correlation" else " correlations", 
                " couldn't be computed and", if (nNA == 1) " is" else " are", " NA")
    }
    class(result) <- "hetcor"
    result
}
