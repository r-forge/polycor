# # last modified 2021-11-30 by J. Fox

hetcor.data.frame <- function(data, ML=FALSE, std.err=TRUE, use=c("complete.obs", "pairwise.complete.obs"),
             bins=4, pd=TRUE, parallel=FALSE, ncores=detectCores(logical=FALSE), ...){

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
            }
            else if (inherits(x, c("factor", "logical", "character")) && 
                     inherits(y, c("factor", "logical", "character"))) {
                type <- "Polychoric"
                result <- try(polychor(x, y, ML=ML, std.err=std.err), silent=TRUE)
                error <- inherits(result, "try-error")
                if (error){
                    warning("could not compute polychoric correlation between variables ", i, " and ", j,
                            "\n    Message: ", result, "\n")
                    result <- NA
                }
                if (std.err && !error && !(length(result) == 1 && is.na(result))){
                    r <- result$rho
                    se <- sqrt(result$var[1,1])
                    test <- if (result$df > 0)
                        pchisq(result$chisq, result$df, lower.tail=FALSE)
                    else NA
                }
                else {
                    r <- result
                    test <- se <- NA
                }
            }
            else {
                if (inherits(x, c("factor", "logical", "character")) && 
                    inherits(y, c("numeric", "integer")))
                    result <- try(polyserial(y, x, ML=ML, std.err=std.err, bins=bins), silent=TRUE)
                else if (inherits(x, c("numeric", "integer")) && 
                         inherits(y, c("factor", "logical", "character")))
                    result <- try(polyserial(x, y, ML=ML, std.err=std.err, bins=bins), silent=TRUE)
                else {
                    stop("columns must be numeric, factors, logical, or character.")
                }
                type <- "Polyserial"
                error <- inherits(result, "try-error")
                if (error){
                    warning("could not compute polyserial correlation between variables ", i, " and ", j,
                            "\n    Message: ", result, "\n")
                    result <- NA
                }
                if (std.err && !error && !(length(result) == 1 && is.na(result))){
                    r <- result$rho
                    se <- sqrt(result$var[1,1])
                    test <- pchisq(result$chisq, result$df, lower.tail=FALSE)
                }
                else {
                    r <- result
                    se <- test <- NA
                }
            }
            list(n=n, r=r, Type=type, SE=se, Test=test)
        }
        
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
        }
        if (pd && !any(is.na(R)) && min(eigen(R, only.values=TRUE)$values) < 0){
            cor <- Matrix::nearPD(R, corr=TRUE)
            if (!cor$converged) warning("attempt to make correlation matrix positive-definite failed")
            else warning("the correlation matrix has been adjusted to make it positive-definite")
            R <- as.matrix(cor$mat)
        }
        rownames(R) <- colnames(R) <- names(data)
        result <- list(correlations=R, type=Type, NA.method=use, ML=ML)
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

# "hetcor.data.frame" <-
#     function(data, ML=FALSE, std.err=TRUE, use=c("complete.obs", "pairwise.complete.obs"),
#              bins=4, pd=TRUE, ...){
#         se.r <- function(r, n){
#             rho <- r*(1 + (1 - r^2)/(2*(n - 3))) # approx. unbiased estimator
#             v <- (((1 - rho^2)^2)/(n + 6))*(1 + (14 + 11*rho^2)/(2*(n + 6)))
#             sqrt(v)
#         }
#         if (any(sapply(data, function(x) inherits(x, "character")))){
#           message("data contain one or more character variables",
#                   "\nthe values of which are ordered alphabetically")
#         }
#         use <- match.arg(use)
#         if (use == "complete.obs") data <- na.omit(data)
#         p <- length(data)
#         if (p < 2) stop("fewer than 2 variables.")
#         R <- matrix(1, p, p)
#         Type <- matrix("", p, p)
#         SE <- matrix(0, p, p)
#         N <- matrix(0, p, p)
#         Test <- matrix(0, p, p)
#         diag(N) <- if (use == "complete.obs") nrow(data)
#         else sapply(data, function(x) sum(!is.na(x)))
#         if (all(diag(N) == 0)) stop("no non-missing cases")
#         for (i in 2:p) {
#             for (j in 1:(i-1)){
#                 x <- data[[i]]
#                 y <- data[[j]]
#                 n <- sum(complete.cases(x, y))
#                 if (n == 0) {
#                     Test[i, j] <- Test[j, i] <- R[i, j] <- R[j, i] <- SE[i, j] <- SE[j, i] <- NA
#                     N[i, j] <- N[j, i] <- 0
#                     warning("no cases for pair ", j, ", ", i)
#                     next
#                 }
#                 if (inherits(x, c("numeric", "integer")) && inherits(y, c("numeric", "integer"))) {
#                     r <- cor(x, y, use="complete.obs")
#                     Type[i, j] <- Type[j, i] <- "Pearson"
#                     R[i, j] <- R[j, i] <- r
#                     if (std.err) {
#                         SE[i, j] <- SE[j, i] <- se.r(r, n)
#                         N[i, j] <- N[j, i] <- n
#                         Test[i, j] <- pchisq(chisq(x, y, r, bins=bins), bins^2 - 2, lower.tail=FALSE)
#                     }
#                 }
#                 else if (inherits(x, c("factor", "logical", "character")) && 
#                          inherits(y, c("factor", "logical", "character"))) {
#                     Type[i, j] <- Type[j, i] <- "Polychoric"
#                     result <- try(polychor(x, y, ML=ML, std.err=std.err), silent=TRUE)
#                     error <- inherits(result, "try-error")
#                     if (error){
#                         warning("could not compute polychoric correlation between variables ", i, " and ", j,
#                                 "\n    Message: ", result, "\n")
#                         result <- NA
#                     }
#                     if (std.err && !error){
#                         R[i, j] <- R[j, i] <- result$rho
#                         SE[i, j] <- SE[j, i] <- sqrt(result$var[1,1])
#                         N[i, j] <- N[j, i] <- n
#                         Test[i, j] <- if (result$df > 0)
#                             pchisq(result$chisq, result$df, lower.tail=FALSE)
#                         else NA
#                     }
#                     else R[i, j] <- R[j, i] <- result
#                 }
#                 else {
#                     if (inherits(x, c("factor", "logical", "character")) && 
#                         inherits(y, c("numeric", "integer")))
#                         result <- try(polyserial(y, x, ML=ML, std.err=std.err, bins=bins), silent=TRUE)
#                     else if (inherits(x, c("numeric", "integer")) && 
#                              inherits(y, c("factor", "logical", "character")))
#                         result <- try(polyserial(x, y, ML=ML, std.err=std.err, bins=bins), silent=TRUE)
#                     else {
#                         stop("columns must be numeric, factors, logical, or character.")
#                     }
#                     Type[i, j] <- Type[j, i] <- "Polyserial"
#                     error <- inherits(result, "try-error")
#                     if (error){
#                         warning("could not compute polyserial correlation between variables ", i, " and ", j,
#                                 "\n    Message: ", result, "\n")
#                         result <- NA
#                     }
#                     if (std.err && !error){
#                         R[i, j] <- R[j, i] <- result$rho
#                         SE[i, j] <- SE[j, i] <- sqrt(result$var[1,1])
#                         N[i, j] <- N[j, i] <- n
#                         Test[i, j] <- pchisq(result$chisq, result$df, lower.tail=FALSE)
#                     }
#                     else R[i, j] <- R[j, i] <- result
#                 }
#             }
#         }
#         if (pd && !any(is.na(R)) && min(eigen(R, only.values=TRUE)$values) < 0){
#             cor <- Matrix::nearPD(R, corr=TRUE)
#             if (!cor$converged) warning("attempt to make correlation matrix positive-definite failed")
#             else warning("the correlation matrix has been adjusted to make it positive-definite")
#             R <- as.matrix(cor$mat)
#         }
#         rownames(R) <- colnames(R) <- names(data)
#         result <- list(correlations=R, type=Type, NA.method=use, ML=ML)
#         if (std.err) {
#             rownames(SE) <- colnames(SE) <- names(data)
#             rownames(N) <- colnames(N) <- names(N)
#             rownames(Test) <- colnames(Test) <- names(data)
#             result$std.errors <- SE
#             result$n <- if (use == "complete.obs") n else N
#             result$tests <- Test
#         }
#         class(result) <- "hetcor"
#         result
#     }
