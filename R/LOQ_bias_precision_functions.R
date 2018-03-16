#' Function for calculating LOQ
#' @param MRD The measured MRD values
#' @param input_cancer_cells The pre-specified input.  Note: can be frequency or any other input value
#' @param prob Detection probability for LOD, default 0.95
#' @param conf.level Confidence level for calculating upper and lower bounds, default 0.95
#' @param log Boolean for log10-transforming input_cancer_cells, default TRUE
#' @param log_offset log offset for transformation to avoid log(0), default 1e-6
#' @param TE.plot  If TRUE, returns a ggplot of the Total Error fit
#' @param sadler_params  A list of initial values for sadler's fit
#' @param extra_sadler_log  Whether to take an extra log to improve Sadler fit stability; i.e.: log(y) ~ log(b1+b2*x)^b3)
#' @return Returns a list containing LOQ, LCL.LOQ and UCL.LOQ.  Note that LCL.LOQ will almost always be NA; use LOD as the lower limit. If TE.plot == TRUE, also contains TE_fit, fitted_data and TE_plot. 
find_LOQ <- function(MRD, input_cancer_cells, threshold = 0.7, conf.level=0.95, log=TRUE, log_offset=1e-6, TE.plot=TRUE, sadler_params=list(b1=2,b2=.1,b3=-1.8), extra_sadler_log=F) {
  bias_results <- calc_bias(MRD=MRD, input_cancer_cells=input_cancer_cells, bias.plot=F, log=log, log_offset=log_offset)$results
  precision_results <- calc_precision(MRD=MRD, input_cancer_cells=input_cancer_cells, precision.plot=F, log=log, 
                                     log_offset=log_offset, sadler_params=sadler_params)$results
  
  results <- merge(bias_results, precision_results, by=c("input_cancer_cells", "n", "mean_MRD"))
  results$TE <- sqrt(results$bias^2+results$sd_MRD^2)
  results$rel.TE <- sqrt(results$bias^2+results$sd_MRD^2) / results$input_cancer_cells
  results <- results[complete.cases(results),]

  if(log) {
    pdata <- data.frame(x=log10(results$input_cancer_cells + log_offset), y=results$TE)
  } else {
    pdata <- data.frame(x=results$input_cancer_cells, y=results$TE)
  }
  if(extra_sadler_log==F) {
    fit.nls <- repeated_nlsLM(y ~ (b1+b2*x)^b3, start=list(b1=sadler_params$b1,b2=sadler_params$b2,b3=sadler_params$b3),
                              data=pdata, control = list(maxiter = 1000))
  } else {
    fit.nls <- repeated_nlsLM(log(y) ~ log((b1+b2*x)^b3), start=list(b1=sadler_params$b1,b2=sadler_params$b2,b3=sadler_params$b3),
                              data=pdata, control = list(maxiter = 1000))
    predict <- function(fit.nls, ...) {
      exp(stats::predict(fit.nls,...))
    }
  }
  critval <- qnorm(conf.level/2 + 0.5)
  se.fit <- summary(fit.nls)$sigma
  umin <- min(pdata$x)
  umax <- max(pdata$x)
  f1 <- function(x, threshold = threshold) {
    pred <- predict(fit.nls,data.frame(x=x))
    if(log) x <- 10^x - log_offset
    pred / x - threshold
  }
  f2 <- function(x, critval = critval, threshold = threshold) {
    pred <- predict(fit.nls,data.frame(x=x))+critval*se.fit
    if(log) x <- 10^x - log_offset
    pred / x - threshold
  }
  f3 <- function(x, critval = critval,  threshold = threshold) {
    pred <- predict(fit.nls,data.frame(x=x))-critval*se.fit
    if(log) x <- 10^x - log_offset
    pred / x - threshold
  }
  LOQ <- rootSolve::uniroot.all(f1,c(umin, umax), threshold = threshold)[1]
  UCL.LOQ <- rootSolve::uniroot.all(f2,c(umin, umax*2), critval = critval, threshold = threshold)[1]
  LCL.LOQ <- rootSolve::uniroot.all(f3,c(umin, umax), critval = critval, threshold = threshold)[1]
  if(log) {
    LOQ <- 10^LOQ - log_offset
    LCL.LOQ <- 10^LCL.LOQ - log_offset
    UCL.LOQ <- 10^UCL.LOQ - log_offset
  }
  if(TE.plot) {
    newdata = data.frame(x = seq(min(pdata$x), max(pdata$x), length.out=1000))
    newdata$TE <- predict(fit.nls, newdata=newdata)
    newdata$lower <- newdata$TE - critval*se.fit
    newdata$upper <- newdata$TE + critval*se.fit
    newdata$rel.TE <- newdata$TE/10^(newdata$x + log_offset)
    newdata$rel.TE.lower <- newdata$lower/10^(newdata$x + log_offset)
    newdata$rel.TE.upper <- newdata$upper/10^(newdata$x + log_offset)
    if(log) newdata$x <- 10^(newdata$x) - log_offset
    p <- ggplot() + 
      geom_point(data = results, mapping = aes(x = input_cancer_cells, y = rel.TE), size = 1, col="red") +
      geom_line(data = newdata, mapping = aes(x=x, y=rel.TE)) + 
      scale_y_continuous(labels = scales::percent)
    if(log) {
      log10Fun <- function(x){log10(x + log_offset)}
      log10Fun_inv <- function(x){10^(x) - log_offset}
      break_min <- round(log10Fun(min(results$input_cancer_cells)))
      break_max <- round(log10Fun(max(results$input_cancer_cells)))
      p <- p + scale_x_continuous(trans = trans_new("log10Offset", log10Fun, log10Fun_inv), breaks=log10Fun_inv(seq(break_min, break_max, by=1)) + log_offset)
    }
  }
  ret <- list(LOQ = LOQ, LCL.LOQ = LCL.LOQ, UCL.LOQ = UCL.LOQ, results=results)
  if(TE.plot) ret$TE_plot <- p
  if(TE.plot) ret$fitted_data <- newdata
  if(TE.plot) ret$TE_fit <- fit.nls
  return(ret)
}

#' Function for calculating LOQ
#' @param MRD The measured MRD values
#' @param input_cancer_cells The pre-specified input.  Note: can be frequency or any other input value
#' @param bias.plot  If TRUE, returns a ggplot of the bias
#' @param log Boolean for log10-transforming input_cancer_cells, default TRUE
#' @param log_offset log offset for transformation to avoid log(0), default 1e-6
#' @return Returns a list containing results: a data.frame with bias and relative bias.  If bias.plot == TRUE, also contains a ggplot of bias vs. input cancer cells. 
calc_bias <- function(MRD, input_cancer_cells, bias.plot=TRUE, log=TRUE, log_offset=1e-6) {
  # data %>% group_by(input_cancer_cells) %>% summarize(mean=mean(MRD), sd=sd(MRD), var=var(MRD)) %>% mutate(bias=mean-input_cancer_cells, rel.bias=bias/input_cancer_cells)
  results <- lapply(unique(input_cancer_cells), function(icc) {
    mi <- MRD[input_cancer_cells == icc]
    # if(length(mi) < 2) stop(paste0("no replicates for input_cancer_cells: ", icc))
    mean_MRD <- mean(mi)
    bias <- mean_MRD - icc
    rel.bias <- bias/icc
    return(data.frame(input_cancer_cells=icc, n=length(mi), mean_MRD=mean_MRD, bias=bias, rel.bias=rel.bias))
  })
  results <- as.data.frame(data.table::rbindlist(results))
  if(bias.plot) {
    p <- ggplot() + 
      geom_point(data = results, mapping = aes(x = input_cancer_cells, y = rel.bias), size = 1, col="red") +
      geom_hline(yintercept = 0, col = "black", lty = 2, lwd = 1)  +
      scale_y_continuous(labels = scales::percent)
    if(log) {
      log10Fun <- function(x){log10(x + log_offset)}
      log10Fun_inv <- function(x){10^(x) - log_offset}
      break_min <- round(log10Fun(min(results$input_cancer_cells)))
      break_max <- round(log10Fun(max(results$input_cancer_cells)))
      p <- p + scale_x_continuous(trans = trans_new("log10Offset", log10Fun, log10Fun_inv), breaks=log10Fun_inv(seq(break_min, break_max, by=1)) + log_offset)
    }
    return(list(results=results, bias_plot=p))
  }
  return(list(results=results))
}

#' Function for calculating MRD precision
#' @param MRD The measured MRD values
#' @param input_cancer_cells The pre-specified input.  Note: can be frequency or any other input value
#' @param precision.plot  If TRUE, returns a ggplot of the Precision and Sadler's fit
#' @param log Boolean for log10-transforming input_cancer_cells, default TRUE
#' @param log_offset log offset for transformation to avoid log(0), default 1e-6
#' @param precision_type The type of precision.  Either percent CV ("CV") or relative S.D. ("rel_SD") , default "CV"
#' @param sadler_params  A list of initial values for sadler's fit
#' @return Returns a list containing LOQ, LCL.LOQ and UCL.LOQ.  Note that LCL.LOQ will almost always be NA; use LOD as the lower limit. If TE.plot == TRUE, also contains TE_fit, fitted_data and TE_plot. 
calc_precision <- function(MRD, input_cancer_cells, precision.plot=TRUE, log=TRUE, log_offset=1e-6, precision_type="CV", sadler_params=list(b1=2,b2=.1,b3=-1.8)) {
  results <- lapply(unique(input_cancer_cells), function(icc) {
    mi <- MRD[input_cancer_cells == icc]
    if(length(mi) < 2) stop(paste0("no replicates for input_cancer_cells: ", icc))
    mean_MRD <- mean(mi)
    sd_MRD <- sd(mi)
    if(precision_type == "CV") {
      precision = sd_MRD/mean_MRD
    } else {
      precision = sd_MRD/icc
    }
    return(data.frame(input_cancer_cells=icc, n=length(mi), mean_MRD=mean_MRD, sd_MRD=sd_MRD, precision=precision))
  })
  results <- as.data.frame(data.table::rbindlist(results))
  if(precision.plot) {
    if(log) {
      pdata <- data.frame(x=log10(results$input_cancer_cells + log_offset), y=results$precision)
    } else {
      pdata <- data.frame(x=results$input_cancer_cells, y=results$sd_MRD)
    }
    fit.nls = repeated_nlsLM(y ~ (b1+b2*x)^b3, start=list(b1=sadler_params$b1,b2=sadler_params$b2,b3=sadler_params$b3),data=pdata, control = list(maxiter = 1000))
    newdata = data.frame(x = seq(min(pdata$x), max(pdata$x), length.out=1000))
    newdata$y = predict(fit.nls, newdata = newdata, type="response", se=FALSE)
    if(log) newdata$x <- 10^(newdata$x) - log_offset
    p <- ggplot() + 
      geom_point(data = results, mapping = aes(x = input_cancer_cells, y = precision), size = 1, col="red") +
      geom_line(data = newdata, mapping = aes(x=x, y=y))
      scale_y_continuous(labels = scales::percent)
    if(log) {
      log10Fun <- function(x){log10(x + log_offset)}
      log10Fun_inv <- function(x){10^(x) - log_offset}
      break_min <- round(log10Fun(min(results$input_cancer_cells)))
      break_max <- round(log10Fun(max(results$input_cancer_cells)))
      p <- p + scale_x_continuous(trans = trans_new("log10Offset", log10Fun, log10Fun_inv), breaks=log10Fun_inv(seq(break_min, break_max, by=1)) + log_offset)
    }
    return(list(results=results, precision_plot=p, precision_fit = fit.nls))
  }
  return(list(results=results))
}

repeated_nlsLM <- function(..., start, niter=100, seed=7) {
  set.seed(seed)
  best_fit <- NULL
  best_rse <- Inf
  max_param <- abs(max(unlist(start)))
  
  try({
    fit <- minpack.lm::nlsLM(start=start, ...)
    rse <- sum(resid(fit)^2)
    if(rse < best_rse) {
      best_fit <- fit
      best_rse <- rse
    }
  }, silent=T)
  
  for(i in 2:niter) {
    start_i <- lapply(start, function(sp) rnorm(1, mean=sp, sd=max_param*10))
    try({
      fit <- minpack.lm::nlsLM(start=start_i, ...)
      rse <- sum(resid(fit)^2)
      if(rse < best_rse) {
        best_fit <- fit
        best_rse <- rse
      }
    }, silent=T)
  }
  if(is.null(best_fit)) stop("No possible fit for minpack::nlsLM")
  return(best_fit)
}



