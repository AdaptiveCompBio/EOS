#' Function for finding linear MRD range
#' @param MRD The measured MRD values
#' @param input_cancer_cells The pre-specified input.  Note: can be frequency or any other input value
#' @param max_deviance_threshold Maximum deviation threshold from linearity
#' @param log Boolean for log10-transforming input_cancer_cells, default TRUE
#' @param log_offset log offset for transformation to avoid log(0), default 1e-6
#' @param linearity.plot  If TRUE, returns a ggplot of the Total Error fit
#' @param log_deviation  If TRUE and if log==T, calculate deviation in log space.  Otherwise, convert deviance to linear space.  
#' @return Returns a list containing linear_range_min, linear_range_max and deviation within the linear range.  If linearity.plot == TRUE, also contains a ggplot of the linear range. 
#' @seealso PRO-00091
find_linearity <- function(MRD, input_cancer_cells, max_deviation_threshold=0.05, log=TRUE, log_offset=1e-6, linearity.plot=TRUE, log_deviation=T) {
  # input_cancer_cells <- Datain$Input.Cancer.Cells
  # MRD <- Datain$Sample.MRD.Frequency
  ord <- order(input_cancer_cells)
  MRD <- MRD[ord]
  input_cancer_cells <- input_cancer_cells[ord]
  while(length(input_cancer_cells) > 3) {
    deviation <- calc_deviation(MRD, input_cancer_cells, log, log_offset, log_deviation)
    nl_points <- abs(deviation) > max_deviation_threshold
    if(all(!nl_points)) {
      break;
    } else {
      # PRO-00091: If any DLi is greater than the stated criterion, then clone frequencies at extremes 
      # where DLi is too extreme will be removed, thereby reducing the linear range of the assay, 
      # and re-analysis performed until linearity is achieved.
      # There is no instructions if deviation > 0.05 only occurs at non-extrema (i.e. in the middle).  
      # Assumption: if deviation only in middle, chew from the bottom end
      input_nonlinear_points <- unique(input_cancer_cells[nl_points])
      input_nonlinear_points_extrema <- intersect(input_nonlinear_points, range(input_cancer_cells))
      if(length(input_nonlinear_points_extrema) != 0) {
        MRD <- MRD[!input_cancer_cells %in% input_nonlinear_points_extrema]
        input_cancer_cells <- input_cancer_cells[!input_cancer_cells %in% input_nonlinear_points_extrema]
      } else {
        MRD <- MRD[!input_cancer_cells %in% min(input_cancer_cells)]
        input_cancer_cells <- input_cancer_cells[!input_cancer_cells %in% min(input_cancer_cells)]
      }
    }
  }
  
  ret <- list(linear_range_min=min(input_cancer_cells), linear_range_max=max(input_cancer_cells), max_deviation=max(abs(deviation)))
  if(linearity.plot) {
    if(log) {
      df <- data.frame(x=log10(input_cancer_cells+log_offset), y=log10(MRD+log_offset))
      y <- log10(MRD + log_offset)
      x <- log10(input_cancer_cells + log_offset)
    } else {
      df <- data.frame(x=input_cancer_cells, y=MRD)
      y <- MRD
      x <- input_cancer_cells
    }
    lin_fit <- lm(y ~ x)
    quad_fit <- lm(y ~ x + I(x^2))
    cubic_fit <- lm(y ~ x + I(x^2) + I(x^3))
    xseq <- data.frame(x=seq(min(x), max(x), length.out = 1000))
    xseq$linear <- predict(lin_fit, newdata=xseq)
    xseq$quadratic <- predict(quad_fit, newdata=xseq)
    xseq$cubic <- predict(cubic_fit, newdata=xseq)
    
    if(linearity.plot) {
      p <- ggplot() + geom_point(data=df, aes(x=x, y=y)) + geom_line(data=xseq, aes(x=x, y=linear)) + geom_line(data=xseq, aes(x=x, y=cubic))
      ret$linearity_plot <- p
      ret$plot_data <- xseq
    }
    ret$lin_fit <- lin_fit
    ret$quad_fit <- quad_fit
    ret$cubic_fit <- cubic_fit
  }
  return(ret)
}


calc_deviation <- function(MRD, input_cancer_cells, log, log_offset, log_deviation) {
  if(log) {
    y <- log10(MRD + log_offset)
    x <- log10(input_cancer_cells + log_offset)
  } else {
    y <- MRD
    x <- input_cancer_cells
  }
  lin_fit <- lm(y ~ x)
  quad_fit <- lm(y ~ x + I(x^2))
  cubic_fit <- lm(y ~ x + I(x^2) + I(x^3))
  poly_p_values <- c(coef(summary(quad_fit))[,"Pr(>|t|)"][c(-1,-2)],coef(summary(cubic_fit))[,"Pr(>|t|)"][c(-1,-2)])
  
  # PRO-00091: "If any of the nonlinear coefficients, b2 or b3, is significant (p < 0.05) then signal response is nonlinear and analyses proceed as follows. For the selected reference standard, the degree of nonlinearity in signal response will be assessed by examining the standard error of the regression and selecting the higher-order (nonlinear) polynomial model with the best fit"
  # The lowest standard error will always be cubic, hence if any polynomial coefficient is significant, deviance will always be determined through the cubic fit
  
  if(any(poly_p_values < 0.05)) {
    if(log & log_deviation) {
      deviation <- (10^predict(cubic_fit) - 10^predict(lin_fit)) / (10^predict(lin_fit)+log_offset)
    } else {
      deviation <- (predict(cubic_fit) - predict(lin_fit)) / predict(lin_fit) # note: denominator at 100% input: log(1) = 0 gives kind of wacky results
    }
  } else {
    deviation <- rep(0, length(x))
  }
  return(deviation)
}
