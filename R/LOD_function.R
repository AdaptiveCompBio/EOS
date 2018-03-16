#' Function for calculating LOD
#' @param MRD The measured MRD values
#' @param input_cancer_cells The pre-specified input.  Note: can be frequency or any other input value
#' @param prob Detection probability for LOD, default 0.95
#' @param conf.level Confidence level for calculating upper and lower bounds, default 0.95
#' @param log Boolean for log10-transforming MRD and input_cancer_cells, default TRUE
#' @param log_offset If log transforming, offset value for the MRD to avoid log(0) issues, default 0
#' @param probit.plot If TRUE, returns a ggplot of the probit fit
#' @return Returns a list containing LOD, LCL.LOD and UCL.LOD.  If probit.plot == TRUE, also contains probit_fit, fitted_data and probit_plot.   
#' @examples
#' Datain <- read.table("processed_data/20171119_MVP-00147_lod_data-table1.tsv", head=TRUE, as.is=TRUE, sep="\t")
#' find_LOD(Datain$Sample.MRD.Frequency, Datain$Input.Cancer.Cells, probit.plot=F)$LOD
#' [1] 2.289792
find_LOD <- function(MRD, input_cancer_cells, prob = 0.95, conf.level = 0.95, log=TRUE, log_offset=1e-6, probit.plot = TRUE) {
  detected <- ifelse(MRD > 0, 1, 0)
  if(all(detected == 1)) stop("No MRD Negatives")
  if(all(detected == 0)) stop("No MRD Positives")
  if(log) input_cancer_cells <- log10(input_cancer_cells + log_offset)
  data <- data.frame(detected=detected, input_cancer_cells=input_cancer_cells)

  fit <- glm(detected ~ input_cancer_cells, family = binomial(link="probit"), data=data, control = glm.control(epsilon=1e-16, maxit=1000))
  intercept <- coef(fit)[1]
  slope <- coef(fit)[2]
  
  ## optimization function for finding LOD
  f1 <- function(x, intercept=intercept, slope=slope, prob=prob) {
    pnorm(intercept + slope*x) - prob
  }
  
  ## optimization function for roots of lower confidence interval
  f2 <- function(x, intercept=intercept, slope=slope, prob=prob, conf.level=conf.level) {
    se.fit <- predict(fit, newdata=data.frame(input_cancer_cells=x), type="link", se=TRUE)$se.fit
    pnorm(intercept + slope*x + qnorm(conf.level/2 + 0.5)*se.fit) - prob
  }
  
  ## optimization function for roots of upper confidence interval
  f3 <- function(x, intercept=intercept, slope=slope, prob=prob, conf.level=conf.level) {
    se.fit <- predict(fit, newdata=data.frame(input_cancer_cells=x), type="link", se=TRUE)$se.fit
    pnorm(intercept + slope*x - qnorm(conf.level/2 + 0.5)*se.fit) - prob
  }
  
  lod <- uniroot(f1, c(-10,200), intercept=intercept, slope=slope, prob=prob)$root
  lcl.lod <- uniroot(f2, c(-10,200), intercept=intercept, slope=slope, prob=prob, conf.level=conf.level)$root
  ucl.lod <- uniroot(f3, c(-10,200), intercept=intercept, slope=slope, prob=prob, conf.level=conf.level)$root
  if(log) {
    lod <- 10^lod - log_offset
    lcl.lod <- 10^lcl.lod - log_offset
    ucl.lod <- 10^ucl.lod - log_offset
  }
  
  # Bare bones probit fit plot
  if(probit.plot){
    #	Calculate predicted values and confidence intervals
    newdata <- data.frame(input_cancer_cells = seq(min(input_cancer_cells), max(input_cancer_cells),length.out = 1000))
    predicted <- as.data.frame(predict(fit, newdata = newdata, type="link", se=TRUE))
    newdata <- cbind(newdata, predicted)
    predicted <- as.data.frame(predict(fit, newdata = newdata, type="link", se=TRUE))
    critval <- qnorm(conf.level/2 + 0.5)
    newdata$lower <- fit$family$linkinv(newdata$fit - critval*newdata$se.fit)
    newdata$upper <- fit$family$linkinv(newdata$fit + critval*newdata$se.fit)
    newdata$fit  <- fit$family$linkinv(newdata$fit)  # Rescale to 0-1
    if(log) {
      newdata$input_cancer_cells <- 10^(newdata$input_cancer_cells) - log_offset
      data$input_cancer_cells <- 10^(data$input_cancer_cells) - log_offset
    }
    
    p <- ggplot() + 
      geom_point(data = data, aes(x = input_cancer_cells, y = detected), size = 2) +
      geom_ribbon(data = newdata, mapping=aes(x=input_cancer_cells, ymin=lower, ymax=upper), alpha=0.5, fill = "darkgreen") +
      geom_line(data=newdata, aes(x=input_cancer_cells, y=fit), col = "black") + 
      geom_hline(yintercept=prob, lty=2, col="red", alpha=0.67) + 
      geom_vline(xintercept = c(lod, lcl.lod, ucl.lod), lty=2, col="red", alpha=0.67) + 
      scale_y_continuous(labels = scales::percent)
    if(log) {
      p <- p + scale_x_log10()
      # log10Fun <- function(x){log10(x + log_offset)}
      # log10Fun_inv <- function(x){10^(x) - log_offset}
      # break_min <- round(log10Fun(min(data$input_cancer_cells)))
      # break_max <- round(log10Fun(max(data$input_cancer_cells)))
      # p <- p + scale_x_continuous(trans = trans_new("log10Offset", log10Fun, log10Fun_inv), breaks=log10Fun_inv(seq(break_min, break_max, by=1)) + log_offset)
    }
  }
  
  ret <- list(LOD = lod, LCL.LOD = lcl.lod, UCL.LOD = ucl.lod)
  if(probit.plot) ret$probit_plot <- p
  if(probit.plot) ret$fitted_data <- newdata
  if(probit.plot) ret$probit_fit <- fit
  return(ret)
}
