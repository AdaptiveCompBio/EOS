#' Function for selecting stability models
#' @param MRD The measured MRD values
#' @param Time Time (numeric)
#' @param Sample.ID Sample ID used for model selection
#' @param UL Upper limit of MRD (in MRD raw scale)
#' @param LL Lower limit of MRD (in MRD raw scale)
#' @param stability.plot  If TRUE, returns a ggplot of the data
#' @return Returns a list containing a table of the ANCOVA analyses, the selected model, the estimated shelf_life and a plot if specified.  
stability_model_selection <- function(MRD, Time, Sample.ID, stability.plot=TRUE, UL = 1.3, LL = 0.7) {
  model_labels <- c("Separate Intercept / Separate Slopes vs.\nCommon Intercept / Common Slope\n(SISS-CICS)",
             "Separate Intercepts / Common Slope vs.\nCommon Intercept / Common Slope\n(SICS-CICS)",
             "Separate Intercepts / Separate Slopes vs.\nSeparate Intercepts / Common Slope\n(SISS-SICS)",
             "Whole Model")
  data <- data.frame(MRD=as.numeric(MRD), Time=as.numeric(Time), Sample.ID=factor(Sample.ID))
  FIT_SISS <- lm(MRD ~ -1 + Sample.ID + Sample.ID:Time, data = data) 
  FIT_RESIDUAL <- lm(MRD ~ -1 + Time, data = data)
  FIT_NULL <- lm(MRD ~ -1, data = data)
  FIT_SICS <- lm(MRD ~ -1 + Sample.ID + Time, data = data)
  FIT_CICS <- lm(MRD ~ Time, data = data)
  ANOVA_EQ_SLOPE <- as.data.frame( anova(FIT_SICS, FIT_SISS) )[2,3:6]
  ANOVA_EQ_INTERCEPT <- as.data.frame( anova(FIT_CICS, FIT_SICS) )[2,3:6]
  ANOVA_SISS_CICS <- as.data.frame( anova(FIT_CICS, FIT_SISS) )[2,3:6]
  ANOVA_WHOLE_MODEL <- as.data.frame( anova(FIT_SISS) )[3,]
  ANOVA_WHOLE_MODEL <- data.frame(Df=ANOVA_WHOLE_MODEL$Df, `Sum of Sq`=ANOVA_WHOLE_MODEL$`Sum Sq`, F=NA, `Pr(>F)`=NA, check.names=F)
  Table <- rbind(ANOVA_SISS_CICS, ANOVA_EQ_INTERCEPT, ANOVA_EQ_SLOPE, ANOVA_WHOLE_MODEL)
  Table <- cbind(Model=model_labels, Table)
  colnames(Table)[c(3,5)] <- c("SS", "p.value")
  if (ANOVA_EQ_INTERCEPT$`Pr(>F)` >= 0.25 & ANOVA_EQ_SLOPE$`Pr(>F)` >= 0.25) {
    model <- "Common Intercept / Common Slope"
    model_fun <- CICS_model
  } else if (ANOVA_EQ_SLOPE$`Pr(>F)` >= 0.25 & ANOVA_EQ_INTERCEPT$`Pr(>F)` < 0.25) {
    model <- "Separate Intercepts / Common Slope"
    model_fun <- SICS_model
  } else if (ANOVA_EQ_SLOPE$`Pr(>F)` < 0.25) {
    model <- "Separate Intercepts / Separate Slopes"
    model_fun <- SISS_model
  } else {
    stop("\nCould not choose model, something is wrong \n")
  }
  results <- with(data, model_fun(MRD, Time, Sample.ID, stability.plot=stability.plot, UL = UL, LL = LL))
  results$model <- model
  results$table <- Table
  return(results)
}

#' Separate Intercepts Separate Slopes
#' Separate Intercepts Separate Slopes for stability analysis
#' @param MRD The measured MRD values
#' @param Time Time (numeric)
#' @param Sample.ID Sample ID used for model selection
#' @param UL Upper limit of MRD (in MRD raw scale)
#' @param LL Lower limit of MRD (in MRD raw scale)
#' @param stability.plot  If TRUE, returns a ggplot of the data
#' @return Returns a list the estimated shelf_life and a plot if specified.  
SISS_model <- function(MRD, Time, Sample.ID, UL = 1.3, LL = 0.7, stability.plot=TRUE) {
  data <- data.frame(MRD=as.numeric(MRD), Time=as.numeric(Time), Sample.ID=factor(Sample.ID))
  FIT <- lm(MRD ~ -1 + Sample.ID + Sample.ID:Time, data = data)
  k <- length(unique(data$Sample.ID))
  shelf_life <- Inf
  intercept_range <- c(0, max(Time)*10)
  for(i in seq_len(k)){
    sid <- unique(data$Sample.ID)[i]
    upper <- function(x) {
      predict.lm(FIT, data.frame(Sample.ID = sid, Time = x), interval="confidence")[,"upr"] - UL
    }
    lower <- function(x) {
      LL - predict.lm(FIT, data.frame(Sample.ID = sid, Time = x), interval="confidence")[,"lwr"]
    }
    Upper.Int <- tail(rootSolve::uniroot.all(upper,intercept_range),1)
    Lower.Int <- tail(rootSolve::uniroot.all(lower,intercept_range),1)
    shelf_life <- min(shelf_life, Lower.Int, Upper.Int)
  }
  if(!stability.plot) return(list(shelf_life=shelf_life))
  xseq <- seq(0,max(shelf_life, data$Time),0.05)
  df <- data.table::rbindlist(lapply(seq_len(k), function(i) {
    sid <- unique(data$Sample.ID)[i]
    data.frame(predict.lm(FIT, data.frame(Time=xseq, Sample.ID=sid), interval="confidence"), Sample.ID=sid, x=xseq)
  }))
  p <- stability_base_plot(data, "Separate Intercepts / Separate Slopes", shelf_life) + 
    geom_ribbon(data=df, mapping=aes(x=x, ymin=lwr, ymax=upr, fill=Sample.ID), alpha=0.3, show.legend=F)
  return(list(shelf_life=shelf_life, stability_plot=p))
}

#' Separate Intercepts Common Slope
#' Separate Intercepts Common Slopes for stability analysis
#' @param MRD The measured MRD values
#' @param Time Time (numeric)
#' @param Sample.ID Sample ID used for model selection
#' @param UL Upper limit of MRD (in MRD raw scale)
#' @param LL Lower limit of MRD (in MRD raw scale)
#' @param stability.plot  If TRUE, returns a ggplot of the data
#' @return Returns a list the estimated shelf_life and a plot if specified.  
SICS_model <- function(MRD, Time, Sample.ID, UL = 1.3, LL = 0.7, stability.plot=TRUE) {
  data <- data.frame(MRD=as.numeric(MRD), Time=as.numeric(Time), Sample.ID=factor(Sample.ID))
  FIT <- lm(MRD ~ -1 + Sample.ID + Time, data = data)
  k <- length(unique(data$Sample.ID))
  shelf_life <- Inf
  intercept_range <- c(0, max(Time)*10)
  for(i in seq_len(k)){
    sid <- unique(data$Sample.ID)[i]
    upper <- function(x) {
      predict.lm(FIT, data.frame(Sample.ID = sid, Time = x), interval="confidence")[,"upr"] - UL
    }
    lower <- function(x) {
      LL - predict.lm(FIT, data.frame(Sample.ID = sid, Time = x), interval="confidence")[,"lwr"]
    }
    Upper.Int <- tail(rootSolve::uniroot.all(upper,intercept_range),1)
    Lower.Int <- tail(rootSolve::uniroot.all(lower,intercept_range),1)
    shelf_life <- min(shelf_life, Lower.Int, Upper.Int)
  }
  if(!stability.plot) return(list(shelf_life=shelf_life))
  xseq <- seq(0,max(shelf_life, data$Time),0.05)
  df <- data.table::rbindlist(lapply(seq_len(k), function(i) {
    sid <- unique(data$Sample.ID)[i]
    data.frame(predict.lm(FIT, data.frame(Time=xseq, Sample.ID=sid), interval="confidence"), Sample.ID=sid, x=xseq)
  }))
  p <- stability_base_plot(data, "Separate Intercepts / Common Slope", shelf_life) + 
    geom_ribbon(data=df, mapping=aes(x=x, ymin=lwr, ymax=upr, fill=Sample.ID), alpha=0.3, show.legend=F) 
  return(list(shelf_life=shelf_life, stability_plot=p))
}

#' Common Intercepts Common Slope
#' Common Intercepts Common Slopes for stability analysis
#' @param MRD The measured MRD values
#' @param Time Time (numeric)
#' @param Sample.ID Sample ID used for model selection
#' @param UL Upper limit of MRD (in MRD raw scale)
#' @param LL Lower limit of MRD (in MRD raw scale)
#' @param stability.plot  If TRUE, returns a ggplot of the data
#' @return Returns a list the estimated shelf_life and a plot if specified.  
CICS_model <- function(MRD, Time, Sample.ID, UL = 1.3, LL = 0.7, stability.plot=TRUE) {
  data <- data.frame(MRD=as.numeric(MRD), Time=as.numeric(Time), Sample.ID=factor(Sample.ID))
  FIT <- lm(MRD ~ Time, data = data)
  intercept_range <- c(0, max(Time)*10)
  upper <- function(x) {
    predict.lm(FIT, data.frame( Time = x), interval="confidence")[,"upr"] - UL
  }
  lower <- function(x) {
    LL - predict.lm(FIT, data.frame(Time = x), interval="confidence")[,"lwr"]
  }
  Upper.Int <- tail(rootSolve::uniroot.all(upper,intercept_range),1)
  Lower.Int <- tail(rootSolve::uniroot.all(lower,intercept_range),1)
  shelf_life <- min(Lower.Int, Upper.Int)
  if(!stability.plot) return(list(shelf_life=shelf_life))
  xseq <- seq(0,max(shelf_life, data$Time),0.05)
  df <- data.frame(predict.lm(FIT, data.frame(Time=xseq), interval="confidence"), x=xseq)
  p <- stability_base_plot(data, "Common Intercept / Common Slope", shelf_life) + 
    geom_ribbon(data=df, mapping=aes(x=x, ymin=lwr, ymax=upr), alpha=0.3, show.legend=F)
  return(list(shelf_life=shelf_life, stability_plot=p))
}

stability_base_plot <- function(data, model, shelf_life) {
  ggplot() + 
  ggthemes::theme_stata() +
  geom_point(data=data, mapping=aes(x=Time, y=MRD, group = Sample.ID, color = Sample.ID)) +
  labs(subtitle=model) +
  coord_cartesian(xlim=c(-0.05,max(shelf_life, data$Time)), expand=T) + ylim(0,NA) + 
  geom_vline(xintercept=shelf_life, col = "red", lty = 2, lwd = 1.52) +
  geom_hline(yintercept=1.3, col = "black", lty = 2, lwd = 1.5) +
  geom_hline(yintercept=0.7, col = "black", lty = 2, lwd = 1.5) +  
  geom_hline(yintercept=1, col = "black", lty = 1, lwd = 1)
}