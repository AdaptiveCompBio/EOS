#' Cox-PH model assumptions
#' This function outputs a bunch of plots related to Cox-PH diagnostics and model assumptions.  Based on "Amendment 2 DFCI Analyses 3November2017.R"
#' @param model A cox-ph fitted model object
#' @param data The data.frame containing the covariates used to generate the model
#' @param continuous_vars A character vector of the names of the continuous variables used in the model
#' @param save_prefix The file prefix for saving the output.  E.g., "Results/coxph_"
#' @return This function doesn't return anything, but writes four plots to disk: the scaled schoenfeld residuals, dfBetas, Deviance and Martingale residuals.  These plots can be used to e.g. diagnose issues with the model, find outliers or find time-varying trends.  
#' @examples
#' data(veteran)
#' model <- coxph(Surv(time, status) ~ ., data=veteran)
#' cvars <- c("karno", "diagtime", "age")
#' CoxPH_Diagnostics(model, data=veteran, continuous_vars=cvars, save_prefix="~/Desktop/coxph_veteran_")
#' @seealso 
#' "Checking the Cox model with cumulative sums of martingale-based residuals." Biometrika 80.3 (1993): 557-572.
#' "Using Time Dependent Covariates..." (https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf)
#' "Proportional Hazards Regression Diagnostics" (http://www.ics.uci.edu/~dgillen/STAT255/Handouts/lecture10.pdf)
CoxPH_Diagnostics <- function(model, data, continuous_vars, save_prefix) {
  
  vars <- all.vars(model$formula)
  data <- data[,vars]
  data <- data[complete.cases(data),]

# Grapical testing the proportional hazards assumption: Scaled Schoenfeld Residuals
p <- survminer::ggcoxzph(survival::cox.zph(model), ggtheme = ggthemes::theme_stata(),xlab = "Time")
save_file <- paste0(save_prefix, "schoenfeld.png")
png(save_file, units='in', width=10, height=10, res=96)
print(p)
dev.off()

# Graphical test for influential observations: dfBetas
p <- survminer::ggcoxdiagnostics(model, type = "dfbeta", linear.predictions = T, ggtheme = ggthemes::theme_stata(), ox.scale="observation.id")
save_file <- paste0(save_prefix, "dfbeta.png")
png(save_file, units='in', width=10, height=10, res=96)
print(p)
dev.off()

# Graphical test for influential observations: Deviance
p <- survminer::ggcoxdiagnostics(model, type = "deviance", linear.predictions = T, title = "Deviance Residuals", 
                     subtitle = " ", ggtheme = ggthemes::theme_stata(), ox.scale="observation.id")
save_file <- paste0(save_prefix, "deviance.png")
png(save_file, units='in', width=10, height=10, res=96)
print(p)
dev.off()

#### TC: The martingale residual plot Drew used only accepted continuous predictors; but not all predictors are continuous.  
#### Let's use a more general approach

# Graphical test for nonlinearity
# p <- ggcoxfunctional(model, data=complete.data, point.size=2, point.col = "blue",  ggtheme = theme_stata())
# save_file <- paste0(save_prefix, "functionalform.png")
# png(save_file, units='in', width=10, height=10, res=96)
# print(p)
# dev.off()


if(length(continuous_vars) == 0) return(NULL) # can't do functional form if no continuous vars

# If right censored, use built in function from goftte library
# Otherwise manually plot out residuals
save_file <- paste0(save_prefix, "functionalform.png")
if(ncol(model$y) == 2) {
  fm <- goftte::fcov(model)
  png(save_file, units='in', width=10, height=10, res=96)
  par(mfrow=c(2,2))
  for(cvar in continuous_vars) {
    plot(fm, idx=which(fm$variable == cvar))
  }
  dev.off()
} else {
  resid <- residuals(model)
  png(save_file, units='in', width=10, height=10, res=96)
  par(mfrow=c(2,2))
  for(cvar in continuous_vars) {
    x <- data[,cvar]
    plot(x, resid, main = paste0("Martingale residuals: ", cvar), ylab="residuals", xlab=cvar)
  }
  dev.off()
}
}
