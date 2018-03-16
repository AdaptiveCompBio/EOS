library(dplyr)
library(EOS) # devtools::install_git("http://bitbucket.dna.corp.adaptivebiotech.com:7990/scm/~tching/eos.git")
library(trqwe)

setwd("~/Documents/PRO-00161")

Datain <- read.table("processed_data/20180130_PRO-00161_long-term-storage-stability_processed-data.tsv", 
                     sep="\t", stringsAsFactors=F, header=T)
Datain <- Datain %>% group_by(Sample.ID, Time.Point..Months.in.Storage., Sample.Source) %>% 
  summarize(Normalized_MRD=unique(Normalized.Sample.MRD.Frequency)) %>% as.data.frame

for(st in c("PB", "BMA")) {
  data <- subset(Datain, Sample.Source == st)
  
  # stability_model_selection function returns shelf life, ANCOVA results table and plot
  results <- with(data, stability_model_selection(Normalized_MRD, Time.Point..Months.in.Storage., Sample.ID, UL=1.3, LL=0.7))
  
  # Decorate the plot a little bit with labels etc.
  p <- results$stability_plot + labs(title=sprintf("Frozen Sample Stability: %s", st)) +
    geom_text(aes(x=0, y=1.55, label = paste("Shelf Life =", round(results$shelf_life, 2),"months")), size=5, hjust=0)
  
  #Save plot results
  ggsave(p,file=paste0("Results/frozen_stability_",st,".png"), width=7, height=4, dpi=300)
  
  # Save ANCOVA results
  g <- trqwe:::drew_table_plot(results$table)
  ggsave(g, file=sprintf("Results/%s_ancova.png", st), width=8, height=5)
}
