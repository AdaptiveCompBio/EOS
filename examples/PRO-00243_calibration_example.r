# This file shows an example of finding the consistancy of calibration between replicates.  
# First, a patient sample (with replicates) match "Dx" sequences to consensus sequences
# Second, the consensus sequence are matched back to each replicate
# 
# See also PRO-00243
###########################################################################

library(EOS)
setwd("~/Documents/PRO-00243_data_prep/")

subject_files <- c("HGK5VBGX3_0_JandJ-Bald_FFPE4_O1_I2_L1_a.adap.txt.results.tsv.gz", 
                   "HGK5VBGX3_0_JandJ-Bald_FFPE4_O1_I2_L1_b.adap.txt.results.tsv.gz", 
                   "HGK5VBGX3_0_JandJ-Bald_FFPE4_O2_I1_L1_a.adap.txt.results.tsv.gz", 
                   "HGK5VBGX3_0_JandJ-Bald_FFPE4_O2_I1_L1_b.adap.txt.results.tsv.gz", 
                   "HGK5VBGX3_0_JandJ-Bald_FFPE4_O2_I2_L2_a.adap.txt.results.tsv.gz", 
                   "HGK5VBGX3_0_JandJ-Bald_FFPE4_O2_I2_L2_b.adap.txt.results.tsv.gz", 
                   "HGK7MBGX3_0_JandJ-Bald_FFPE4_O1_I1_L2_a.adap.txt.results.tsv.gz", 
                   "HGK7MBGX3_0_JandJ-Bald_FFPE4_O1_I1_L2_b.adap.txt.results.tsv.gz")

sample_ids <- c("FFPE4_O1_I2_L1_a", "FFPE4_O1_I2_L1_b",
                "FFPE4_O2_I1_L1_a", "FFPE4_O2_I1_L1_b",
                "FFPE4_O2_I2_L2_a", "FFPE4_O2_I2_L2_b",
                "FFPE4_O1_I1_L2_a", "FFPE4_O1_I1_L2_b")

Datain <- getReportSeqsByTag(subject_files, directory="Raw/", tag = "Dx")
obj <- new("EOS_Report", data=Datain)
obj <- processReport(obj, seq_edge_trim=2)
consensus <- findConsensusSequences(obj, no_shm_allowance=T, max_adj_mut_igkl=0, max_adj_mut_igh=2, ndn_cost = 1.05)

calibrated_list <- lapply(1:length(subject_files), function(i) {
  Datain <- getReportSeqsByTag(subject_files[i], directory="Raw/", tag = "Dx")
  obj <- new("EOS_Report", data=Datain)
  obj <- processReport(obj)
  matches <- matchReportToConsensus(obj, consensus)
  is_calibrated <- ifelse(names(consensus) %in% matches$Consensus_Sequence, "yes", "no")
  data.frame(Subject_Id = sample_ids[i], is_calibrated = is_calibrated, Calibrating_Sequence = names(consensus))
})
calibration_table <- do.call(rbind, calibrated_list)
write.table(calibration_table, sep="\t", file="FFPE4_calibration.tsv", row.names=F)

###########################################################################
# Example function for extracting meta data from EOS report files

getReportMetadata("HGK5VBGX3_0_JandJ-Bald_FFPE4_O1_I2_L1_a.adap.txt.results.tsv.gz", 
                  directory="Raw/", property = "estTotalNucleatedCells")
# [1] "8542.05"


