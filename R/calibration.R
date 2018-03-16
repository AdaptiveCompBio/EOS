report_col_classes <- structure(c("character","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character","character","character","character","character","character","character","character","character","character","character","character","character","character","character","character","character","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character","character","character","character","character","character","character","character","character","character","character","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character","character","character","character","character","character","character","character","character","character","character","character","numeric","numeric","character","numeric","character","numeric","numeric","numeric","character","character","numeric"),.Names=c("nucleotide","aminoAcid","copy","copyNormalized","count","frequency","frequencyNormalized","frequencyCount","cdr3Length","vFamilyName","vGeneName","vGeneAllele","vFamilyTies","vGeneNameTies","vGeneAlleleTies","dFamilyName","dGeneName","dGeneAllele","dFamilyTies","dGeneNameTies","dGeneAlleleTies","jFamilyName","jGeneName","jGeneAllele","jFamilyTies","jGeneNameTies","jGeneAlleleTies","vDeletion","d5Deletion","d3Deletion","jDeletion","n2Insertion","n1Insertion","vIndex","n2Index","dIndex","n1Index","jIndex","vdNormalizationFactor","jNormalizationFactor","inputTemplateEstimate","frequencyInputTemplateEstimate","sequenceStatus","vMaxResolved","d2MaxResolved","dMaxResolved","jMaxResolved","cloneResolved","d2FamilyName","d2GeneName","d2GeneAllele","d2FamilyTies","d2GeneNameTies","d2GeneAlleleTies","vScore","dScore","jScore","vAlignLength","vAlignSubstitutionCount","dAlignLength","dAlignSubstitutionCount","jAlignLength","jAlignSubstitutionCount","vOrphon","dOrphon","jOrphon","vFunction","dFunction","jFunction","vAlignSubstitutionIndexes","dAlignSubstitutionIndexes","jAlignSubstitutionIndexes","vAlignSubstitutionGeneThreePrimeIndexes","dAlignSubstitutionGeneThreePrimeIndexes","jAlignSubstitutionGeneThreePrimeIndexes","overlapCount","overlapCountReads","locus","cloneProbability","sequenceTags","ndnMutationWeight","maxAdjustedMutations","diseaseLoadMultiplier","cellularSensitivityBin","cellFreeSensitivityBin","inputTemplateEstimateFloat"))

#' Object for calibrating and matching sequences
#'
#' @slot trimmed_sequence trimmed sequence
#' @slot locus locus determined through getLocus function
#' @slot cdr3Index calculated cdr3 index
#' @slot cdr3Length calculated cdr3 length
#' @slot ndnIndex calculated ndn index
#' @slot ndnLength calculated ndn length
#' @slot left_sequence left split sequence based on cdr3 index
#' @slot right_sequence right split sequence based on cdr3 index
#' @slot left_ndn_start left ndn base start
#' @slot left_ndn_end left ndn base end
#' @slot right_ndn_start left ndn base start
#' @slot right_ndn_end left ndn base end
#' @slot left_hash tree hash of left sequence
#' @slot right_hash tree hash of right sequence
#' @slot data The original CSV report data
#' @examples
#' new("EOS_Report", data=Datain)
#' @seealso 
#' See examples/PRO-00243_calibration_example.r
#' @name EOS_Report-class
#' @rdname EOS_Report-class
#' @export
setClass("EOS_Report", representation(trimmed_sequence="character",
                                      locus = "character",
                                      cdr3Index = "numeric",
                                      cdr3Length = "numeric",
                                      ndnIndex = "numeric",
                                      ndnLength = "numeric",
                                      left_sequence = "character",
                                      right_sequence = "character",
                                      left_ndn_start="numeric",
                                      left_ndn_end="numeric", 
                                      right_ndn_start="numeric", 
                                      right_ndn_end="numeric",
                                      left_hash = "list",
                                      right_hash = "list",
                                      data = "data.frame"))

setMethod("initialize", "EOS_Report", function(.Object, ...) {
  .Object <- callNextMethod()
  for(i in seq_len(ncol(.Object@data))) {
    mode(.Object@data[[i]]) <- report_col_classes[colnames(.Object@data)[i]] # mode not class, because R
  }
  .Object
})

#' Function for calculating NDN base pairs
#' @param obj The EOS report object
#' @return The object with ndnIndex and ndnLength slots filled
setGeneric("determineNdnValues", function(obj) {standardGeneric("determineNdnValues")})

#' Function for calculating cdr3 index
#' @param obj The EOS report object
#' @return The object with cdr3Index slot filled
setGeneric("determineCdr3Index", function(obj) {standardGeneric("determineCdr3Index")})

#' Function for calculating gene locus
#' @param obj The EOS report object
#' @return The object with locus slot filled
setGeneric("getLocus", function(obj) {standardGeneric("getLocus")})

#' Function for generating trimmed sequences
#' @param obj The EOS report object
#' @param seq_edge_trim how many base pairs to trim on each side
#' @return The object with trimmed_sequence slot filled.  "N" base pairs are also removed.
setGeneric("trimCloneSequence", function(obj, seq_edge_trim) {standardGeneric("trimCloneSequence")})

#' Function for splitting sequence into subsequences
#' This function splits sequences into left and right sequences based on cdr3Index.  Must have cdr3Index, ndnIndex, ndnLength and trimmed_sequence slots filled.  Call determineNdnValues, determineCdr3Index and trimCloneSequence functions first.  
#' @param obj The EOS report object
#' @return The object with left_seq and right_seq slots filled.  
setGeneric("splitSequenceIntoSubstrings", function(obj) {standardGeneric("splitSequenceIntoSubstrings")})

#' Function for creating keyword trees
#' Create keyword trees used for finding compatible sequences between reports.  Call getLocus and spltiSequenceIntoSubstrings first.  
#' @param obj The EOS report object
#' @return The object with left_hash and right_hash slots filled.  
setGeneric("createKeywordTrees", function(obj) {standardGeneric("createKeywordTrees")})

#' determineNdnValues
#' 
setMethod("determineNdnValues", signature("EOS_Report"), function(obj) {
  n2Index <- obj@data$n2Index
  dIndex <- obj@data$dIndex
  n1Index <- obj@data$n1Index
  jIndex <- obj@data$jIndex
  
  ndnIndex <- n2Index
  ndnIndex <- ifelse(ndnIndex == -1, dIndex, ndnIndex)
  ndnIndex <- ifelse(ndnIndex == -1, n1Index, ndnIndex)
  
  select <- jIndex != -1 & jIndex > ndnIndex
  ndnIndex <- ifelse(!select, 0, dplyr::case_when(
    ndnIndex < -1 ~ 0,
    ndnIndex == -1 ~ 0,
    ndnIndex > -1 ~ ndnIndex
  ))
  ndnLength<- ifelse(!select, 0, dplyr::case_when(
    ndnIndex < -1 ~ jIndex,
    ndnIndex == -1 ~ 0,
    ndnIndex > -1 ~ jIndex - ndnIndex
  ))
  obj@ndnIndex <- ndnIndex
  obj@ndnLength <- ndnLength
  return(obj);
})

#' determineCdr3Index
setMethod("determineCdr3Index", signature("EOS_Report"), function(obj) {
  cloneLength <- nchar(obj@data$nucleotide)
  cdr3Index <- C_determineCdr3Index(cloneLength, obj@data$vIndex, obj@data$dIndex, obj@data$cdr3Length)
  cdr3Index <- ifelse(cdr3Index == -1, nchar(obj@data$nucleotide), cdr3Index)
  obj@cdr3Index <- cdr3Index
  return(obj)
})

#' getLocus
setMethod("getLocus", signature("EOS_Report"), function(obj) {
  obj@locus <- C_getLocus(obj@data$cloneResolved, obj@data$vFamilyName,
             obj@data$vFamilyTies, obj@data$dFamilyName,
             obj@data$dFamilyTies, obj@data$jFamilyName,
             obj@data$jFamilyTies)
  return(obj)
})

#' trimCloneSequence
setMethod("trimCloneSequence", signature("EOS_Report", "numeric"), function(obj, seq_edge_trim=2) {
  clone_sequence <- obj@data$nucleotide
  cdr3Index <- obj@cdr3Index
  ndnIndex <- obj@ndnIndex
  ndnLength <- obj@ndnLength
  
  trimmed_clone_sequence <- clone_sequence
  seq_lens <- nchar(clone_sequence)
  if(seq_edge_trim > 0) trimmed_clone_sequence <- substr(trimmed_clone_sequence, seq_edge_trim+1, seq_lens-2) # assumes max sequence < 1e9
  trimmed_clone_sequence <- gsub("N", "", trimmed_clone_sequence)
  n_trimmed_left_bases <- nchar(clone_sequence) - nchar(trimmed_clone_sequence) - seq_edge_trim
  
  # TO DO: go through logic and simplify
  {
    cdr3Index <- ifelse(n_trimmed_left_bases > 0, cdr3Index - n_trimmed_left_bases, cdr3Index)
    stopifnot(all(cdr3Index > 0)) #where did all the cdr3 go?
    
    ndnIndex <- ifelse(n_trimmed_left_bases > 0 & ndnLength > 0, ndnIndex - n_trimmed_left_bases, ndnIndex)
    
    select <- n_trimmed_left_bases > 0 & ndnLength > 0 & ndnIndex < 0
    ndnLength <- ifelse(select, ndnLength + ndnIndex, ndnLength)
    ndnIndex <- ifelse(select, 0, ndnIndex)
  }
  
  stopifnot(all(cdr3Index > 0 & cdr3Index <= seq_lens)) #bad cdr3 index
  stopifnot(all(ndnIndex >= 0 & (ndnIndex + ndnLength <= seq_lens + seq_edge_trim))) #bad ndnIndex
  
  obj@trimmed_sequence <- trimmed_clone_sequence
  obj@cdr3Index=cdr3Index
  obj@ndnLength=ndnLength
  obj@ndnIndex=ndnIndex
  return(obj)
})

#' splitSequenceIntoSubstrings
setMethod("splitSequenceIntoSubstrings", signature("EOS_Report"), function(obj) {
  
  clone_sequence <- obj@trimmed_sequence
  cdr3Index <- obj@cdr3Index
  ndnIndex <- obj@ndnIndex
  ndnLength <- obj@ndnLength
  
  left_seq <- substr(clone_sequence, 1, cdr3Index)
  right_seq <- substr(clone_sequence, cdr3Index + 1, 1e9)
  
  select1 <- ndnLength > 0 & ndnIndex < cdr3Index
  select2 <- ndnLength > 0 & ndnIndex +ndnLength > cdr3Index
  left_ndn_start <- ifelse(select1, ndnIndex, -1)
  left_ndn_end <- ifelse(select1, pmin(ndnIndex+ndnLength-1, cdr3Index-1), -1)
  right_ndn_start <- ifelse(select2, pmax(cdr3Index, ndnIndex), -1)
  right_ndn_end <- ifelse(select2, ndnIndex+ndnLength-1, -1)
  
  obj@left_sequence <- left_seq
  obj@right_sequence <- right_seq
  obj@left_ndn_start <- left_ndn_start
  obj@left_ndn_end <- left_ndn_end
  obj@right_ndn_start <- right_ndn_start
  obj@right_ndn_end <- right_ndn_end
  return(obj)
})

UpdateKeywordTree <- function(clone_sequence, locus_vec, ndn_start_vec, ndn_end_vec, is_left_seq_vec) {
  locus_index_vec <- as.integer(ave(locus_vec, locus_vec, FUN=seq_along))
  is_ndn_vec <- ifelse(ndn_end_vec > ndn_start_vec, 1, 0)
  
  # TO DO: optimize this loop
  for(ci in 1:length(clone_sequence)) {
    is_left_seq <- is_left_seq_vec[ci]
    is_ndn <- is_ndn_vec[ci]
    if(is_left_seq) {
      clone_sequence[[ci]] <- stringi::stri_reverse(clone_sequence[[ci]])
      if(is_ndn) {
        ndn_start_vec_temp <- ndn_start_vec[ci]
        ndn_start_vec[ci] <- nchar(clone_sequence[[ci]]) - ndn_end_vec[ci] - 1
        ndn_end_vec[ci] <- nchar(clone_sequence[[ci]]) - ndn_start_vec_temp - 1
      }
    }
    if(!is_ndn) {
      ndn_end_vec[ci] <- -1
      ndn_start_vec[ci] <- -1
    }
  }
  
  sapply(unique(locus_vec), function(locus) {
    lx <- locus_vec == locus
    C_updatekeyWordTree(clone_sequence[lx], as.integer(ndn_start_vec[lx]), as.integer(ndn_end_vec[lx]), locus_index_vec[lx], vector(mode="raw"))
  }, simplify=F)
}

#' createKeywordTrees
setMethod("createKeywordTrees", signature("EOS_Report"), function(obj) {
  obj@left_hash <- UpdateKeywordTree(obj@left_sequence, obj@locus, obj@left_ndn_start, obj@left_ndn_end, rep(T, length(obj@left_sequence)))
  obj@right_hash <- UpdateKeywordTree(obj@right_sequence, obj@locus, obj@right_ndn_start, obj@right_ndn_end, rep(F, length(obj@right_sequence)))
  return(obj)
})

#' processReport
#' Calls functions in the following order: getLocus, determineNdnValues, determineCdr3Index, trimCloneSequence, splitSequenceIntoSubstrings, createKeywordTrees
#' @param seq_edge_trim Value passed to trimCloneSequence
#' @return The processed EOS report object
processReport <- function(obj, seq_edge_trim=2L) {
  # obj %>% getLocus %>% 
  #   determineNdnValues %>% 
  #   determineCdr3Index %>% 
  #   trimCloneSequence(seq_edge_trim=seq_edge_trim) %>% 
  #   splitSequenceIntoSubstrings %>% 
  #   createKeywordTrees
  obj <- getLocus(obj)
  obj <- determineNdnValues(obj)
  obj <- determineCdr3Index(obj)
  obj <- trimCloneSequence(obj, seq_edge_trim=seq_edge_trim)
  obj <- splitSequenceIntoSubstrings(obj)
  obj <- createKeywordTrees(obj)
  return(obj)
}


#' Match sequences between objects
#' Match sequences in one object to sequences in another object using keyword hash trees
#' @param obj_query The query EOS report object
#' @param obj_target The target EOS report object (i.e., the one with the trees)
#' @param no_shm_allowance Should mutations be allowed?  
#' @param max_adj_mut_igkl If the EOS report query doesn't have a pre-specified max adjusted mutation value for a sequence, use this default for IGKL sequences
#' @param max_adj_mut_igh If the EOS report query doesn't have a pre-specified max adjusted mutation value for a sequence, use this default for IGH/IGH_D sequences
#' @param ndn_cost Modified cost for an ndn base pair mismatch
#' @param include_report_data If TRUE, append all EOS fields to the data.frame of matched sequences
#' @return A data.frame of matching sequences, locus, superlocus, total mutations and optionally the EOS report fields
matchQuerySequenceWithTargetSequence <- function(obj_query, obj_target, no_shm_allowance=T,
                                                 max_adj_mut_igkl=0, max_adj_mut_igh=2,
                                                 ndn_cost = 1.05, include_report_data=F) {
  # obj_target is the object with the hash trees to search
  
  max_adj_mut <- obj_query@data$maxAdjustedMutations
  if(no_shm_allowance) {
    max_adj_mut[] <- 0
  } else {
    max_adj_mut <- dplyr::case_when(is.na(max_adj_mut) ~ {
      ifelse(grepl("IGH",obj_query@locus), max_adj_mut_igh, max_adj_mut_igkl)
    },
    T ~ max_adj_mut)
  }
  
  # to do: this is kind of messy, should simplify somehow later
  loci <- intersect(unique(obj_query@locus), unique(obj_target@locus))
  data.table::rbindlist(lapply(loci, function(locus) {
    lx <- obj_query@locus == locus
    mleft <- C_findMatchingIdx(stringi::stri_reverse(obj_query@left_sequence[lx]), max_adj_mut[lx], ndn_cost, obj_target@left_hash[[locus]])
    mright <- C_findMatchingIdx(obj_query@right_sequence[lx], max_adj_mut[lx], ndn_cost, obj_target@right_hash[[locus]])
    data.table::rbindlist(lapply(1:length(mleft), function(i) {
      idx <- intersect(names(mleft[[i]]), names(mright[[i]])) # idx are the indices of the target (i.e., the ones with trees)
      msum <- sapply(idx, function(j) unname(mleft[[i]][j] + mright[[i]][j]))
      select <- idx[msum <= max_adj_mut[i]]
      if(length(select) == 0) return(NULL);
      superlocus <- ifelse(grepl("IGH", locus), "IGH","IGKL")
      querySeq <- obj_query@data$nucleotide[obj_query@locus == locus][i]
      targetSeq <- obj_target@data$nucleotide[obj_target@locus == locus][as.numeric(select)]
      if(!include_report_data) {
        data.frame(querySequence=querySeq, targetSequence=targetSeq, locus=locus, superlocus=superlocus, adjustedMutations = msum[select], stringsAsFactors = F)
      } else {
        ms <- match(targetSeq, obj_target@data$nucleotide)
        cbind(querySequence=querySeq, targetSequence=targetSeq, locus=locus, superlocus=superlocus, adjustedMutations=msum[select], obj_target@data[ms,], stringsAsFactors=F)
      }
    }))
  }))
}

#' Find consensus sequences
#' Find consensus sequences within EOS report(s)
#' @param ... Any number of EOS report objects
#' @param report_list A list of EOS report objects
#' @param no_shm_allowance Should mutations be allowed?  
#' @param max_adj_mut_igkl If the EOS report query doesn't have a pre-specified max adjusted mutation value for a sequence, use this default for IGKL sequences
#' @param max_adj_mut_igh If the EOS report query doesn't have a pre-specified max adjusted mutation value for a sequence, use this default for IGH/IGH_D sequences
#' @param ndn_cost Modified cost for an ndn base pair mismatch
#' @return A list with S3 class "Consensus_List".  Each element in the list is a vector of matching sequences and corresponding locus in the "locus" attribute field.  The names of the list correspond to the dominant sequence of the consensus cluster (i.e., the shortest).  
findConsensusSequences <- function(..., report_list=NULL, no_shm_allowance=T, max_adj_mut_igkl=0, max_adj_mut_igh=2, ndn_cost = 1.05) {
  report_list <- c(report_list, list(...))
  sequence_matches <- lapply(seq_along(report_list), function(i) {
    matchQuerySequenceWithTargetSequence(report_list[[i]], report_list[[i]],
                                         no_shm_allowance=no_shm_allowance, max_adj_mut_igkl=max_adj_mut_igkl, 
                                         max_adj_mut_igh=max_adj_mut_igh, ndn_cost=ndn_cost)
    })
  
  # Use igraph to find consensus sequence by looking for clusters
  # This approach is more robust to edge cases that will probably never happen (i.e., if not all sequences match to all other sequences in a consensus)
  sequence_matches <- data.table::rbindlist(sequence_matches)
  sequence_matches <- sequence_matches[,c("querySequence", "targetSequence", "locus")]
  sequence_matches <- unique(sequence_matches)
  consensi <- lapply(unique(sequence_matches$locus), function(loc) {
    seq_mat_loq <- sequence_matches[sequence_matches$locus == loc, ]
    seq_mat_loq <- as.matrix(seq_mat_loq[ ,c("querySequence", "targetSequence")])
    g <- igraph::graph_from_edgelist(seq_mat_loq)
    cl <- igraph::clusters(g)$membership
    consensus <- lapply(unique(cl), function(i) {
      names(cl[cl == i])
    })
    names(consensus) <- sapply(consensus, function(cons) {
      nc <- nchar(cons)
      cons[which.min(nc)]
    })
    for(i in seq_along(consensus)) {
      attr(consensus[[i]], "locus") <- loc
    }
    return(consensus)
  })
  consensi <- do.call(c, consensi)
  class(consensi) <- "Consensus_List"
  return(consensi)
}

#' Match Sequences to a Consensus
#' Find exact matching sequences in an EOS report object that match to any sequences in a "Consensus_List"
#' @param obj_query The EOS report query
#' @param consensus The Consensus List
#' @param include_report_data If TRUE, append all EOS fields to the data.frame of matched sequences
#' @return A data.frame of matching sequence, consensus sequence locus and optionally the EOS report fields
matchReportToConsensus <- function(obj_query, consensus, include_report_data=F) {
  data.table::rbindlist(lapply(seq_along(consensus), function(i) {
    loc <- attr(consensus[[i]], "locus")
    cons_seq <- names(consensus)[i]
    select <- obj_query@locus == loc & obj_query@data$nucleotide %in% consensus[[i]]
    if(sum(select) > 0) {
      if(include_report_data) {
        cbind(Consensus_Sequence=cons_seq, Locus=loc, obj_query@data[select])
      } else {
        data.frame(Consensus_Sequence=cons_seq, Locus=loc, nucleotide=obj_query@data$nucleotide[select], stringsAsFactors=F)
      }
    } else {
      return(NULL)
    }
  }))
}

#' Extract metadata from report
#' Extract metadata from the header of a report file
#' @param file_name The name of the EOS report file.  Can be gzipped; will be determined through file extension.  
#' @param directory Optional directory location of file
#' @param property The name of the property to extract
#' @return The metadata value.  Note: the return will be a string value.  
getReportMetadata <- function(file_name, directory="", property = "estTotalNucleatedCells") {
  file_name <- paste0(directory, file_name)
  if(grepl(".gz$", file_name)) {
    gz <- gzfile(file_name)
    lines <- readLines(gz,100)
    close(gz)
  } else {
    lines <- readLines(con = file_name, 100)
  }
  ret <- grep("^\\#", lines, value=T)
  ret <- gsub("#", "", ret) 
  ret <- strsplit(ret, "=")
  ret <- data.frame(do.call(rbind,ret), stringsAsFactors=F)
  colnames(ret) <- c("property", "value")
  ret <- ret$value[ret$property == property]
  stopifnot(length(ret) == 1) # more than one property
  return(ret)
}

#' Read in data and select by tag
#' Helper function for reading in data.  Given a report file name, reads in the data and selects rows based on the sequenceTags column.  
#' @param file_names A vector of file names
#' @param directory Optional directory location of file
#' @param tag The tag to select (by exact string matching the field); if `all` return all data.
#' @return A data.table of the report.  If more than one file was named, the reports are combined.  
getReportSeqsByTag <- function(file_names, directory="", tag="Dx") {
  data.table::rbindlist(lapply(file_names, function(fi) {
    df <- data.table::fread(sprintf("zcat < '%s%s'", directory, fi), sep="\t", colClasses = report_col_classes)
    if(tag == "all") {
      return(df)
    } else {
      return( df[df$sequenceTags == tag,] )
    }
  }))
}
