withDiag <- FALSE
corrMethod <- "pearson"
require(foreach)
require(flux)
# Xreg with tad
# Xreg_genes with tad_genes
# Xgene2tadDT with gt2_DT
# Xnorm_rnaseqDT with expr_DT
# tad_coexprDistDT with tad_coexprDist_DT

#*************************************************************************************************************
### IDEA 1:
# pour chaque TAD, je peux fitter une coexpr ~ dist
# et faire de même sur le même distance range avec des paires de gènes qui sont dans des TADs différents
# ensuite le score qqh comme AUC
#*************************************************************************************************************



tad_score_trial1 <- function(
  tad_coexprDist_DT,
  tadRange_coexprDist_DT,
  distRad=1000
) {
  
  #### inTAD
  # fit models
  inTAD_lm <- lm(coexpr ~ dist, data = tad_coexprDistDT)
  inTAD_mod <- loess(coexpr ~ dist, data = tad_coexprDistDT)
  smooth_vals_inTAD <- predict(inTAD_mod, sort(tad_coexprDistDT$dist))
  auc_inTAD_obsDist <- auc(x = sort(tad_coexprDistDT$dist), y = smooth_vals_inTAD)
  slope_inTAD <- as.numeric(coef(inTAD_lm)["dist"])
  
  #### diffTAD
  # fit models
  smooth_vals_diffTAD <- predict(diffTAD_mod, sort(tadRange_coexprDistDT$dist))
  auc_diffTAD_obsDist <- auc(x = sort(tadRange_coexprDistDT$dist), y = smooth_vals_diffTAD)
  diffTAD_lm <- lm(coexpr ~ dist, data = tadRange_coexprDistDT)
  diffTAD_mod <- loess(coexpr ~ dist, data = tadRange_coexprDistDT)
  slope_diffTAD <- as.numeric(coef(diffTAD_lm)["dist"])
  ### ratios
  auc_inTAD_diffTAD_ratio <- auc_inTAD_obsDist/auc_diffTAD_obsDist
  inTAD_diffTAD_lmSlopeRatio <- slope_inTAD/slope_diffTAD
  stopifnot(!is.na(auc_inTAD_diffTAD_ratio))
  stopifnot(!is.na(inTAD_diffTAD_lmSlopeRatio))
  
  # TAD score is then:
  auc_inTAD_diffTAD_ratio
  # or
  inTAD_diffTAD_lmSlopeRatio  # or 1-slopeRatio to have higher score better
  
  return(list(auc_inTAD_diffTAD_ratio = auc_inTAD_diffTAD_ratio, 
              inTAD_diffTAD_lmSlopeRatio = inTAD_diffTAD_lmSlopeRatio))
  
  return(list(
    auc_inTAD_obsDist = auc_inTAD_obsDist,
    auc_diffTAD_obsDist = auc_diffTAD_obsDist,
    slope_inTAD = slope_inTAD,
    slope_diffTAD = slope_diffTAD
  ))
  
  
}





#*************************************************************************************************************
### IDEA 2:
# pour chaque TAD, pour chaque paire de gènes, trouver une paire qui est dans une distance similaire (+- un radius)
# et calculer qqh comme un log2FC et prendre la moyenne (pour chaque paire, puis pour le TAD)
#*************************************************************************************************************

tad_score_trial2 <- function(
  tad_coexprDist_DT,
  tadRange_coexprDist_DT,
  distRad=1000
){
  i = 1
  inTAD_diffTAD_pairCoexpr <- foreach(i = seq_len(nrow(tad_coexprDist_DT))) %do% {
    i_dist <- tad_coexprDist_DT$dist[i]
    diffGenes_coexprDistDT <- tadRange_coexprDist_DT[
      abs(tadRange_coexprDist_DT$dist - i_dist) <= distRad,
      ]
    obs_coexpr <- tad_coexprDist_DT$coexpr[i]
    diffAvg_coexpr <- mean(diffGenes_coexprDistDT$coexpr)
    
    list(obs_coexpr=obs_coexpr, diffAvg_coexpr=diffAvg_coexpr)
    #obs_coexpr/diffAvg_coexpr
  }
  # inTAD_diffTAD_pairCoexprRatio <- unlist(lapply(inTAD_diffTAD_pairCoexpr, function(x)
  #   x[["obs_coexpr"]]/x[["diffAvg_coexpr"]]))
  # # TAD score is then:
  # mean(inTAD_diffTAD_pairCoexprRatio, na.rm=TRUE)
  # # OR
  # mean(inTAD_diffTAD_pairCoexprRatio > 1)
  # #### z-score transformation ??????
  return(inTAD_diffTAD_pairCoexpr)
}
  
  
  
  





#*************************************************************************************************************
### IDEA 4:
# soit g-1, g-2, g-3, ... les gènes à gauche du TAD, et g+1, g+2, g+3 les gènes à droite du TAD
# et g1, g2, g3, ... les gènes dans le TAD
# calculer les meanCorr pour les sets
# m1={g1,g2,g3,..gn} m2={g2,g3,..gn, g+1} mk={gi, ..., g+j} tant que g+j est dans le même intervalle de TAD size
# ensuite faire qqh comme un cumdiff
#*************************************************************************************************************

tad_score_trial4 <- function(
  tad,
  gt2_DT,
  expr_DT
){
  
  tad_genes <- gt2_DT$entrezID[gt2_DT$region == tad]
  tad_nGenes <- length(tad_genes)
  
  
  tad_chromo <- gsub("(^chr.+?)_TAD.+$", "\\1", unique(tad_sameTADdt$region))
  
  chromo_gene2tadDT <- gt2_DT[as.character(gt2_DT$chromo) == as.character(tad_chromo),]
  stopifnot(nrow(chromo_gene2tadDT) > 0)
  stopifnot(is.numeric(chromo_gene2tadDT$start))
  stopifnot(is.numeric(chromo_gene2tadDT$end))
  chromo_gene2tadDT <- chromo_gene2tadDT[order(chromo_gene2tadDT$start, chromo_gene2tadDT$end),]
  head(chromo_gene2tadDT)
  
  tad_genes <- chromo_gene2tadDT$entrezID[chromo_gene2tadDT$region == tad]
  stopifnot(tad_genes %in% chromo_gene2tadDT$entrezID)
  
  idx_tadStart <- min(which(chromo_gene2tadDT$entrezID %in% tad_genes))
  idx_tadEnd <- max(which(chromo_gene2tadDT$entrezID %in% tad_genes))
  
  tad_startPos <- chromo_gene2tadDT$start[idx_tadStart]
  tad_endPos <- chromo_gene2tadDT$end[idx_tadEnd]
  
  tad_size <- tad_endPos - tad_startPos + 1
  stopifnot(tad_size > 0)
  
  ### CORR ON LEFT 
  # do not windowing exceed tad size [but gene set size can be larger than tad size]
  left_gene2tadDT <- chromo_gene2tadDT[(tad_startPos - chromo_gene2tadDT$start + 1) <= tad_size & 
                                         tad_startPos > chromo_gene2tadDT$start,
                                       ]
  # do not take more genes
  if(nrow(left_gene2tadDT) > (tad_nGenes-1))
    left_gene2tadDT <- left_gene2tadDT[1:(tad_nGenes-1),]
  
  if(nrow(left_gene2tadDT) == 0) {
    left_window_corrs <- NA
  } else{
    
    left_genes <- left_gene2tadDT$entrezID
    left_window_corrs <- foreach(i = seq_len(nrow(left_gene2tadDT)), .combine='c') %do% {
      
      gene_set <- c(tad_genes[-c((length(tad_genes)-i+1):length(tad_genes))], 
                    left_genes[length(left_genes)-i+1])
      
      stopifnot(length(gene_set) == tad_nGenes)
      # stopifnot(max(chromo_gene2tadDT$end[chromo_gene2tadDT$entrezID %in% gene_set]) - 
      #   min(chromo_gene2tadDT$start[chromo_gene2tadDT$entrezID %in% gene_set]) + 1 <= tad_size)
      # -> do not impose this constraint (too stringent ?)
      stopifnot(min(chromo_gene2tadDT$start[chromo_gene2tadDT$entrezID %in% gene_set]) -
                  min(chromo_gene2tadDT$start[chromo_gene2tadDT$entrezID %in% tad_genes]) + 1 <= tad_size)
      stopifnot(gene_set %in% geneList)
      rowsToKeep <- which(geneList %in% gene_set)
      subData <- as.data.frame(t(expr_DT[rowsToKeep,,drop=F]))
      stopifnot(colnames(subData) == names(geneList[rowsToKeep]))
      stopifnot(ncol(subData) == length(gene_set))
      corrMatrix_all <- cor(subData, method = corrMethod)
      stopifnot(nrow(corrMatrix_all) == length(gene_set))
      stopifnot(ncol(corrMatrix_all) == length(gene_set))
      mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiag)], na.rm=T)
    }
  }
  
  # CORR ON RIGHT
  # do not windowing exceed tad size [but gene set size can be larger than tad size]
  right_gene2tadDT <- chromo_gene2tadDT[(chromo_gene2tadDT$end - tad_endPos + 1) <= tad_size & 
                                          tad_endPos < chromo_gene2tadDT$end,
                                        ]
  # do not take more genes
  if(nrow(right_gene2tadDT) > (tad_nGenes-1))
    right_gene2tadDT <- right_gene2tadDT[1:(tad_nGenes-1),]
  
  if(nrow(right_gene2tadDT) == 0) {
    right_window_corrs <- NA
  } else{
    right_genes <- right_gene2tadDT$entrezID
    right_window_corrs <- foreach(i = seq_len(nrow(right_gene2tadDT)), .combine='c') %do% {
      gene_set <- c(tad_genes[-c(1:i)], right_genes[c(1:i)])
      stopifnot(length(gene_set) == tad_nGenes)
      # stopifnot(max(chromo_gene2tadDT$end[chromo_gene2tadDT$entrezID %in% gene_set]) - 
      #   min(chromo_gene2tadDT$start[chromo_gene2tadDT$entrezID %in% gene_set]) + 1 <= tad_size)
      # -> do not impose this constraint (too stringent ?)
      stopifnot(max(chromo_gene2tadDT$end[chromo_gene2tadDT$entrezID %in% gene_set]) -
                  max(chromo_gene2tadDT$end[chromo_gene2tadDT$entrezID %in% tad_genes]) + 1 <= tad_size)
      stopifnot(gene_set %in% geneList)
      rowsToKeep <- which(geneList %in% gene_set)
      subData <- as.data.frame(t(expr_DT[rowsToKeep,,drop=F]))
      stopifnot(colnames(subData) == names(geneList[rowsToKeep]))
      stopifnot(ncol(subData) == length(gene_set))
      corrMatrix_all <- cor(subData, method = corrMethod)
      stopifnot(nrow(corrMatrix_all) == length(gene_set))
      stopifnot(ncol(corrMatrix_all) == length(gene_set))
      mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiag)], na.rm=T)
    }
  }
  
  rowsToKeep <- which(geneList %in% tad_genes)
  subData <- as.data.frame(t(expr_DT[rowsToKeep,,drop=F]))
  stopifnot(colnames(subData) == names(geneList[rowsToKeep]))
  stopifnot(ncol(subData) == length(tad_genes))
  corrMatrix_all <- cor(subData, method = corrMethod)
  stopifnot(nrow(corrMatrix_all) == length(tad_genes))
  stopifnot(ncol(corrMatrix_all) == length(tad_genes))
  tad_corr <- mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiag)], na.rm=T)
  
  # all_corrs <- na.omit(c(tad_corr, left_window_corrs, right_window_corrs))
  # stopifnot(length(all_corrs) > 1)
  # #TAD score is then 
  # tad_score <- as.numeric(scale(all_corrs, center = TRUE, scale=TRUE))
  # return(tad_score)
  
  all_corrs <- list(tad_corr=tad_corr, 
                    left_window_corrs=left_window_corrs, 
                    right_window_corrs=right_window_corrs)
  return(all_corrs)
  
  
  
}
