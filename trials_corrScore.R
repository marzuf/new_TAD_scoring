SSHFS <- TRUE
setDir <- ifelse(SSHFS, "/media/electron", "")

suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, 40))

source("utils_fct.R")

dataFold <- file.path(setDir, "/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(dataFold))

hicds <- "GSE105318_DLD1_40kb"
exprds <- "TCGAcoad_msi_mss"
stopifnot(dir.exists(file.path(dataFold, hicds)))

pipOutFold <- file.path(dataFold, "PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
stopifnot(dir.exists(file.path(pipOutFold)))

outFold <- file.path("TRIALS_CORRSCORE")
dir.create(outFold, recursive = TRUE)

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script_name <- "4v2_runScoreTADCorr"
script4_name <- "4_runMeanTADCorr"
script4v2_name <- "4v2_runConcordTADCorr"
script4v3_name <- "4v3_runConcordTADCorr"
script4vAll_name <- "4vAll_runConcordTADCorr"


diffTADcol <- "blue"
inTADcol <- "red"
plotCex <- 1.2

norm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_qqnorm_rnaseqDT.Rdata"))))
initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

meanTADcorr <- eval(parse(text = load(file.path(pipOutFold,  script4_name, "all_meanCorr_TAD.Rdata"))))
meanTADconcord_v2 <- eval(parse(text = load(file.path(pipOutFold,  script4v2_name, "all_concordCorr_TAD.Rdata"))))
meanTADconcord_v3 <- eval(parse(text = load(file.path(pipOutFold,  script4v3_name, "all_concordCorr_TAD.Rdata"))))
meanTADconcord_vAll <- eval(parse(text = load(file.path(pipOutFold,  script4vAll_name, "all_concordCorr_TAD.Rdata"))))

ratioMagnitude <- unlist(lapply(meanTADconcord_vAll, function(x)x[["ratioMagnitude"]]))
ratioNbr <- unlist(lapply(meanTADconcord_vAll, function(x)x[["ratioNbr"]]))

densplot(x=ratioMagnitude, y = ratioNbr,
         xlab="ratioMagnitude", ylab="ratioNbr")
abline(v=0.5, h=0.5, lty=2)


f1 <- function(x,y) x*y
cf_func(f1, xlim = c(0,1),ylim = c(0, 1))


result <- function(x, y,d){
  ydata=numeric(length(seq(y[1],y[2],d)))
  mat=matrix(NA,length(seq(y[1],y[2],d)),length(seq(x[1],x[2],d)))
  yy=seq(y[1],y[2],d)
  xx=seq(x[1],x[2],d)
  for(i in 1:length(xx)){
    for(n in 1:length(yy)){
      ydata[n]=xx[i]*yy[n]
    }
    mat[,i]=ydata
  }
  return(mat)}


rotate <- function(x) t(apply(x, 2, rev))


mydata=result(x=c(0,1),y=c(0,1),d=10)
image(rotate(mydata))
par(new=TRUE)
contour(mydata)





stopifnot(setequal(names(meanTADcorr), names(meanTADconcord_vAll)))

stopifnot(setequal(names(meanTADcorr), names(meanTADconcord_v2)))
stopifnot(setequal(names(meanTADcorr), names(meanTADconcord_v3)))
commonDS <- intersect(names(meanTADcorr), names(meanTADconcord_v2))

plot(x=meanTADcorr[commonDS], y=meanTADconcord_v2[commonDS])
plot(x=meanTADcorr[commonDS], y=meanTADconcord_v3[commonDS])
plot(x=meanTADconcord_v2[commonDS], y=meanTADconcord_v3[commonDS])

plot_multiDens(
  list(
    meanTADcorr=meanTADcorr,
    meanTADconcord_v2=meanTADconcord_v2,
    meanTADconcord_v3 =meanTADconcord_v3
    
  )
)


norm_rnaseqDT <- norm_rnaseqDT[names(geneList),]     
stopifnot(all(rownames(norm_rnaseqDT) == names(geneList)))
stopifnot(!any(duplicated(names(geneList))))

gene2tadDT_file <- file.path(dataFold, hicds, "genes2tad", "all_genes_positions.txt")
stopifnot(file.exists(gene2tadDT_file))
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(geneList),]

corrMethodCoexpr <- "pearson"

genecoexpr_file <- file.path(dataFold,"CREATE_COEXPR_SORTNODUP", hicds, paste0(exprds, "_", corrMethodCoexpr), "coexprDT.Rdata")
stopifnot(file.exists(genecoexpr_file))
cat("... load geneCoexpr data\n")
geneCoexpr_DT <- eval(parse(text = load(genecoexpr_file)))

chr20_geneCoexpr_DT <- geneCoexpr_DT[geneCoexpr_DT$chromo == "chr20",]
save(chr20_geneCoexpr_DT, file = "chr20_geneCoexpr_DT.Rdata")
geneCoexpr_DT <- eval(parse(text = load("chr20_geneCoexpr_DT.Rdata")))

genedist_file <- file.path(dataFold, "CREATE_DIST_SORTNODUP", hicds,"all_dist_pairs.Rdata")
stopifnot(file.exists(genedist_file))
cat("... load geneDist data\n")
geneDist_DT <- eval(parse(text = load(genedist_file)))

chr20_geneDist_DT <- geneDist_DT[geneDist_DT$chromo == "chr20",]
save(chr20_geneDist_DT, file = "chr20_geneDist_DT.Rdata")
geneDist_DT <- eval(parse(text = load("chr20_geneDist_DT.Rdata")))

genesametad_file <- file.path(dataFold, "CREATE_SAME_TAD_SORTNODUP", hicds,"all_TAD_pairs.Rdata")
stopifnot(file.exists(genesametad_file))
cat("... load geneSameTAD data\n")
geneSameTAD_DT <- eval(parse(text = load(genesametad_file)))

chr20_geneSameTAD_DT <- geneSameTAD_DT[geneSameTAD_DT$chromo == "chr20",]
save(chr20_geneSameTAD_DT, file = "chr20_geneSameTAD_DT.Rdata")
geneSameTAD_DT <- eval(parse(text = load("chr20_geneSameTAD_DT.Rdata")))

### take only the filtered data according to initial settings
pipeline_regionList <- eval(parse(text = load(file.path(pipOutFold, script0_name, "pipeline_regionList.Rdata"))))

gene2tadDT <- gene2tadDT[gene2tadDT$region %in% pipeline_regionList,]

all_regions <- unique(as.character(gene2tadDT$region))
stopifnot(grepl("_TAD", all_regions))

reg <- all_regions[grepl("chr20_", all_regions)][1]

reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
rowsToKeep <- which(geneList %in% reg_genes)
subData <- as.data.frame(norm_rnaseqDT[rowsToKeep,,drop=F])
stopifnot(rownames(subData) == names(geneList[rowsToKeep]))

#*************************************************************************************************************
# prepare inTAD and diffTAD data

head(geneCoexpr_DT)
head(geneDist_DT)
head(geneSameTAD_DT)

### inTAD
tad_nGenes <- length(reg_genes)

tad_coexprDT <- geneCoexpr_DT[geneCoexpr_DT$gene1 %in% reg_genes &
                                geneCoexpr_DT$gene2 %in% reg_genes, ]
head(tad_coexprDT)

tad_distDT <- geneDist_DT[geneDist_DT$gene1 %in% reg_genes &
                            geneDist_DT$gene2 %in% reg_genes, ]
head(tad_distDT)

tad_sameTADdt <- geneSameTAD_DT[geneSameTAD_DT$gene1 %in% reg_genes &
                                  geneSameTAD_DT$gene2 %in% reg_genes, ]
head(tad_sameTADdt)
tad_chromo <- gsub("(^chr.+?)_TAD.+$", "\\1", unique(tad_sameTADdt$region))
stopifnot(length(tad_chromo) == 1)

tad_coexprDistDT <- merge(tad_coexprDT, tad_distDT, 
                          by=c("gene1","gene2"), 
                          all.x=TRUE, all.y=TRUE)
stopifnot(!is.na(tad_coexprDistDT))
stopifnot(tad_nGenes * (tad_nGenes-1)/2 == nrow(tad_coexprDistDT) )
head(tad_coexprDistDT)
nrow(tad_coexprDistDT)

### diffTAD
tad_dist <- max(tad_coexprDistDT$dist)
stopifnot(is.numeric(tad_dist))

tadRange_distDT <- geneDist_DT[geneDist_DT$dist <= tad_dist,]
range_genes <- unique(c(tadRange_distDT$gene2, tadRange_distDT$gene1))

tadRange_coexprDT <- geneCoexpr_DT[geneCoexpr_DT$gene1 %in% range_genes &
                                     geneCoexpr_DT$gene2 %in% range_genes,]

tadRange_coexprDistDT <- merge(tadRange_coexprDT, tadRange_distDT, 
                               by=c("gene1","gene2"), 
                               all.x=FALSE, all.y=FALSE)
nrow(tadRange_coexprDistDT)
stopifnot(nrow(tadRange_coexprDistDT) > 0)
tadRange_coexprDistDT <- tadRange_coexprDistDT[! tadRange_coexprDistDT$gene1 %in% tad_sameTADdt$gene1 &
                                                 ! tadRange_coexprDistDT$gene2 %in% tad_sameTADdt$gene2,
                                               ]
nrow(tadRange_coexprDistDT)
stopifnot(tadRange_coexprDistDT$chromo == tad_chromo)
stopifnot( nrow(tadRange_coexprDistDT) > 0 )

#*************************************************************************************************************
### IDEA 1:
# pour chaque TAD, je peux fitter une coexpr ~ dist
# et faire de même sur le même distance range avec des paires de gènes qui sont dans des TADs différents
# ensuite le score qqh comme AUC
#*************************************************************************************************************

my_xlab <- "dist (bp) [TAD range]"
my_ylab <- paste0("coexpr (", corrMethodCoexpr, ")")


#### inTAD

# fit models
inTAD_lm <- lm(coexpr ~ dist, data = tad_coexprDistDT)
inTAD_mod <- loess(coexpr ~ dist, data = tad_coexprDistDT)
smooth_vals_inTAD <- predict(inTAD_mod, sort(tad_coexprDistDT$dist))
auc_inTAD_obsDist <- auc(x = sort(tad_coexprDistDT$dist), y = smooth_vals_inTAD)

# plot
my_tit <- paste0(reg)
plot(coexpr ~ dist, 
     data=tad_coexprDistDT,
     pch=16,
     cex=0.7,
     cex.lab=plotCex,
     cex.axis = plotCex,
     xlab = my_xlab,
     ylab = my_ylab, 
     main = my_tit)
abline(inTAD_lm, col = inTADcol)
lines( x = sort(tad_coexprDistDT$dist), y = smooth_vals_inTAD, col = inTADcol)

#### diffTAD
# fit models
smooth_vals_diffTAD <- predict(diffTAD_mod, sort(tadRange_coexprDistDT$dist))
auc_diffTAD_obsDist <- auc(x = sort(tadRange_coexprDistDT$dist), y = smooth_vals_diffTAD)
diffTAD_lm <- lm(coexpr ~ dist, data = tadRange_coexprDistDT)
diffTAD_mod <- loess(coexpr ~ dist, data = tadRange_coexprDistDT)

# plot
my_tit <- paste0(reg, " - diffTAD gene pairs")
plot(coexpr ~ dist, 
     data=tadRange_coexprDistDT,
     pch=16,
     cex=0.7,
     cex.lab=plotCex,
     cex.axis = plotCex,
     xlab = my_xlab,
     ylab = my_ylab, 
     main = my_tit)
abline(diffTAD_lm, col = diffTADcol)
lines( x = sort(tadRange_coexprDistDT$dist), y = smooth_vals_diffTAD, col = diffTADcol)

auc_inTAD_diffTAD_ratio <- auc_inTAD_obsDist/auc_diffTAD_obsDist

inTAD_diffTAD_lmSlopeRatio <- as.numeric(coef(inTAD_lm)["dist"])/as.numeric(coef(diffTAD_lm)["dist"])
stopifnot(!is.na(inTAD_diffTAD_lmSlopeRatio))

my_tit <- paste0(reg, " - inTAD/diffTAD")
plot(NULL,
     xlim = range(tad_coexprDistDT$dist), 
     ylim = range(c(smooth_vals_diffTAD, smooth_vals_inTAD)),
     xlab=my_xlab, 
     ylab=my_ylab,
     main=my_tit)
mtext(text = paste0("(AUC ratio=", round(auc_inTAD_diffTAD_ratio,2), "); slope ratio=",  round(inTAD_diffTAD_lmSlopeRatio,2), ")"), side = 3)
lines( x = sort(tad_coexprDistDT$dist), y = smooth_vals_inTAD, col = inTADcol, lty=2)
lines( x = sort(tadRange_coexprDistDT$dist), y = smooth_vals_diffTAD, col = diffTADcol, lty=2)
abline(inTAD_lm, col = inTADcol)
abline(diffTAD_lm, col = diffTADcol)

legend("topright",
       legend = c(reg, "diffTAD"),
       col = c(inTADcol, diffTADcol),
       lty=1, bty="n")

# score is then:
auc_inTAD_diffTAD_ratio
# or
inTAD_diffTAD_lmSlopeRatio  # or 1-slopeRatio to have higher score better

#*************************************************************************************************************
### IDEA 2:
# pour chaque TAD, pour chaque paire de gènes, trouver une paire qui est dans une distance similaire (+- un radius)
# et calculer qqh comme un log2FC et prendre la moyenne (pour chaque paire, puis pour le TAD)
#*************************************************************************************************************
head(tadRange_coexprDistDT)
head(tad_coexprDistDT)

distRad <- 1000

i = 1
inTAD_diffTAD_pairCoexpr <- foreach(i = seq_len(nrow(tad_coexprDistDT))) %dopar% {
  i_dist <- tad_coexprDistDT$dist[i]
  diffGenes_coexprDistDT <- tadRange_coexprDistDT[
    abs(tadRange_coexprDistDT$dist - i_dist) <= distRad,
    ]
  obs_coexpr <- tad_coexprDistDT$coexpr[i]
  diffAvg_coexpr <- mean(diffGenes_coexprDistDT$coexpr)
  
  list(obs_coexpr=obs_coexpr, diffAvg_coexpr=diffAvg_coexpr)
  #obs_coexpr/diffAvg_coexpr
}
inTAD_diffTAD_pairCoexprRatio <- unlist(lapply(inTAD_diffTAD_pairCoexpr, function(x)
                                        x[["obs_coexpr"]]/x[["diffAvg_coexpr"]]))

# TAD score is then:
mean(inTAD_diffTAD_pairCoexprRatio, na.rm=TRUE)
# OR
mean(inTAD_diffTAD_pairCoexprRatio > 1)

#### z-score transformation ??????


#*************************************************************************************************************
### IDEA 4:
# soit g-1, g-2, g-3, ... les gènes à gauche du TAD, et g+1, g+2, g+3 les gènes à droite du TAD
# et g1, g2, g3, ... les gènes dans le TAD
# calculer les meanCorr pour les sets
# m1={g1,g2,g3,..gn} m2={g2,g3,..gn, g+1} mk={gi, ..., g+j} tant que g+j est dans le même intervalle de TAD size
# ensuite faire qqh comme un cumdiff
#*************************************************************************************************************
chromo_gene2tadDT <- gene2tadDT[as.character(gene2tadDT$chromo) == as.character(tad_chromo),]
stopifnot(nrow(chromo_gene2tadDT) > 0)
stopifnot(is.numeric(chromo_gene2tadDT$start))
stopifnot(is.numeric(chromo_gene2tadDT$end))
chromo_gene2tadDT <- chromo_gene2tadDT[order(chromo_gene2tadDT$start, chromo_gene2tadDT$end),]
head(chromo_gene2tadDT)

reg_genes <- chromo_gene2tadDT$entrezID[chromo_gene2tadDT$region == reg]
stopifnot(reg_genes %in% chromo_gene2tadDT$entrezID)

idx_tadStart <- min(which(chromo_gene2tadDT$entrezID %in% reg_genes))
idx_tadEnd <- max(which(chromo_gene2tadDT$entrezID %in% reg_genes))

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
    
    gene_set <- c(reg_genes[-c((length(reg_genes)-i+1):length(reg_genes))], 
                  left_genes[length(left_genes)-i+1])
    
    stopifnot(length(gene_set) == tad_nGenes)
    # stopifnot(max(chromo_gene2tadDT$end[chromo_gene2tadDT$entrezID %in% gene_set]) - 
    #   min(chromo_gene2tadDT$start[chromo_gene2tadDT$entrezID %in% gene_set]) + 1 <= tad_size)
    # -> do not impose this constraint (too stringent ?)
    stopifnot(min(chromo_gene2tadDT$start[chromo_gene2tadDT$entrezID %in% gene_set]) -
                min(chromo_gene2tadDT$start[chromo_gene2tadDT$entrezID %in% reg_genes]) + 1 <= tad_size)
    stopifnot(gene_set %in% geneList)
    rowsToKeep <- which(geneList %in% gene_set)
    subData <- as.data.frame(t(norm_rnaseqDT[rowsToKeep,,drop=F]))
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
    gene_set <- c(reg_genes[-c(1:i)], right_genes[c(1:i)])
    stopifnot(length(gene_set) == tad_nGenes)
    # stopifnot(max(chromo_gene2tadDT$end[chromo_gene2tadDT$entrezID %in% gene_set]) - 
    #   min(chromo_gene2tadDT$start[chromo_gene2tadDT$entrezID %in% gene_set]) + 1 <= tad_size)
    # -> do not impose this constraint (too stringent ?)
    stopifnot(max(chromo_gene2tadDT$end[chromo_gene2tadDT$entrezID %in% gene_set]) -
      max(chromo_gene2tadDT$end[chromo_gene2tadDT$entrezID %in% reg_genes]) + 1 <= tad_size)
    stopifnot(gene_set %in% geneList)
    rowsToKeep <- which(geneList %in% gene_set)
    subData <- as.data.frame(t(norm_rnaseqDT[rowsToKeep,,drop=F]))
    stopifnot(colnames(subData) == names(geneList[rowsToKeep]))
    stopifnot(ncol(subData) == length(gene_set))
    corrMatrix_all <- cor(subData, method = corrMethod)
    stopifnot(nrow(corrMatrix_all) == length(gene_set))
    stopifnot(ncol(corrMatrix_all) == length(gene_set))
    mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiag)], na.rm=T)
  }
}


rowsToKeep <- which(geneList %in% reg_genes)
subData <- as.data.frame(t(norm_rnaseqDT[rowsToKeep,,drop=F]))
stopifnot(colnames(subData) == names(geneList[rowsToKeep]))
stopifnot(ncol(subData) == length(reg_genes))
corrMatrix_all <- cor(subData, method = corrMethod)
stopifnot(nrow(corrMatrix_all) == length(reg_genes))
stopifnot(ncol(corrMatrix_all) == length(reg_genes))
tad_corr <- mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiag)], na.rm=T)

all_corrs <- na.omit(c(tad_corr, left_window_corrs, right_window_corrs))
stopifnot(length(all_corrs) > 1)

#TAD score is then 
as.numeric(scale(all_corrs, center = TRUE, scale=TRUE))


#*************************************************************************************************************
#*************************************************************************************************************
#*************************************************************************************************************
#*************************************************************************************************************





 <- foreach(reg=all_regions, .combine='c') %dopar% {
  # cat(paste0("... doing region: ", reg, "/", length(all_regions), "\n"))
  reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
  rowsToKeep <- which(geneList %in% reg_genes)
  subData <- as.data.frame(t(norm_rnaseqDT[rowsToKeep,,drop=F]))
  stopifnot(colnames(subData) == names(geneList[rowsToKeep]))
  # columns => the genes
  # rows => the samples
  ######################################################## BECAUSE I COULD HAVE DUPLICATED GENE ENTREZ ID FOR DIFFERENT NAMES !!!
  ### UPDATE: this should not happen in the latest version !!!
  # stopifnot(ncol(subData) == length(reg_genes))
  # stopifnot(ncol(subData) == length(geneList[which(geneList %in% reg_genes)]))
  stopifnot(ncol(subData) == length(reg_genes))
  #### CORRELATION
  corrMatrix_all <- cor(subData, method = corrMethod)
  # should be correlation of the genes
  ######################################################## BECAUSE I COULD HAVE DUPLICATED GENE ENTREZ ID FOR DIFFERENT NAMES !!!
  ### UPDATE: this should not happen in the latest version !!!
  # stopifnot(nrow(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
  # stopifnot(ncol(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
  stopifnot(nrow(corrMatrix_all) == length(reg_genes))
  stopifnot(ncol(corrMatrix_all) == length(reg_genes))
  mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiag)], na.rm=T)
}

cat(paste0("... end intra TAD correlation\n"))

names(all_meanCorr_TAD) <- all_regions
stopifnot(length(all_meanCorr_TAD) == length(all_regions))

save(all_meanCorr_TAD, file= paste0(curr_outFold, "/all_meanCorr_TAD.Rdata"))
cat(paste0("... written: ", curr_outFold, "/all_meanCorr_TAD.Rdata", "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))

