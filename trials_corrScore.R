SSHFS <- TRUE
setDir <- ifelse(SSHFS, "~/media/electron", "")

suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, 40))

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

norm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_qqnorm_rnaseqDT.Rdata"))))
initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

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
### IDEA 1:
# pour chaque TAD, je peux fitter une coexpr ~ dist
# et faire de même sur le même distance range avec des paires de gènes qui sont dans des TADs différents
# ensuite le score qqh comme AUC
#*************************************************************************************************************

head(geneCoexpr_DT)
head(geneDist_DT)
head(geneSameTAD_DT)

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
plot(coexpr ~ dist, 
     data=tad_coexprDistDT)

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

tadRange_coexprDistDT <- tadRange_coexprDistDT[! tadRange_coexprDistDT$gene1 %in% tad_sameTADdt$gene1 &
                                                 ! tadRange_coexprDistDT$gene2 %in% tad_sameTADdt$gene2,
                                                 ]
nrow(tadRange_coexprDistDT)

stopifnot(tadRange_coexprDistDT$chromo == tad_chromo)

stopifnot( nrow(tadRange_coexprDistDT) > 0 )

plot(coexpr ~ dist, 
     data=tadRange_coexprDistDT)



diffTADcol ="blue"
inTADcol ="red"

diffTAD_mod <- loess(coexpr ~ dist, data = tadRange_coexprDistDT)
sameTAD_mod <- loess(coexpr ~ dist, data = tad_coexprDistDT)

smooth_vals_inTAD <- predict(sameTAD_mod, sort(tad_coexprDistDT$dist))
smooth_vals_diffTAD <- predict(diffTAD_mod, sort(tadRange_coexprDistDT$dist))

auc_diffTAD_obsDist <- auc(x = sort(tadRange_coexprDistDT$dist), y = smooth_vals_diffTAD)
auc_sameTAD_obsDist <- auc(x = sort(tad_coexprDistDT$dist), y = smooth_vals_inTAD)

lines( x = sort(tad_coexprDistDT$dist), y = smooth_vals_inTAD, col = inTADcol)
lines( x = sort(tadRange_coexprDistDT$dist), y = smooth_vals_diffTAD, col = diffTADcol)


my_xlab <- paste0("")
my_ylab <- paste0("")
my_tit <- paste0("")
plot(NULL,
     xlim = range(tad_coexprDistDT$dist), 
     ylim = range(c(smooth_vals_diffTAD, smooth_vals_inTAD)),
     xlab=my_xlab, 
     ylab=my_ylab,
     main=my_tit)
mtext(text = "observed distance values", side = 3)
lines( x = sort(diffTAD_mod$dist), y = smooth_vals_inTAD, col = inTADcol)
lines( x = sort(tadRange_coexprDistDT$dist), y = smooth_vals_diffTAD, col = diffTADcol)



require(flux)
auc_diffTAD_obsDist <- auc(x = sort(tadRange_coexprDistDT$dist), y = smooth_vals_diffTAD)
auc_diffTAD_obsDist
auc_inTAD_obsDist <- auc(x = sort(tad_coexprDistDT$dist), y = smooth_vals_inTAD)
auc_inTAD_obsDist




all_meanCorr_TAD <- foreach(reg=all_regions, .combine='c') %dopar% {
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




#*************************************************************************************************************
### IDEA 1:
# pour chaque TAD, je peux fitter une coexpr ~ dist
# et faire de même sur le même distance range avec des paires de gènes qui sont dans des TADs différents
# ensuite le score qqh comme AUC
#*************************************************************************************************************

#*************************************************************************************************************
### IDEA 2:
# pour chaque TAD, pour chaque paire de gènes, trouver une paire qui est dans une distance similaire (+- un radius)
# et calculer qqh comme un log2FC et prendre la moyenne (pour chaque paire, puis pour le TAD)
#*************************************************************************************************************

#*************************************************************************************************************
### IDEA 3:
# similaire à IDEA 2 mais travailler avec des classes de distances ?
#*************************************************************************************************************

#*************************************************************************************************************
### IDEA 4:
# soit g-1, g-2, g-3, ... les gènes à gauche du TAD, et g+1, g+2, g+3 les gènes à droite du TAD
# et g1, g2, g3, ... les gènes dans le TAD
# calculer les meanCorr pour les sets
# m1={g1,g2,g3,..gn} m2={g2,g3,..gn, g+1} mk={gi, ..., g+j} tant que g+j est dans le même intervalle de TAD size
# ensuite faire qqh comme un cumdiff
#*************************************************************************************************************

