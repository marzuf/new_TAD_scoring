SSHFS <- FALSE
setDir <- ifelse(SSHFS, "~/media/electron", "")

startTime <- Sys.time()

suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, 40))

source("corrScore_fcts.R")
source("../Cancer_HiC_data_TAD_DA/utils_plot_fcts.R")
plot_multiDens_setcols(size_list, plotTit="", legTxt=NULL, legPos="topright", my_ylab="density", my_xlab="", my_cols = NULL) 

dataFold <- file.path(setDir, "/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(dataFold))

hicds <- "GSE105318_DLD1_40kb"
exprds <- "TCGAcoad_msi_mss"

chromo <- "chr20"

stopifnot(dir.exists(file.path(dataFold, hicds)))

pipOutFold <- file.path(dataFold, "PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
stopifnot(dir.exists(file.path(pipOutFold)))

outFold <- file.path(paste0("TRIALS_CORRSCORE_", chromo))
dir.create(outFold, recursive = TRUE)

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"


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

geneCoexpr_DT <- eval(parse(text = load(file.path(paste0(chromo, "_geneCoexpr_DT.Rdata")))))

geneDist_DT <- eval(parse(text = load(file.path(paste0(chromo, "_geneDist_DT.Rdata")))))

geneSameTAD_DT <- eval(parse(text = load(file.path(paste0(chromo, "_geneSameTAD_DT.Rdata")))))

### take only the filtered data according to initial settings
pipeline_regionList <- eval(parse(text = load(file.path(pipOutFold, script0_name, "pipeline_regionList.Rdata"))))

gene2tadDT <- gene2tadDT[gene2tadDT$region %in% pipeline_regionList,]

all_regions <- unique(as.character(gene2tadDT$region))
stopifnot(grepl("_TAD", all_regions))

all_regions <- all_regions[grepl(paste0(chromo, "_"), all_regions)]

reg <- all_regions[1]
#*************************************************************************************************************


outFile <- file.path(outFold, "all_scores_all_TADs.Rdata")
all_scores_all_TADs <- eval(parse(text=load(outFile)))

head(names(all_scores_all_TADs))
all_scores_tad1 <- all_scores_all_TADs[[1]]
names(all_scores_tad1)
tad1_score1 <- all_scores_tad1[["score1"]]
# "auc_inTAD_obsDist"   "auc_diffTAD_obsDist" "slope_inTAD"         "slope_diffTAD"    
names(tad1_score1)
tad1_score2 <- all_scores_tad1[["score2"]]
names(tad1_score2[[1]])
# "obs_coexpr"     "diffAvg_coexpr"
tad1_score4 <- all_scores_tad1[["score4"]]
names(tad1_score4)
# "tad_corr"           "left_window_corrs"  "right_window_corrs"


### SCORES 1

all_scores_all_TADs[[1]][[1]][["auc_inTAD_obsDist"]]

score1_all_auc_inTAD <- lapply(all_scores_all_TADs, function(tad_scores_all_TADs) tad_scores_all_TADs[["score1"]][["auc_inTAD_obsDist"]])
score1_all_auc_diffTAD <- lapply(all_scores_all_TADs, function(tad_scores_all_TADs) tad_scores_all_TADs[["score1"]][["auc_diffTAD_obsDist"]])
score1_all_slope_inTAD <- lapply(all_scores_all_TADs, function(tad_scores_all_TADs) tad_scores_all_TADs[["score1"]][["slope_inTAD"]])
score1_all_slope_diffTAD <- lapply(all_scores_all_TADs, function(tad_scores_all_TADs) tad_scores_all_TADs[["score1"]][["slope_diffTAD"]])


plot_multiDens_setcols(
  list(auc_inTAD = unlist(score1_all_auc_inTAD),
       auc_diffTAD = unlist(score1_all_auc_diffTAD))
)

plot_multiDens_setcols(
  list(slope_inTAD = unlist(score1_all_slope_inTAD),
       slope_diffTAD = unlist(score1_all_slope_diffTAD))
)


### SCORES 2

all_scores_all_TADs[[1]][[2]][[1]][["obs_coexpr"]]

score2_all_obs_coexpr <- lapply(all_scores_all_TADs, function(tad_scores_all_TADs) 
                                        lapply(tad_scores_all_TADs[["score2"]], function(tad_idx) 
                                            tad_idx[["obs_coexpr"]]))

score2_all_diffAvg_coexpr <- lapply(all_scores_all_TADs, function(tad_scores_all_TADs) 
  lapply(tad_scores_all_TADs[["score2"]], function(tad_idx) 
    tad_idx[["diffAvg_coexpr"]]))


plot_multiDens_setcols(
  list(obs_coexpr = unlist(score2_all_obs_coexpr),
       diffAvg_coexpr = unlist(score2_all_diffAvg_coexpr))
)


### SCORES 4

all_scores_all_TADs[[4]][[1]][["auc_inTAD_obsDist"]]

score4_all_tad_corr <- lapply(all_scores_all_TADs, function(tad_scores_all_TADs) tad_scores_all_TADs[["score4"]][["tad_corr"]])
score4_all_left_window_corrs <- lapply(all_scores_all_TADs, function(tad_scores_all_TADs) tad_scores_all_TADs[["score4"]][["left_window_corrs"]])
score4_all_right_window_corrs <- lapply(all_scores_all_TADs, function(tad_scores_all_TADs) tad_scores_all_TADs[["score4"]][["right_window_corrs"]])

plot_multiDens_setcols(
  list(tad_corr = unlist(score4_all_tad_corr),
       left_window_corrs = unlist(score4_all_left_window_corrs),
       right_window_corrs = unlist(score4_all_right_window_corrs))
)


#*************************************************************************************************************
#*************************************************************************************************************
#*************************************************************************************************************

txt <- paste0(startTime, "\n", Sys.time(), "\n")

cat(paste0("*** DONE\n"))
cat(txt)
