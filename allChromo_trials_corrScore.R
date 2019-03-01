SSHFS <- FALSE
setDir <- ifelse(SSHFS, "~/media/electron", "")

startTime <- Sys.time()

suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, 40))

source("corrScore_fcts.R")

dataFold <- file.path(setDir, "/mnt/etemp/marie/Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(dataFold))

hicds <- "GSE105318_DLD1_40kb"
exprds <- "TCGAcoad_msi_mss"

# chromo <- "chr20"

stopifnot(dir.exists(file.path(dataFold, hicds)))

pipOutFold <- file.path(dataFold, "PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
stopifnot(dir.exists(file.path(pipOutFold)))

outFold <- file.path(paste0("TRIALS_CORRSCORE"))
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

genecoexpr_file <- file.path(dataFold,"CREATE_COEXPR_SORTNODUP", hicds, paste0(exprds, "_", corrMethodCoexpr), "coexprDT.Rdata")
stopifnot(file.exists(genecoexpr_file))
cat(paste0("... load geneCoexpr data\t", Sys.time(), " - "))
geneCoexpr_DT <- eval(parse(text = load(genecoexpr_file)))
cat(paste0(Sys.time(), "\n"))

genedist_file <- file.path(dataFold, "CREATE_DIST_SORTNODUP", hicds,"all_dist_pairs.Rdata")
stopifnot(file.exists(genedist_file))
cat(paste0("... load geneDist data\t", Sys.time(), " - "))
geneDist_DT <- eval(parse(text = load(genedist_file)))
cat(paste0(Sys.time(), "\n"))

genesametad_file <- file.path(dataFold, "CREATE_SAME_TAD_SORTNODUP", hicds,"all_TAD_pairs.Rdata")
stopifnot(file.exists(genesametad_file))
cat(paste0("... load geneSameTAD data\t", Sys.time(), " - "))
geneSameTAD_DT <- eval(parse(text = load(genesametad_file)))
cat(paste0(Sys.time(), "\n"))

# geneCoexpr_DT <- eval(parse(text = load(file.path(paste0(chromo, "_geneCoexpr_DT.Rdata")))))
# geneDist_DT <- eval(parse(text = load(file.path(paste0(chromo, "_geneDist_DT.Rdata")))))
# geneSameTAD_DT <- eval(parse(text = load(file.path(paste0(chromo, "_geneSameTAD_DT.Rdata")))))

### take only the filtered data according to initial settings
pipeline_regionList <- eval(parse(text = load(file.path(pipOutFold, script0_name, "pipeline_regionList.Rdata"))))

gene2tadDT <- gene2tadDT[gene2tadDT$region %in% pipeline_regionList,]

all_regions <- unique(as.character(gene2tadDT$region))
stopifnot(grepl("_TAD", all_regions))

# all_regions <- all_regions[grepl(paste0(chromo, "_"), all_regions)]

reg <- all_regions[35]

# all_regions <- all_regions[1:3]

all_scores_all_TADs <- foreach(reg = all_regions) %dopar% {
  
  cat(paste0("... START: ", reg, "\n"))
  
  reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
  
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
  
  # stopifnot(tadRange_coexprDistDT$chromo == tad_chromo)
  # => take all chromos
  
  stopifnot( nrow(tadRange_coexprDistDT) > 0 )
  
  #*************************************************************************************************************
  ### IDEA 1
  #*************************************************************************************************************
  
  
  score1 <- tad_score_trial1(
    tad_coexprDist_DT=tad_coexprDistDT,
    tadRange_coexprDist_DT=tadRange_coexprDistDT
  ) 
  
    
    
  #*************************************************************************************************************
  ### IDEA 2
  #*************************************************************************************************************
  
    score2 <- tad_score_trial2(
      tad_coexprDist_DT=tad_coexprDistDT,
      tadRange_coexprDist_DT=tadRange_coexprDistDT,
      distRad=1000
    )
    
      
  
  #*************************************************************************************************************
  ### IDEA 4
  #*************************************************************************************************************
  
      score4 <- tad_score_trial4(
        tad=reg,
        gt2_DT=gene2tadDT,
        expr_DT=norm_rnaseqDT
      )
      
      list(score1=score1, score2=score2, score4=score4)
      
}

#*************************************************************************************************************
#*************************************************************************************************************
#*************************************************************************************************************
#*************************************************************************************************************

names(all_scores_all_TADs) <- all_regions

outFile <- file.path(outFold, "all_scores_all_TADs.Rdata")
save(all_scores_all_TADs, file= outFile)
cat(paste0("... written: ", outFile, "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")

cat(paste0("*** DONE\n"))
cat(txt)
