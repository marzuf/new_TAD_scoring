SSHFS <- TRUE
setDir <- ifelse(SSHFS, "~/media/electron", "")

suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
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

outFold <- file.path("TRIALS_NEWPERMUT")
dir.create(outFold, recursive = TRUE)

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script_name <- "5v2_runPermutations"

diffTADcol <- "blue"
inTADcol <- "red"
plotCex <- 1.2

# norm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_qqnorm_rnaseqDT.Rdata"))))
initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

gene2tadDT_file <- file.path(dataFold, hicds, "genes2tad", "all_genes_positions.txt")
stopifnot(file.exists(gene2tadDT_file))
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(geneList),]


nGenesByTAD <- setNames(as.numeric(table(gene2tadDT$region)), as.character(names(table(gene2tadDT$region))))

tad_set <- names(nGenesByTAD)
random_tad_set <- sample(x=tad_set, size=length(tad_set), replace=FALSE)

tad="chr10_TAD1"

i=0
while(nrow(random_g2t) < nrow(gene2tadDT)){
  first_gene_idx <- 1
  last_gene_idx <- NA
  
  cat("... ", i, "\n")
  random_tad_set <- sample(x=tad_set, size=length(tad_set), replace=FALSE)
  
  random_g2t <- data.frame(entrezID = character(0),
                           chromo = character(0),
                           start = numeric(0),
                           end = numeric(0),
                           region = character(0),
                           stringsAsFactors = FALSE)
  for(tad in random_tad_set)  {
    # for(tad in "chr10_TAD1") {
    
    stopifnot(tad %in% names(nGenesByTAD))
    tad_nGenes <- as.numeric(nGenesByTAD[tad])
    stopifnot(!is.na(nGenes))
    last_gene_idx <- first_gene_idx + tad_nGenes - 1
    if(last_gene_idx > nrow(gene2tadDT))
      break
    tad_geneDT <- gene2tadDT[first_gene_idx:last_gene_idx,,drop=FALSE]
    stopifnot(nrow(tad_geneDT) == tad_nGenes)
    tad_geneDT$region <- tad
    
    obs_tad <- unique(gene2tadDT$region[gene2tadDT$entrezID %in% tad_geneDT$entrezID])
    
    obs_g2t <- gene2tadDT[gene2tadDT$region %in% obs_tad,]
    stopifnot(nrow(obs_g2t) > 0)
    # # discard if exactly the same TAD
    if(length(unique(obs_g2t$region)) == 1 & nrow(tad_geneDT) == nrow(obs_g2t) ) {
      break
    }
    # if(length(unique(obs_g2t$region)) == 1  ) {
    #   
    #   cat("break")
    # } 
    random_g2t <- rbind(random_g2t, tad_geneDT)
    
    if( nrow(random_g2t) == nrow(gene2tadDT)) break
    
    first_gene_idx <- last_gene_idx + 1

  }
  
  cat(paste0("...... nrow(random_g2t) = ", nrow(random_g2t), "\n"))
  
  cat(paste0("...... nrow(gene2tadDT) = ", nrow(gene2tadDT), "\n"))
  
  i=i+1
}
random_nGenesByTAD <- setNames(as.numeric(table(random_g2t$region)), as.character(names(table(random_g2t$region))))

stopifnot(!duplicated(random_g2t$entrezID))
stopifnot(random_nGenesByTAD == nGenesByTAD)






gene2tadDT



