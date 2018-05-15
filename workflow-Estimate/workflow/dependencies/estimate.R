library(estimate)
library(GSVA)
library(GSEABase)
library(limma)

args <- commandArgs(trailingOnly=TRUE)

# inputs
expr_data <- args[1]
outdir <- args[2]
ensFile <- args[3]
gmtFile <- args[4]

# source RODiC preProcessing function
source("convert_rsem_results_zscore.r")

# fix the functions (or else names will have very annoying periods):
filterCommonGenes <- function (input.f, output.f, id = c("GeneSymbol", "EntrezID")) 
 {
    stopifnot((is.character(input.f) && length(input.f) == 1 && 
        nzchar(input.f)) || (inherits(input.f, "connection") && 
        isOpen(input.f, "r")))
    stopifnot((is.character(output.f) && length(output.f) == 
        1 && nzchar(output.f)))
    id <- match.arg(id)
    input.df <- read.table(input.f, header = TRUE, row.names = 1, 
        sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
    merged.df <- merge(common_genes, input.df, by.x = id, by.y = "row.names")
    rownames(merged.df) <- merged.df$GeneSymbol
    merged.df <- merged.df[, -1:-ncol(common_genes)]
    print(sprintf("Merged dataset includes %d genes (%d mismatched).", 
        nrow(merged.df), nrow(common_genes) - nrow(merged.df)))
    outputGCT(merged.df, output.f)
 }

outputGCT <- function (input.f, output.f) 
 {
    if (is.data.frame(input.f) == TRUE) {
        exp.data <- input.f
    }
    else {
        exp.data <- read.table(input.f, header = TRUE, row.names = 1, 
            sep = "\t", quote = "", check.names = FALSE)
    }
    exp.data1 <- data.frame(NAME = rownames(exp.data), Description = rownames(exp.data), 
        exp.data, check.names=FALSE)
    column1 <- colnames(exp.data1)
    column1[1] <- "NAME"
    column1[2] <- "Description"
    exp.data1$NAME <- factor(exp.data1$NAME)
    exp.data1$Description <- factor(exp.data1$Description)
    levels(exp.data1[, 1]) <- c(levels(exp.data1[, 1]), "NAME")
    levels(exp.data1[, 2]) <- c(levels(exp.data1[, 2]), "Description")
    exp.data2 <- rbind(column1, exp.data1)
    row1 <- rep("", length(1:ncol(exp.data)))
    row1_2 <- data.frame(row1, row1, check.names=FALSE)
    row1_2 <- t(row1_2)
    No_gene <- nrow(exp.data1)
    No_sample <- (ncol(exp.data1) - 2)
    GCT <- matrix(c("#1.2", No_gene, "", No_sample), nrow = 2, 
        ncol = 2)
    gct <- cbind(GCT, row1_2)
    colnames(gct) <- colnames(exp.data2)
    tmp <- rbind(gct, exp.data2)
    write.table(tmp, output.f, sep = "\t", row.names = FALSE, 
        col.names = FALSE, quote = FALSE)
    invisible(NULL)
 }

estimateScore <- function (input.ds, output.ds, platform = c("affymetrix", "agilent", 
    "illumina")) 
 {
    stopifnot(is.character(input.ds) && length(input.ds) == 1 && 
        nzchar(input.ds))
    stopifnot(is.character(output.ds) && length(output.ds) == 
        1 && nzchar(output.ds))
    platform <- match.arg(platform)
    ds <- read.delim(input.ds, header = TRUE, sep = "\t", skip = 2, 
        row.names = 1, blank.lines.skip = TRUE, as.is = TRUE, 
        na.strings = "", check.names=FALSE)
    descs <- ds[, 1]
    ds <- ds[-1]
    row.names <- row.names(ds)
    names <- names(ds)
    dataset <- list(ds = ds, row.names = row.names, descs = descs, 
        names = names)
    m <- data.matrix(dataset$ds)
    gene.names <- dataset$row.names
    sample.names <- dataset$names
    Ns <- length(m[1, ])
    Ng <- length(m[, 1])
    temp <- strsplit(input.ds, split = "/")
    s <- length(temp[[1]])
    input.file.name <- temp[[1]][s]
    temp <- strsplit(input.file.name, split = ".gct")
    input.file.prefix <- temp[[1]][1]
    for (j in 1:Ns) {
        m[, j] <- rank(m[, j], ties.method = "average")
    }
    m <- 10000 * m/Ng
    gs <- as.matrix(SI_geneset[, -1], dimnames = NULL)
    N.gs <- 2
    gs.names <- row.names(SI_geneset)
    score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
    for (gs.i in 1:N.gs) {
        gene.set <- gs[gs.i, ]
        gene.overlap <- intersect(gene.set, gene.names)
        print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", 
            length(gene.overlap)))
        if (length(gene.overlap) == 0) {
            score.matrix[gs.i, ] <- rep(NA, Ns)
            next
        }
        else {
            ES.vector <- vector(length = Ns)
            for (S.index in 1:Ns) {
                gene.list <- order(m[, S.index], decreasing = TRUE)
                gene.set2 <- match(gene.overlap, gene.names)
                correl.vector <- m[gene.list, S.index]
                TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
                no.TAG <- 1 - TAG
                N <- length(gene.list)
                Nh <- length(gene.set2)
                Nm <- N - Nh
                correl.vector <- abs(correl.vector)^0.25
                sum.correl <- sum(correl.vector[TAG == 1])
                P0 <- no.TAG/Nm
                F0 <- cumsum(P0)
                Pn <- TAG * correl.vector/sum.correl
                Fn <- cumsum(Pn)
                RES <- Fn - F0
                max.ES <- max(RES)
                min.ES <- min(RES)
                if (max.ES > -min.ES) {
                  arg.ES <- which.max(RES)
                }
                else {
                  arg.ES <- which.min(RES)
                }
                ES <- sum(RES)
                EnrichmentScore <- list(ES = ES, arg.ES = arg.ES, 
                  RES = RES, indicator = TAG)
                ES.vector[S.index] <- EnrichmentScore$ES
            }
            score.matrix[gs.i, ] <- ES.vector
        }
    }
    score.data <- data.frame(score.matrix)
    names(score.data) <- sample.names
    row.names(score.data) <- gs.names
    estimate.score <- apply(score.data, 2, sum)
    if (platform != "affymetrix") {
        score.data <- rbind(score.data, estimate.score)
        rownames(score.data) <- c("StromalScore", "ImmuneScore", 
            "ESTIMATEScore")
    }
    else {
        convert_row_estimate_score_to_tumor_purity <- function(x) {
            stopifnot(is.numeric(x))
            cos(0.6049872018 + 0.0001467884 * x)
        }
        est.new <- NULL
        for (i in 1:length(estimate.score)) {
            est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
            est.new <- rbind(est.new, est_i)
            if (est_i >= 0) {
                next
            }
            else {
                message(paste(sample.names[i], ": out of bounds", 
                  sep = ""))
            }
        }
        colnames(est.new) <- c("TumorPurity")
        estimate.t1 <- cbind(estimate.score, est.new)
        x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] < 
            0
        estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
        score.data <- rbind(score.data, t(estimate.t1))
        rownames(score.data) <- c("StromalScore", "ImmuneScore", 
            "ESTIMATEScore", "TumorPurity")
    }
    outputGCT(score.data, output.ds)
 }

# preProcessing
df <- preProcRNA(expr_data, ensFile)
df <- log2(df+1)
write.table(df, file=paste0(outdir, "/", basename(expr_data), ".pre"), sep="\t", quote=FALSE)
filterCommonGenes(input.f=paste0(outdir, "/", basename(expr_data), ".pre"), output.f=paste0(outdir, "/", basename(expr_data), ".gct"), id="GeneSymbol")

# get the score (illumina)
estimateScore(paste0(outdir, "/", basename(expr_data), ".gct"), paste0(outdir, "/", basename(expr_data), ".estimate.gct"), platform="illumina")

### ssGSEA ###
#gmt_data <- getGmt(gmtFile)
gmt_data <- as.matrix(read.delim(gmtFile, sep="\t", header=FALSE, stringsAsFactors=FALSE, row.names=1)[,-1])
gmt_data <- setNames(split(gmt_data, seq(nrow(gmt_data))), rownames(gmt_data))
z <- as.matrix(df)
z <- log2(z+1)
gsva_es <- gsva(z, gmt_data, mx.diff=1, method="ssgsea")
write.table(gsva_es, file=paste0(outdir, "/", basename(expr_data), ".ssGSEA.txt"), sep="\t", quote=FALSE)
