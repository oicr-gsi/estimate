# ===================
# preprocess function
# ===================
preProcRNA <- function(gepFile, ensFile){
 
 # read in data
 print("reading gep")
 gepData <- read.csv(gepFile, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
 print("reading ensemble conversion")
 ensConv <- read.csv(ensFile, sep="\t", header=FALSE, stringsAsFactors=FALSE)

 # rename columns
 print("renaming columns")
 colnames(ensConv) <- c("gene_id", "Hugo_Symbol")

 # merge in Hugo's, re-order columns, deduplicate
 print("merge HUGOs")
 df <- merge(x=gepData, y=ensConv, by="gene_id", all.x=TRUE)
 print("subset HUGOs, deduplicate")
 df <- subset(df[,c(ncol(df),2:(ncol(df)-1))], !duplicated(df[,c(ncol(df),2:(ncol(df)-1))][,1]))

 print("setting rownames")
 row.names(df) <- df$Hugo_Symbol
 print("dropping rownames")
 df$Hugo_Symbol <- NULL
 # return the data frame
 return(df)
}

# ======================
# simple zscore function
# ======================
compZ <- function(df) {
  df_zscore <- df

  # require at least one sample has a gene with > 3
  df_zscore <- df_zscore[apply(df_zscore, 1, function(x) any(x > 3)), ]
  
  # log the data
  df_zscore <- log(df_zscore+1)

  # apply mean and SD
  df_zscore$MEAN <- apply(df_zscore, 1, mean)
  df_zscore$STDV <- apply(df_zscore, 1, sd)

  # compute zscore for each column
  for (i in c(1:(ncol(df_zscore)-2)))
   {
	  df_zscore[,i] <- (df_zscore[,i] - as.numeric(df_zscore$MEAN)) / as.numeric(df_zscore$STDV)
   }
 
  # drop the columns that we don't need anymore
  df_zscore <- df_zscore[,c(1:(ncol(df_zscore)-2))]
  return(df_zscore)
 }

# better zscore function:
compZ2 <- function(df, dft) {
  # define the mahalanobis function
  signMahal <- function(x) {
  # x: column 1 tumour, column 2 normal
  ST=sd(x[,1]); SN=sd(x[,2]); SNT=cov(x)[1,2]; MET=sqrt((SN^2)/((SN^2*ST^2)-(SNT)^2))
  MD <- (x[,1] - x[,2]) * MET
  return(MD)
 }
 # first, make a matrix of the common genes
 comGenes <- intersect(as.character(row.names(df)), as.character(row.names(dft)))
 GE <- df[row.names(df) %in% comGenes,][ order(row.names(df[row.names(df) %in% comGenes,])), ]
 SE <- dft[row.names(dft) %in% comGenes,][ order(row.names(dft[row.names(dft) %in% comGenes,])), ]

 # apply mean and SD over diploid samples
 dfz <- data.frame()
 pb <- txtProgressBar(min=0, max=length(comGenes), style=3)
 for (i in 1:length(comGenes)) {
   df_s <- GE[i,][colnames(GE)[colnames(GE) %in% colnames(SE[,SE[i,]==0])]]     # subset samples when CN=0 (diploid)
   medmad=data.frame(cbind(MEAN=apply(df_s, 1, mean), STDV=apply(df_s, 1, sd))) # apply mean and SD over these samples
   dfz=rbind(dfz, medmad)
   setTxtProgressBar(pb, i)
 }
 close(pb)

 # apply mean and SD over all samples without CN:
 GE2 <- df[!(row.names(df) %in% comGenes),]
 medmad=data.frame(cbind(MEAN=apply(GE2, 1, mean), STDV=apply(GE2, 1, sd)))

 # construct final dfz
 dfz=rbind(dfz, medmad)
 dfz$STDV <- ifelse(dfz$STDV == 0, 1, dfz$STDV)
 dfz <- merge(df, dfz, by="row.names", all.x=TRUE)
 row.names(dfz) <- dfz[,1]
 dfz <- dfz[,-1]
 
 # for each column, transform to zscore, and apply mahalanobis equation
 df_zscore <- dfz
 df_mahalo <- dfz
 for (i in c(1:(ncol(dfz)-3)))
  {
	  df_zscore[,i] <- (df_zscore[,i] - as.numeric(df_zscore$MEAN)) / as.numeric(df_zscore$STDV)
    df_mahalo[,i] <- signMahal(df_mahalo[,c(i,ncol(df_mahalo))])
  }
 # drop the columns that we don't need anymore
 df_zscore <- df_zscore[,c(1:(ncol(df_zscore)-2))]
 df_mahalo <- df_mahalo[,c(1:(ncol(df_mahalo)-2))]
 
 # make list and return
 ZSCORES=list()
 ZSCORES[[1]] <- df_zscore
 ZSCORES[[2]] <- df_mahalo
 return(ZSCORES)
}
