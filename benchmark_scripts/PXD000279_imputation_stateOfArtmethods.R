message("Preprocessing PXD000279")

library(limma)
library(pcaMethods)
library(PEMM)
library(missForest)
library(mice)



# ups_data <- read.delim("/stornext/General/data/academic/lab_davis/prot/benchmarking/PXD000279/dynamicrangebenchmark/peptides.txt", stringsAsFactors = FALSE)
# #ups_data <- as.data.frame(fread("PXD000279/dynamicrangebenchmark/modificationSpecificPeptides.txt", stringsAsFactors = FALSE))
# 
# ups_data <- ups_data[grep("CON__|REV__", ups_data$Leading.razor.protein, invert=TRUE),]
# 
# # Intensities are already summed up for each peak (i.e. sequence/charge combination)
# # ups_data$PeptideID <- paste(ups_data$Sequence, ups_data$Charges, ups_data$Leading.razor.protein, sep="_")
# ups_data$PeptideID <- paste(ups_data$Sequence, ups_data$Charges, ups_data$Leading.razor.protein, sep="_")
# y <- ups_data[,grep("Intensity.UPS", colnames(ups_data))]
# rownames(y) <- ups_data$PeptideID
# 
# # confirm UPS peptides are in the expression matrix
# 
# 
# 
# y[y == 0] <- NA
# keep <- rowSums(!is.na(y)) >= 4
# 
# 
# y <- y[keep,]
# y <- log2(as.matrix.data.frame(y))
# # y <- normalizeBetweenArrays(y)


devtools::load_all("/stornext/General/data/academic/lab_davis/prot/benchmarking/msImpute/")
data <- read.delim("/stornext/General/data/academic/lab_davis/prot/spikeIn_maxLFQpaper/txt/evidence.txt", stringsAsFactors = FALSE)


# data <- data[grep("CON_|REV_", data$Leading.razor.protein, invert=TRUE),]
data <- data[data$Charge > 1,]
data$PeptideID <- paste0(data$Modified.sequence, data$Charge)
data$matrix.row.id <- paste(data$PeptideID, data$Leading.razor.protein, sep ="_")


genes <- data[,c("PeptideID","matrix.row.id", "Leading.razor.protein","Modifications")]
genes <- genes[!duplicated(genes),]


y <- evidenceToMatrix(data)
rownames(y) <- genes$matrix.row.id[match(rownames(y), genes$PeptideID)]



genes <- genes[match(rownames(y), genes$matrix.row.id),]

keep1 <- (!grepl("CON__|REV__", genes$Leading.razor.protein))
keep2 <- (genes$Modifications %in% "Unmodified")



# y <- normalizeBetweenArrays(log2(y[keep1&keep2,]), method = "quantile")
y <- log2(y[keep1&keep2,])
print(dim(y))
y <- y[rowSums(!is.na(y)) >= 4,]

colnames(y) <- gsub("-","_", colnames(y))




message("LLS")
Sys.time()
y_lls <- llsImpute(t(y))
y_lls <- t(completeObs(y_lls))

message("Imputation in progress")
message("RF")
Sys.time()
y_rf <- missForest(y)
y_rf <- y_rf$ximp



message("Imputation in progress (column-wise imputation)")
message("mice")
Sys.time()
# From mice -----



y_mice <- mice(data.frame(y))
y_mice <- complete(y_mice)



message("EM")
Sys.time()
# from PEMM ----
y_em <- PEMM_fun(y, phi=0)$Xhat


message("bpca")
Sys.time()
# From pcaMethods --------
y_bpca <- pca(y, nPcs=2, method="bpca")
y_bpca <- completeObs(y_bpca)

saveRDS(list(y_rf, y_mice,y_em,y_bpca,y_lls), file="~/softImpute_low_rank_experimentation/impute_PXD000279_sOa_filter4obs_noNorm_WEHI.rds")

# v2  does column-wise imputation for all except LLE
