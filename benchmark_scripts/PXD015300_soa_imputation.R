message("Preprocessing PXD015300")

library(limma)
library(pcaMethods)
library(PEMM)
library(missForest)
library(mice)
library(tidyr)
devtools::load_all("/stornext/General/data/academic/lab_davis/prot/benchmarking/msImpute/")


dda <- read.delim("/stornext/General/data/academic/lab_davis/prot/benchmarking/PXD015300/evidence.tsv",
                  stringsAsFactors = FALSE)



keep1 <- (!grepl("CON__|REV__", dda$Leading.razor.protein))
keep2 <- (dda$Modifications %in% "Unmodified")



dda <- dda[keep1&keep2,]
dda <- dda[dda$Charge > 1,]
dda$PeptideID <- paste0(dda$Modified.sequence, dda$Charge)
dda$matrix.row.id <- paste(dda$PeptideID, dda$Leading.razor.protein, sep ="_")


genes <- dda[,c("PeptideID","matrix.row.id", "Leading.razor.protein")]
genes <- genes[!duplicated(genes),]


y_dda <- evidenceToMatrix(dda)

keep3 <- rowSums(!is.na(y_dda)) >= 4
table(keep3)



y_dda <- y_dda[keep3,]

#ensure column names are allowed names
colnames(y_dda) <- gsub("-","_", colnames(y_dda))
genes <- genes[match(rownames(y_dda), genes$PeptideID),]
rownames(y_dda) <- genes$matrix.row.id




y_dda <- log2(y_dda)
# y_dda <- normalizeBetweenArrays(y_dda, method = "quantile")
colnames(y_dda) <- gsub("-","_", colnames(y_dda))


message("LLS")
Sys.time()
y_lls <- llsImpute(t(y_dda))
y_lls <- t(completeObs(y_lls))

message("Imputation in progress")
message("RF")
Sys.time()
y_rf <- missForest(y_dda)
y_rf <- y_rf$ximp



message("Imputation in progress (column-wise imputation)")
message("mice")
Sys.time()
# From mice -----



y_mice <- mice(y_dda)
y_mice <- complete(y_mice)



message("EM")
Sys.time()
# from PEMM ----
y_em <- PEMM_fun(y_dda, phi=0)$Xhat


message("bpca")
Sys.time()
# From pcaMethods --------
y_bpca <- pca(y_dda, nPcs=2, method="bpca")
y_bpca <- completeObs(y_bpca)

saveRDS(list(y_rf, y_mice,y_em,y_bpca,y_lls), file="~/softImpute_low_rank_experimentation/impute_PXD015300_sOa_filter4obs_noNorm.rds")
