message("Preprocessing PXD016647")

library(limma)
library(pcaMethods)
library(PEMM)
library(missForest)
library(mice)

devtools::load_all("/stornext/General/data/academic/lab_davis/prot/benchmarking/msImpute/")

data <- read.delim("/stornext/General/data/academic/lab_davis/prot/benchmarking/PXD016647/Spike-in-bio-var-OT/combined/txt/evidence.txt", 
                   stringsAsFactors = FALSE)


# data <- data[grep("CON_|REV_", data$Leading.razor.protein, invert=TRUE),]
data <- data[data$Charge > 1,]
data$PeptideID <- paste0(data$Modified.sequence, data$Charge)
data$matrix.row.id <- paste(data$PeptideID, data$Leading.Razor.Protein, sep ="_")


genes <- data[,c("PeptideID","matrix.row.id", "Leading.razor.protein")]
genes <- genes[!duplicated(genes),]


y <- evidenceToMatrix(data)

genes <- genes[match(rownames(y), genes$PeptideID),]

y <- log2(y)
y <- y[rowSums(!is.na(y)) >=4,]

#ensure column names are allowed names
colnames(y) <- gsub("-","_", colnames(y))

message("LLS")
Sys.time()
y_lls <- llsImpute(t(y))
y_lls <- t(completeObs(y_lls))

message("Imputation in progress")
message("RF")
Sys.time()
y_rf <- missForest(t(y))
y_rf <- t(y_rf$ximp)



message("Imputation in progress (column-wise imputation)")
message("mice")
Sys.time()
# From mice -----

y_mice <- mice(y)
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

saveRDS(list(y_rf, y_mice,y_em,y_bpca,y_lls), file="~/softImpute_low_rank_experimentation/impute_PXD016647_sOa_280621_dda.rds")
