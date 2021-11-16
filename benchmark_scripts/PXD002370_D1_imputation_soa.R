message("Preprocessing PXD002370-D1")

library(limma)
library(pcaMethods)
library(PEMM)
library(missForest)
library(mice)

dda <- read.delim("/stornext/General/data/academic/lab_davis/prot/benchmarking/PXD002370/txt/evidence.txt", stringsAsFactors = FALSE, header = TRUE)


table(grepl("CON__|REV__", dda$Leading.Razor.Protein))
table(grepl("ups", dda$Leading.Razor.Protein))
table(dda$Modifications)
table(dda$Type)


keep1 <- (!grepl("CON__|REV__", dda$Leading.Razor.Protein))
keep2 <- (dda$Modifications %in% "Unmodified")
keep3 <- (dda$Type != "MSMS")

table(keep1&keep2&keep3)


dda <- dda[keep1&keep2,]


# dda$Leading.Razor.Protein <- gsub("sp\\|","", dda$Leading.Razor.Protein)
# dda$Leading.Razor.Protein <- gsub("(.*)\\|(.*)","\\1", dda$Leading.Razor.Protein)

length(unique(dda$Experiment))

# these peptides have more than 1 record per sample
table(table(paste(gsub("_", "" , dda$Modified.sequence), dda$Charge, dda$Leading.Razor.Protein, sep = "_"))>6)


dda$PeptideID <- paste(gsub("_", "" , dda$Modified.sequence), dda$Charge, dda$Leading.Razor.Protein, sep = "_")


y_dda <- aggregate(Intensity ~ Experiment + PeptideID,
                   #FUN=function(x) sum(x[!duplicated(x)], na.rm = TRUE),
                   FUN = function(x) median(x, na.rm = TRUE),
                   na.action = na.pass, data = dda)

y_dda <- tidyr::spread(y_dda, key = Experiment, value = Intensity)
rownames(y_dda) <- y_dda$PeptideID
y_dda$PeptideID <- NULL

# pcv <- plotCV2(log2(y_dda), outlier = FALSE)
# mtext("data")
# text(pcv$mean[pcv$CV > 0.02], pcv$CV[pcv$CV>0.02], col="red", label=rownames(pcv)[pcv$CV> 0.02])

#######################################################################
### lets now completely remove the peptides with more than one record
#####################################################################

# tbl <- table(paste(gsub("_", "" , dda$Modified.sequence), dda$Charge, dda$Leading.Razor.Protein, sep = "_"))
# duppep <- names(tbl)[which(tbl>6)]
# 
# m <- dda[!dda$PeptideID %in% duppep, c("Experiment","PeptideID","Intensity")]
# 
# 
# y_dda <- aggregate(Intensity ~ Experiment + PeptideID,
#                     #FUN=function(x) sum(x[!duplicated(x)], na.rm = TRUE),
#                FUN = function(x) x,
#                na.action = na.pass, data = m)
# 
# y_dda$Intensity <- sapply(y_dda$Intensity, FUN=function(x) ifelse(is.null(x),NA,x))
# 
# 
# y_dda <- spread(y_dda, key = Experiment, value = Intensity)
# head(y_dda)
# 
# pcv <- plotCV2(log2(y_dda[,-1]), outlier = FALSE)
# 
# 
# # annot <- y[,1]
# 
# rownames(y_dda) <- y_dda[,1]
# y_dda <- y_dda[,-1]
#head(y)

dim(y_dda)

table(complete.cases(y_dda))


keep <- rowSums(!is.na(y_dda)) >= 4
table(keep)

y_dda <- y_dda[keep,]
y_dda <- normalizeBetweenArrays(log2(y_dda), method = "quantile")
colnames(y_dda) <- gsub("-","_", colnames(y_dda))


message("LLS")
Sys.time()
y_lls <- llsImpute(t(y_dda))
y_lls <- t(completeObs(y_lls))

message("Imputation in progress")
message("RF")
Sys.time()
y_rf <- missForest(t(y_dda))
y_rf <- t(y_rf$ximp)



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

saveRDS(list(y_rf, y_mice,y_em,y_bpca,y_lls), file="impute_PXD002370_D1_sOa.rds")
