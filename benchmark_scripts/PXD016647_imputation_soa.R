message("Preprocessing PXD016647")

library(limma)
library(pcaMethods)
library(PEMM)
library(missForest)
library(mice)

dia <- read.delim("/stornext/General/data/academic/lab_davis/prot/benchmarking/PXD016647/Spike-in-biol-var-OT-SN-Report.txt", header = TRUE, stringsAsFactors = FALSE)


# table(grepl("CON__|REV__", dda$Leading.razor.protein))
# table(dda$Modifications)
# table(dda$Type)
# 
# table(grepl("ups", dda$Leading.razor.protein))
# table(grepl("iRT", dda$Gene.names))
# table(dda$Raw.file)

dia$experiment <- paste(dia$R.Condition, dia$R.Replicate,sep="_")
table(dia$experiment)

table(dia$R.Replicate)
table(dia$R.Condition)
table(grepl("UPS", dia$PG.ProteinAccessions))
table(grepl("CON|REV", dia$PG.ProteinAccessions))
table(grepl("ups", dia$PG.BGSFASTAHeader))
table(grepl("iRT", dia$PG.BGSUniProtId))
table(dia$SpikeIn=="True")

table(dia$PG.BGSDatabase)
sum(table(dia$PG.BGSDatabase)[grep("sp", names(table(dia$PG.BGSDatabase)), invert = TRUE)])

hist(dia$EG.Qvalue)
table(table(dia$EG.PrecursorId) > 25)

table(dia$IsProteotypic, dia$EG.Qvalue < 0.1)

keep1 <- (dia$IsProteotypic == "True")
keep2 <- (dia$EG.Qvalue < 0.1)

# table(keep1&keep2, dia$SpikeIn)


dia <- dia[keep1&keep2,]


# MS1 --------------
# dia_ms1 <- tidyr::spread(dia[,c("experiment","EG.PrecursorId", "FG.NormalizedMS1PeakArea")], key = experiment, value = FG.NormalizedMS1PeakArea)
# 
# hist(log2(data.matrix(dia_ms1[,-1]))) # the study says the limit of detection is 100, and they removed anything less
# # the limit of detection
# 
# table(log2(data.matrix(dia_ms1[,-1])) < log2(100))


# MS2 --------------

dia_ms2 <- tidyr::spread(dia[,c("experiment","EG.PrecursorId", "FG.NormalizedMS2PeakArea")], key = experiment, value = FG.NormalizedMS2PeakArea)

hist(log2(data.matrix(dia_ms2[,-1])))


# dia_ms1[dia_ms1==1] <- NA
dia_ms2[dia_ms2==1] <- NA

# par(mfrow=c(1,2))
# pcv <- plotCV2(normalizeBetweenArrays(log2(data.matrix(dia_ms2[,-1]))), outlier = FALSE)
# mtext("MS2")
# pcv <- plotCV2(normalizeBetweenArrays(log2(data.matrix(dia_ms1[,-1]))), outlier = FALSE)
# mtext("MS1")

y_dia_ms2 <- log2(data.matrix(dia_ms2[,-1]))
rownames(y_dia_ms2) <- dia_ms2$EG.PrecursorId

keep <- rowSums(!is.na(y_dia_ms2)) >=4
table(keep)

y_dia_ms2 <- y_dia_ms2[keep,]

#keep <- rowSums(!is.na(y_dia_ms1)) >=4
#table(keep)



# y_dia_ms2 <- normalizeBetweenArrays(y_dia_ms2, method = "quantile")

# y_dia_ms2 <- y_dia_ms2/colSums(y_dia_ms2, na.rm = TRUE)
# y_dia_ms2 <- y_dia_ms2*(10^6)


### 28 June 2021 : no normalization, qvalue filter, min 4 observation filter


message("LLS")
Sys.time()
y_lls <- llsImpute(t(y_dia_ms2))
y_lls <- t(completeObs(y_lls))

message("Imputation in progress")
message("RF")
Sys.time()
y_rf <- missForest(t(y_dia_ms2))
y_rf <- t(y_rf$ximp)



message("Imputation in progress (column-wise imputation)")
message("mice")
Sys.time()
# From mice -----



y_mice <- mice(y_dia_ms2)
y_mice <- complete(y_mice)



message("EM")
Sys.time()
# from PEMM ----
y_em <- PEMM_fun(y_dia_ms2, phi=0)$Xhat


message("bpca")
Sys.time()
# From pcaMethods --------
y_bpca <- pca(y_dia_ms2, nPcs=2, method="bpca")
y_bpca <- completeObs(y_bpca)

saveRDS(list(y_rf, y_mice,y_em,y_bpca,y_lls), file="impute_PXD016647dia_sOa_filter4obs_noNorm.rds")