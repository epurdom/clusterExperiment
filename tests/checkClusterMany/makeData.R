library(brainData)
load("../data/L5_NOIMPUTE_FQ_NOWEIGHT_NOBIO_BATCH_Q_3.rda")
fq <- scone_out$exp_mat
fq[fq<0] <- 0

data(cortical)
cortical <- cortical[,colnames(fq)]

qc <- as.matrix(colData(cortical)[,1:14])
batch <- droplevels(cortical$MD_c1_run_id)

stopifnot(all(rownames(qc)==colnames(fq)))

glia0 <-intersect(c("DSJN001_N710_S505", "DSJN001_N702_S510", "DSJN003_N724_S521", "DSJN003_N721_S522", "DSJN003_N726_S522", "DSJN003_N728_S518", "DSJN003_N727_S517", "DSJN003_N728_S520", "DSJN003_N728_S515", "DSJN003_N728_S513", "DS07_N703_S506", "DS07_N703_S502", "DS07_N702_S506", "DS07_N701_S505", "DS07_N703_S505", "DS07_N701_S508", "DS07_N702_S505", "DS07_N704_S506", "DS07_N714_S506", "DS07_N715_S505", "DSJN001_N724_S503", "DSJN002_N703_S507", "DSJN002_N702_S502", "DSJN002_N705_S508", "DSJN002_N702_S510", "DSJN002_N711_S502", "DSJN003_N701_S507", "DSJN002_N706_S511", "DSJN002_N715_S502", "DSJN002_N712_S502", "DSJN003_N705_S506", "DSJN003_N704_S507", "DSJN002_N724_S515", "DSJN002_N722_S513", "DSJN002_N720_S521", "DSJN002_N719_S517", "DSJN002_N718_S521", "DSJN002_N728_S513", "DSJN003_N718_S513", "DSJN002_N728_S516", "DS08_N719_S503", "DS08_N723_S507"), colnames(fq))

glia <- rep("Neuron", length(batch))
names(glia) <- colnames(fq)
glia[glia0] <- "Glia"
glia <- as.factor(glia)

qcpca <- prcomp(qc, scale. = TRUE, center = TRUE)

mar16 <- read.table("results_160302/clustering_layer5/cluster_labels.txt", row.names=1)
mar16 <- mar16[colnames(fq),1]

mar16_merged <- read.table("results_160302/clustering_layer5/cluster_labels_merged.txt", row.names=1)
mar16_merged <- mar16_merged[colnames(fq),1]

l5 <- SummarizedExperiment(fq, colData=DataFrame(Glia=glia,
                                                 Batch=batch,
                                                 Mar_clusters=mar16,
                                                 Mar_merged=mar16_merged))
save(l5,file="L5_sumExp.rda")		
print(sessionInfo())									