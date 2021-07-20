library(DESeq2)

dds <- makeExampleDESeqDataSet(betaSD=1,interceptMean=10)
dds$batch <- factor(rep(c("A","B"),each=6))

vsd <- vst(dds)
plotPCA(vsd, "batch")

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
plotPCA(vsd, "batch")
