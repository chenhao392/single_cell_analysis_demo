library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(DropletUtils)
library(GenomeInfoDb)

# load the local functions
source("./src/util.R")


# set up gene annotations 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"

# load the processed peaks as a seurat object
ccRCC.atac.set1<-peak2seurat("GSM4572187_Control1")

# infer gene activities from peaks
gene.activities <- GeneActivity(ccRCC.atac.set1)

# add the gene activity matrix to the Seurat object as a new assay 
ccRCC.atac.set1[['RNA']] <- CreateAssayObject(counts = gene.activities)

# normalize gene activity matrix
ccRCC.atac.set1 <- NormalizeData(object = ccRCC.atac.set1,
                                  assay = 'RNA',
                                  normalization.method = 'LogNormalize',
                                  scale.factor = median(ccRCC.atac.set1$nCount_RNA)
                                 )



# export the inferred gene activity and peaks in MTX format
write10xCounts(x = ccRCC.atac.set1@assays$RNA@counts, path = "ccRCC_RNA")
write10xCounts(x = ccRCC.atac.set1@assays$peaks@counts, path = "ccRCC_ATAC")

#post processing for generating feature files
genes<-rownames(ccRCC.atac.set1@assays$RNA@counts)
features<-getGenesFeatures(genes)
write.table(features,file = "ccRCC_RNA/features.tsv",sep = "\t",row.names = T,col.names = F,quote = F)

peaks<-rownames(ccRCC.atac.set1@assays$peaks@counts)
features<-getPeaksFeatures(peaks)
write.table(features,file = "ccRCC_ATAC/features.tsv",sep = "\t",row.names = T,col.names = F,quote = F)
