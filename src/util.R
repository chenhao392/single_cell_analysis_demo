library(stringr)

############################################ 
# Read in the cellranger atac format peaks #
# then return a seurat object              #
############################################
peak2seurat<-function(fileBase){
  #files
  countFile<-paste("process/",fileBase,"_filtered_peak_bc_matrix.h5",sep="")
  fragFile<-paste("process/",fileBase,"_fragments.tsv.gz",sep="")
  #load files
  count = Read10X_h5(filename =countFile)
  
  #assay
  chrom_assay <- CreateChromatinAssay(
    counts = count,
    sep = c(":", "-"),
    genome = 'hg19',
    fragments = fragFile,
    min.cells = 10,
    min.features = 200
  )
  
  #seurat object
  seuObj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks"
  )
  
  # add the gene information to the object
  Annotation(seuObj) <- annotations
  return(seuObj)
}

############################################################
# read in a list of official gene symbols, then            #
# return a data frame with modality and genome coordinates #
############################################################
getGenesFeatures<-function(genes){
  
  #select genes
  isgene <- (annotations$gene_name %in% genes)
  annot.sub <- annotations[isgene]
  
  #select isoform with longest cds length
  annot.sub<-subset(annot.sub,type=="cds")
  anno.df<-data.frame(annot.sub)
  anno.df<-anno.df[order(anno.df$width,decreasing = T),]
  anno.uniq<-anno.df[!duplicated(anno.df$gene_name), ]
  
  #add row name
  rownames(anno.uniq)<-anno.uniq$gene_name
  
  #reorder the annotation according to the input gene list
  anno.uniq<-anno.uniq[genes,]
  #change rownames to EMSEMBL ID
  rownames(anno.uniq)<-genes

  #get the data frame according to the 10x feature format
  anno.coor<-anno.uniq[,1:3]
  anno.coor$modality<-"expression"
  anno.coor$ID<-rownames(anno.coor)
  anno.coor<-anno.coor[,c(5,4,1,2,3)]
  colnames(anno.coor)<-c("ID","modality","chr","start","end")
  
  return(anno.coor)
}


################################################################
# read in a list of genome coordinates in compact format, then #
# return a data frame with modality and genome coordinates     #
################################################################

getPeaksFeatures<-function(peaks){
  #string split
  anno.coor<-str_split_fixed(peaks, "-",3)
  
  #get the data frame according to the 10x feature format
  anno.coor<-data.frame(anno.coor)
  anno.coor$ID<-peaks
  anno.coor$modality<-"peaks"
  anno.coor<-anno.coor[,c(4,5,1,2,3)]
  colnames(anno.coor)<-c("ID","modality","chr","start","end")
  rownames(anno.coor)<-anno.coor$ID
  
  return(anno.coor)
}