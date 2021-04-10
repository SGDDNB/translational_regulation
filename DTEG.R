library(DESeq2)
library(ggplot2)
## Input files
# Calculating differential translation genes (DTGs) requires the count matrices from Ribo-seq and RNA-seq.
# These should be the raw counts obtained from feature counts or any other tool for counting reads, 
# they should not be normalized or batch corrected.
# It also requires a sample information file which should be in the same order as samples in the count
# matrices. It should include information on sequencing type, treatment, batch or any other covariate you 
# need to model.

### Get count matrix files and sample information
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
  stop("Ribo-seq counts, RNA-seq counts and Sample Information and batch presence (0/1) should be supplied.n", call.=FALSE)
} 

# Input filenames
ribo_file <- args[1]
rna_file <- args[2]
data_file <- args[3]
batch <- args[4]

# Read and merge count matrices
ribo <- read.delim(ribo_file)  
rna <- read.delim(rna_file)
merge <- cbind(ribo,rna)

head(merge)

# Sample information file
coldata <- read.delim(data_file)
coldata <- as.data.frame(apply(coldata,2,as.factor))

head(coldata)

## Detecting differential translation regulation
### DESeq2 object with batch and interaction term in the design

if(batch == 1){
  ddsMat <- DESeqDataSetFromMatrix(countData = merge,
                                   colData = coldata, design =~ Batch + Condition + SeqType + Condition:SeqType)
}else if(batch == 0){
  ddsMat <- DESeqDataSetFromMatrix(countData = merge,
                                   colData = coldata, design =~ Condition + SeqType + Condition:SeqType)
}else{
  stop("Batch presence should be indicated by 0 or 1 only", call.=FALSE)
}

ddsMat$SeqType = relevel(ddsMat$SeqType,"RNA")
ddsMat <- DESeq(ddsMat)
resultsNames(ddsMat)

system("mkdir Results")
setwd("Results")
system("mkdir fold_changes")
system("mkdir gene_lists")

# Choose the term you want to look at from resultsNames(ddsMat) 
# Condition2.SeqTypeRibo.seq means Changes in Ribo-seq levels in Condition2 vs
# Condition1 accounting for changes in RNA-seq levels in Condition2 vs Condition1
res <- results(ddsMat, contrast=list("Condition2.SeqTypeRIBO"))
summary(res)
length(which(res$padj < 0.05))
write.table(rownames(res)[which(res$padj < 0.05)],"gene_lists/DTEGs.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(res,"fold_changes/deltaTE.txt",quote=F,sep="\t",col.names = T,row.names = T)

pdf("Result_figures.pdf",useDingbats = F)

## Visualisation and interpretation

### DESeq2 object with batch for Ribo-seq
ind = which(coldata$SeqType == "RIBO")
coldata_ribo = coldata[ind,]

# PCA 
if(batch == 1){
  ddsMat_ribo <- DESeqDataSetFromMatrix(countData = ribo,
                                        colData = coldata_ribo, design =~ Condition + Batch)
  vsd <- vst(ddsMat_ribo)
  pcaData <- plotPCA(vsd, intgroup=c("Condition", "Batch"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+ theme(aspect.ratio=1)

  }else if(batch ==0){
  ddsMat_ribo <- DESeqDataSetFromMatrix(countData = ribo,
                                        colData = coldata_ribo, design =~ Condition)
  vsd <- vst(ddsMat_ribo)
  pcaData <- plotPCA(vsd, intgroup="Condition", returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+ theme(aspect.ratio=1)
}

ddsMat_ribo <- DESeq(ddsMat_ribo)
res_ribo <- results(ddsMat_ribo, contrast=c("Condition","2","1"))
res_ribo <- lfcShrink(ddsMat_ribo, coef=2,res=res_ribo,type="apeglm")
write.table(res_ribo,"fold_changes/deltaRibo.txt",quote=F,sep="\t",col.names = T,row.names = T)


### DESeq2 object with batch for RNA-seq
ind = which(coldata$SeqType == "RNA")
coldata_rna = coldata[ind,]
# PCA 
if(batch == 1){
    ddsMat_rna <- DESeqDataSetFromMatrix(countData = rna,
                                       colData = coldata_rna, design =~ Condition + Batch)
    vsd <- vst(ddsMat_rna)
    pcaData <- plotPCA(vsd, intgroup=c("Condition", "Batch"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Batch)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+ theme(aspect.ratio=1)
}else if(batch ==0){
  ddsMat_rna <- DESeqDataSetFromMatrix(countData = rna,
                                       colData = coldata_rna, design =~ Condition)
  vsd <- vst(ddsMat_rna)
  pcaData <- plotPCA(vsd, intgroup="Condition", returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(aspect.ratio=1)
}
ddsMat_rna <- DESeq(ddsMat_rna)


res_rna <- results(ddsMat_rna, contrast=c("Condition","2","1"))
res_rna <- lfcShrink(ddsMat_rna, coef=2,type="apeglm",res=res_rna)
write.table(res_rna,"fold_changes/deltaRNA.txt",quote=F,sep="\t",col.names = T,row.names = T)
write.table(rownames(res_rna)[which(res_rna$padj < 0.05)],"gene_lists/DTG.txt",quote=F,sep="\t",col.names = F,row.names = F)

## Classes of genes
forwarded = rownames(res)[which(res$padj > 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05)]
write.table(forwarded,"gene_lists/forwarded.txt",quote=F,sep="\t",col.names = F,row.names = F)
exclusive = rownames(res)[which(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj > 0.05)]
write.table(exclusive,"gene_lists/exclusive.txt",quote=F,sep="\t",col.names = F,row.names = F)
both = which(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05)
intensified = rownames(res)[both[which(res[both,2]*res_rna[both,2] > 0)]]
write.table(intensified,"gene_lists/intensified.txt",quote=F,sep="\t",col.names = F,row.names = F)
buffered = rownames(res)[both[which(res[both,2]*res_rna[both,2] < 0)]]
buffered = c(rownames(res)[which(res$padj < 0.05 & res_ribo$padj > 0.05 & res_rna$padj < 0.05)],buffered)
write.table(buffered,"gene_lists/buffered.txt",quote=F,sep="\t",col.names = F,row.names = F)

max_val = max(res_ribo[,2],res_rna[,2],na.rm = T)
plot(y=res_ribo[,2],x=res_rna[,2], xlab="RNA-seq log2 fold change",ylab = "Ribo-seq log2 fold change",asp=1,pch=16,col=rgb(128/255,128/255,128/255,0.1),ylim=c(-max_val,max_val),xlim=c(-max_val,max_val),cex=0.4)
abline(a=0,b=1,col="gray")
abline(h=0,v=0,col="gray")
points(y=res_ribo[forwarded,2],x=res_rna[forwarded,2],pch=16,col=rgb(0,0,1,1))
points(y=res_ribo[exclusive,2],x=res_rna[exclusive,2],pch=16,col=rgb(1,0,0,1))
points(y=res_ribo[intensified,2],x=res_rna[intensified,2],pch=16,col=rgb(1,0,1,1))
points(y=res_ribo[buffered,2],x=res_rna[buffered,2],pch=16,col=rgb(1,0,1,1))

### Examples for each class of genes
par(mfrow=c(2,2))
goi = forwarded[1]
y_u = max(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
y_l = min(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
plot(c(1,2),c(0,res[goi,2]),type="l",col="red",xaxt="n",
     xlab="Conditions",ylim=c(y_l,y_u),ylab="Log2 Fold Change", 
     main="Forwarded gene")
lines(c(1,2),c(0,res_ribo[goi,2]),col="gray")
lines(c(1,2),c(0,res_rna[goi,2]),col="blue")
axis(1,at=c(1,2),labels=c(1,2),las=1)
legend("bottomleft",c("RNA","Ribo","RibOnly"), fill=c("blue","gray","red"),
       cex=1, border = NA, bty="n")

goi = exclusive[1]
y_u = max(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
y_l = min(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
plot(c(1,2),c(0,res[goi,2]),type="l",col="red",xaxt="n",
     xlab="Conditions",ylim=c(y_l,y_u),ylab="Log2 Fold Change", 
     main="Exclusive gene")
lines(c(1,2),c(0,res_ribo[goi,2]),col="gray")
lines(c(1,2),c(0,res_rna[goi,2]),col="blue")
axis(1,at=c(1,2),labels=c(1,2),las=1)
legend("bottomleft",c("RNA","Ribo","RibOnly"), fill=c("blue","gray","red"),
       cex=1, border = NA, bty="n")

goi = buffered[1]
y_u = max(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
y_l = min(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
plot(c(1,2),c(0,res[goi,2]),type="l",col="red",xaxt="n",
     xlab="Conditions",ylim=c(y_l,y_u),ylab="Log2 Fold Change", 
     main="Buffered gene")
lines(c(1,2),c(0,res_ribo[goi,2]),col="gray")
lines(c(1,2),c(0,res_rna[goi,2]),col="blue")
axis(1,at=c(1,2),labels=c(1,2),las=1)
legend("bottomleft",c("RNA","Ribo","RibOnly"), fill=c("blue","gray","red"),
       cex=1, border = NA, bty="n")

goi = intensified[1]
y_u = max(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
y_l = min(res[goi,2],res_ribo[goi,2], res_rna[goi,2],0)
plot(c(1,2),c(0,res[goi,2]),type="l",col="red",xaxt="n",
     xlab="Conditions",ylim=c(y_l,y_u),ylab="Log2 Fold Change", 
     main="Intensified gene")
lines(c(1,2),c(0,res_ribo[goi,2]),col="gray")
lines(c(1,2),c(0,res_rna[goi,2]),col="blue")
axis(1,at=c(1,2),labels=c(1,2),las=1)
legend("bottomleft",c("RNA","Ribo","RibOnly"), fill=c("blue","gray","red"),
       cex=1, border = NA, bty="n")

dev.off()





