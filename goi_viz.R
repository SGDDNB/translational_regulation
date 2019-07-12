## Gene of interest visalization

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("DESeq2 output for RPF, mRNA and TE, and gene ID", call.=FALSE)
} 

# Input filenames
ribo <- read.delim(args[1])
rna <- read.delim(args[2])
te <- read.delim(args[3])
gene <- as.character(args[4])

pdf(paste(gene,".pdf",sep=""))

ymax=max(ribo[gene,2],rna[gene,2],te[gene,2],0)
ymin=min(ribo[gene,2],rna[gene,2],te[gene,2],0)
plot(c(0,1), c(0,ribo[gene,2]), type="l", col="gray", ylim=c(ymin,ymax), ylab="Log2 fold change",xlab="",xaxt="n")
lines(c(0,1), c(0,rna[gene,2]), type="l", col="blue")
lines(c(0,1), c(0,te[gene,2]), type="l", col="red")
legend("bottomleft",c("RNA","Ribo","TE"), fill=c("blue","gray","red"),
       cex=1, border = NA, bty="n")
axis(1,at=c(0,1),labels=c(1,2),las=1)
dev.off()

