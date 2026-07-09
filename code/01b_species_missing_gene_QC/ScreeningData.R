gene_missing <- read.delim("species302.MissingCount.txt")
alignment_length <- read.table("seq_length.txt")
colnames(alignment_length) <- c("gene", "length")
single_copy <- scan("6366Genes.txt", "characters")

hist(gene_missing$MissingNumber, breaks = 100)

missing_taxa <- 1
n_alignment <- 1

gene_missing.screen <- gene_missing[gene_missing$MissingNumber < missing_taxa,]
dim(gene_missing.screen);head(gene_missing.screen)
write.table(gene_missing.screen$Sequence, file = "mammal302.noMissing.geneList.txt", quote = F, row.names = F, col.names = F)


alignment_length.screen <- alignment_length[alignment_length$length >= n_alignment,]
dim(alignment_length.screen)

gene_qcPass <- intersect(gene_missing.screen$Sequence, alignment_length.screen$gene)
length(gene_qcPass)
gene_qcPass.singleCopy <- intersect(gene_qcPass, single_copy)
length(gene_qcPass.singleCopy)

write.table(gene_qcPass.singleCopy, file = "mammal302.noMissing_singleCopy.geneList.txt", quote = F, row.names = F, col.names = F)

gene_missing.noMis <- gene_missing[gene_missing$MissingNumber == 0,]
dim(gene_missing.noMis)

read.delim()



