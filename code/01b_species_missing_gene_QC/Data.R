missing0 <- read.delim("mammal377.MissingCount.txt")
SingleCopy <- scan("6366Genes.txt", "characters")
missing0 <- missing0[order (missing0$MissingNumber, decreasing = F),]
dim(missing0[missing0$MissingNumber == 0,])

missing0$order <- c(1:nrow(missing0))

barplot(missing0$MissingNumber)
hist(missing0$MissingNumber, breaks = 100)
hist(missing0$MissingNumber, breaks = 1000)

head(missing0)
missing0.single_copy <- missing0[missing0$Sequence %in% SingleCopy,]
missing0.single_copy <- missing0.single_copy[order (missing0.single_copy$MissingNumber),]
missing0.single_copy$order <- c(1:nrow(missing0.single_copy))
write.table(missing0.single_copy, file = "species410.singleCopy.txt", quote = F, row.names = F, col.names = F)

##########
gene2275 <- scan("gene_qcPass.singleCopy.2275.geneList.txt", "characters")
gene2249 <- scan("all_tree.2249.success.txt", "characters")

failed26 <- setdiff(gene2275,gene2249)
write.table(failed26, file = "failed26.list.txt", quote = F, row.names = F, col.names = F)
