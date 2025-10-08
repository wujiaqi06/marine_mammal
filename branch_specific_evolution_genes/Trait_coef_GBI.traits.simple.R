###
d0 <- read.delim("trait.Pred.allgene.trait3.txt")
d0$diff <- d0$aquatic - d0$marine
branches <- read.delim("../mammal.branch.txt");head(branches)
internal <- branches[branches$Species == "internal",]
terminal <- branches[branches$Species != "internal",]



gene_branch_interaction <- read.delim(gzfile("../gene_branch_interaction.17432genes.txt.gz"))
#gene_branch_interaction[c(1:4),c(1:3)]
marine.gbi <- read.delim("../marine.mammal.FDR0.01.gene_branch_interaction.imputated.trait3.1504.t_test.txt")
marine.gbi <- scale(marine.gbi)
#marine.gbi <- data.frame(t(marine.gbi))
aquatic.gbi <- read.delim("../aquatic2.mammal.FDR0.01.gene_branch_interaction.imputated,trait3.1315.t_test.txt")
aquatic.gbi <- scale(aquatic.gbi)
#aquatic.gbi <- data.frame(t(aquatic.gbi))

marine.coef <- read.delim("Coefficient/marine2.mammal.FDR0.01.n1504.coef.t_test.trait3.txt")
marine.coef <- marine.coef[order(marine.coef$coef),]
rownames(marine.coef) <- marine.coef$gene
aquatic.coef <- read.delim("Coefficient/aquatic2.mammal.FDR0.01.n1315.coef.t_test.trait3.txt")
aquatic.coef <- aquatic.coef[order(aquatic.coef$coef),]
rownames(aquatic.coef) <- aquatic.coef$gene

gene_branch_interaction.marine <- t(marine.gbi[,marine.coef$gene])
dim(gene_branch_interaction.marine);gene_branch_interaction.marine[c(1:3),c(1:4)]

gene_branch_interaction.aquatic <- t(aquatic.gbi[,aquatic.coef$gene])
dim(gene_branch_interaction.aquatic);gene_branch_interaction.aquatic[c(1:3),c(1:4)]

all_genes <- unique(c(marine.coef$gene, aquatic.coef$gene))
all_genes.gbi_orginal <- gene_branch_interaction[all_genes,]
all_genes.gbi_orginal.missing <- matrix(as.numeric(!(is.na(all_genes.gbi_orginal))), nrow = nrow(all_genes.gbi_orginal))
rownames(all_genes.gbi_orginal.missing) <- rownames(all_genes.gbi_orginal)
colnames(all_genes.gbi_orginal.missing) <- colnames(all_genes.gbi_orginal)

marine.gbi_missing <- all_genes.gbi_orginal.missing[marine.coef$gene,]
aquatic.gbi_missing <- all_genes.gbi_orginal.missing[aquatic.coef$gene,]

a0_marine <- -3.527868
a0_aquatic <- -2.468869
par(mfrow = c(3, 3))

####Example, Humpback whale, with branch lable B537 in data
label <- "Humpback whale"
Branch_show <- "B537"
####
plot(d0$marine~d0$aquatic, cex = 1, col="#bababa",
     xlab = "Aquatic", ylab = "Marine",
     main = "Whales + hippopotamus",
     ylim=c(0,1), xlim=c(0,1), xaxs = "i",yaxs = "i")
abline(0,1,col="#9ecae1")

label1 <- paste(label, ", marine", sep = "")
trait_pred.clade <- d0[Branch_show,]
round(trait_pred.clade, 3)
points(trait_pred.clade$aquatic, trait_pred.clade$marine, 
       col='#b2182b', 
       pch=16, 
       cex =2)
text((trait_pred.clade$aquatic-0.1), (trait_pred.clade$marine-0.1),
     label, cex = 1)
trait_pred.clade

node <- Branch_show
gene_branch_interaction.marine.one <- data.frame(gbi = gene_branch_interaction.marine[,node])
rownames(gene_branch_interaction.marine.one) <- rownames(gene_branch_interaction.marine)
gene_branch_interaction.marine.one$coef <- marine.coef[rownames(gene_branch_interaction.marine),"coef"]
gene_branch_interaction.marine.one$coef_bgi <- gene_branch_interaction.marine.one$gbi * gene_branch_interaction.marine.one$coef
gene_branch_interaction.marine.one$missing <- marine.gbi_missing[rownames(gene_branch_interaction.marine),node]
gene_branch_interaction.marine.one$border_col <- "#1f78b4"
gene_branch_interaction.marine.one[gene_branch_interaction.marine.one$missing == 0,"border_col"] <- "#d9d9d9"
gene_branch_interaction.marine.one$col <- "#bdbdbd"
gene_branch_interaction.marine.one[gene_branch_interaction.marine.one$missing == 0,"col"] <- "#d9d9d9"

sum0 <- round(sum(gene_branch_interaction.marine.one$coef_bgi, na.rm = T),2);sum0 + a0_marine
max0 <- 2
min0 <- -1

barplot(gene_branch_interaction.marine.one$coef_bgi, border = gene_branch_interaction.marine.one$border_col, names.arg = "", horiz=F , 
        col = gene_branch_interaction.marine.one$col,
        ylim = c(min0,max0),
        main = label1,
        ylab = "GBI * coef", xlab = "Genes",
        las=2, cex.names = 0.6)

label1 <- paste(label, ", aqutic", sep = "")
gene_branch_interaction.aquatic.one <- data.frame(gbi = gene_branch_interaction.aquatic[,node])
rownames(gene_branch_interaction.aquatic.one) <- rownames(gene_branch_interaction.aquatic)
gene_branch_interaction.aquatic.one$coef <- aquatic.coef[rownames(gene_branch_interaction.aquatic),"coef"]
gene_branch_interaction.aquatic.one$coef_bgi <- gene_branch_interaction.aquatic.one$gbi * gene_branch_interaction.aquatic.one$coef
gene_branch_interaction.aquatic.one$missing <- aquatic.gbi_missing[rownames(gene_branch_interaction.aquatic),node]
gene_branch_interaction.aquatic.one$border_col <- "#1f78b4"
gene_branch_interaction.aquatic.one[gene_branch_interaction.aquatic.one$missing == 0,"border_col"] <- "#d9d9d9"
gene_branch_interaction.aquatic.one$col <- "#bdbdbd"
gene_branch_interaction.aquatic.one[gene_branch_interaction.aquatic.one$missing == 0,"col"] <- "#d9d9d9"

sum0 <- round(sum(gene_branch_interaction.aquatic.one$coef_bgi, na.rm = T),3);sum0 + a0_aquatic
max0 <- 2
min0 <- -1

barplot(gene_branch_interaction.aquatic.one$coef_bgi, border = gene_branch_interaction.aquatic.one$border_col, 
        names.arg = "", horiz=F ,
        col = gene_branch_interaction.aquatic.one$col,
        ylim = c(min0,max0),
        main = label1,
        ylab = "GBI * coef", xlab = "Genes",
        las=2, cex.names = 0.6)

