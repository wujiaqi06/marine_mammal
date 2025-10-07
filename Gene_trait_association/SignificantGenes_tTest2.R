library(ape)
library(motmot)
library(castor)
source("../marine.function.R")

#trait_name <- "marine"
trait_number <- 1 ##1 for marine, 2 for marine prediction (excluding 3 river dolphins), 3 for aquatic
clade_name <- "mammal"

gene_branch_interaction <- read.table(gzfile("../gene_branch_interaction.17432genes.txt.gz"))
dim(gene_branch_interaction);gene_branch_interaction[c(1:4),c(1:3)]
tree <- read.tree("../mammal302.anno.BL_support.nwk")
trait_all <- read.table("../MarineTrait.trait3.txt");head(trait_all)
trait_name <- colnames(trait_all)[trait_number]
trait_name
#### read clade
clade_branch <- read.delim("../mammal.branch.txt", row.names = 2)
terminals <- clade_branch[clade_branch$Species != "internal",1]
phy <- tree
####subdata
gene_branch_interaction[c(1:3),c(1:4)]
gene_branch_interaction <- gene_branch_interaction[,rownames(clade_branch)]
dim(gene_branch_interaction)
####ancestral states
traits <- trait_all[phy$tip.label,trait_name] + 1
number_find <- sum(as.numeric(traits == 2))
Nstates <- length(unique(traits))
anc.results = asr_max_parsimony(phy, traits, Nstates)
node_states = max.col(anc.results$ancestral_likelihoods)

traits.number <- data.frame(traits = traits, node_number = c(1:length(traits)))
anc.node.number <- data.frame(traits = node_states, node_number = c((length(traits)+1):(length(traits)+phy$Nnode)))

traits.states <- rbind(traits.number, anc.node.number)
rownames(traits.states) <- paste("B", traits.states$node_number, sep = "")
dim(traits.states);head(traits.states)

phy.edges <- data.frame(phy$edge)
rownames(phy.edges) <- paste("B",phy.edges[,2], sep ="")
colnames(phy.edges) <- c("anc", "offspring")
phy.edges$node_states <- traits.states[rownames(phy.edges),1]
states_info <- phy.edges$node_states
plot(phy, no.margin = TRUE, edge.width = 2,
     show.tip.label = T, edge.color = (phy.edges$node_states * 2), 
     cex = 0.3)
title(trait_name, line = -1.5)
traits <- traits -1
size <- 2
tiplabels(pch = 21, bg = rgb(0.1333333,0.5921569,0.9019608, 0.7), col = rgb(0,0,0,0.5), cex = traits*size, adj = 1.4)

phy <- read.tree("../mammal302.anno.nwk")
phy.edge <- data.frame(phy$edge, phy$edge.length, (1:nrow(phy$edge)))
colnames(phy.edge) <- c("anc", "offspring", "branch_label", "order")
phy.edge$branch_label <- paste("B",phy.edge$branch_label, sep = "")

Ape_Bl_march <- phy.edge
rownames(Ape_Bl_march) <- paste(Ape_Bl_march$anc, Ape_Bl_march$offspring, sep = "-")
head(Ape_Bl_march)
write.table(Ape_Bl_march, file = paste(trait_name,".",clade_name,".Ape_Bl_march.txt", sep = ""), quote = F, row.names = F, sep = "\t")
rownames(phy.edges) <- paste(phy.edges$anc,phy.edges$offspring, sep = "-")
head(phy.edges)
trait_anc <- data.frame(Ape_Bl_march, node_states = phy.edges[rownames(Ape_Bl_march),3])
trait_anc$node_states <- trait_anc$node_states - 1
rownames(trait_anc) <- trait_anc$branch_label
head(trait_anc);dim(trait_anc)
trait_anc.screen <- trait_anc[trait_anc$node_states != 0.5,]
dim(trait_anc.screen)

gene_branch_interaction_trait <- rbind(trait = trait_anc[colnames(gene_branch_interaction),"node_states"], gene_branch_interaction)
gene_branch_interaction_trait <- gene_branch_interaction_trait[,rownames(trait_anc.screen)]
dim(gene_branch_interaction_trait)
gene_branch_interaction_trait[c(1:4),c(1:3)]
table(as.numeric(gene_branch_interaction_trait[1,]))
####
trait.t <- data.frame()
for (i in 2:nrow(gene_branch_interaction_trait)){
  data1<- data.frame(gene_branch_interaction_trait[c(1,i),])
  data11 <- t(na.omit(t(data1)))
  trait_count <- table(data11[1,])
  if ((dim(trait_count) > 1)&(min(trait_count) > 1)){
    trait.t0 <- t.test(data11[2,]~data11[1,])
    trait.t.value <- trait.t0$statistic
    trait.p.value <- trait.t0$p.value
    trait.mean1 <- trait.t0$estimate[1]
    trait.mean2 <- trait.t0$estimate[2]
    trait.static.line <- data.frame(gene = rownames(gene_branch_interaction_trait)[i], 
                                    tvalue = trait.t.value, pvalue = trait.p.value, 
                                    mean1 = trait.mean1, mean2 = trait.mean2, 
                                    marine = trait_count[2], non_marine = trait_count[1])
    trait.t <- rbind(trait.t, trait.static.line)
  }
}
head(trait.t);dim(trait.t)
rownames(trait.t) <- trait.t$gene
par(mar = c(5, 4, 4, 2) + 0.1)
hist(trait.t$tvalue, breaks = 100, col = 2)
head(trait.t)
trait.t <- trait.t[order(trait.t$pvalue, decreasing = F),]
write.table(trait.t, file = paste(trait_name,".trait.t_test.all.trait3.17241.txt", sep = ""), quote = F, row.names = F, sep = "\t")

####FDR
####  FDR by Benjamini-Hochberg procedure
alpha <- 0.01
trait.pvalue_sort <- sort(trait.t$pvalue,na.last=NA)
trait.pvalue_order <- order(trait.t$pvalue,na.last=NA)  

gene.list <- NULL
p.list <- NULL

m <- length(trait.pvalue_sort)
k <- 1
while (trait.pvalue_sort[k] < k/m*alpha){
  p.list <- c(trait.pvalue_sort[k])
  gene.list <- c(gene.list,trait.pvalue_order[k])
  k <- k+1
}
length(gene.list)
t.gene <- trait.t[gene.list,]
write.table(t.gene, file = paste(trait_name,".",clade_name,".FDR",alpha, ".trait3.t_test.txt", sep = ""), quote = F, sep = "\t", row.names = F)
t.gene.sig <- trait.t[t.gene$gene,]
t.gene.sig <- t.gene.sig[order(t.gene.sig$tvalue),]
hist(t.gene.sig$tvalue, breaks = 100, col = 2, border=F,)

####
col0 = ifelse(trait.t$gene %in% t.gene.sig$gene, rgb(0.2,0.8,0.5,0.5) , ifelse (!(t.gene$gene %in% t.gene.sig$gene), "purple", rgb(0.2,0.2,0.2,0.2) ))

##traiting data
clade_branch <- read.delim("../mammal.branch.txt")
terminal <- clade_branch[clade_branch$Species != "internal",]

gene_branch_interaction_sig <- t(gene_branch_interaction[t.gene.sig$gene,])
gene_branch_interaction_sig[c(1:5),c(1:5)]
dim(gene_branch_interaction_sig)
gene_branch_interaction_sig.terminal <- as.matrix(gene_branch_interaction_sig[terminal$Branch,])
gene_branch_interaction_sig.terminal[c(1:4),c(1:4)];dim(gene_branch_interaction_sig.terminal)

terminal_mean <- apply(gene_branch_interaction_sig.terminal,2,mean_noNA)
head(terminal_mean)
rate.all.imputate <- imputation_by_column_terminal(gene_branch_interaction_sig,terminal_mean)
rate.all.imputate[c(1:5),c(1:5)]
dim(rate.all.imputate)
write.table(rate.all.imputate, file = paste("gene_branch_interaction.trait.imputated.t_test",".txt", sep = ""), quote = F, sep = "\t")
