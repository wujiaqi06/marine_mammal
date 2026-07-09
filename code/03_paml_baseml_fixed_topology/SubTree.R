library("ape")
library("motmot")
read.tree("speciesTree.nwk") -> tree0
species0 <- scan("species.sub.txt","characters")

trait0 <- data.frame(species = species0)
trait0$trait <- sample(c(0,1),length(species0), replace = T)
rownames(trait0) <- species0

if (length(tree0$tip.label) == length(species0)){
  write.tree(tree0, file = "tree.tre")
} else{
  data0 <- sortTraitData(phy = tree0, y = trait0, 
                         data.name = "trait", pass.ultrametric = TRUE)
  phy <- data0$phy
  write.tree(phy, file = "tree.tre")
}


