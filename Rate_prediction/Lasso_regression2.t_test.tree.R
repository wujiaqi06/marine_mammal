library(glmnet)
library(ROCR)
library(minpack.lm)
library(ape)

trait_number <- 2
clade_name <- "mammal"
alpha <- 0.01


trait_all <- read.table("../MarineTrait.trait3.txt");head(trait_all);dim(trait_all)
#trait_all <- unique(trait_all);dim(trait_all)
trait_name <- colnames(trait_all)[trait_number];trait_name

gene_branch_interaction <- read.table("../marine.mammal.FDR0.01.gene_branch_interaction.imputated.trait3.1504.t_test.txt")
gene_branch_interaction[c(1:8),c(1:5)]
gene_branch_interaction <- scale(gene_branch_interaction)
dim(gene_branch_interaction);gene_branch_interaction[c(1:8),c(1:5)]
n <- ncol(gene_branch_interaction)

head(apply(gene_branch_interaction,1,sum))
head(apply(gene_branch_interaction,2,sum))


branches <- read.delim("../mammal.branch.txt");head(branches)
terminal <- branches[branches$Species != "internal",]
head(trait_all)
trait_all <- trait_all[terminal$Species,];head(trait_all)
#trait_all <- unique(trait_all)
rownames(terminal) <- terminal$Species
setdiff(terminal$Species, rownames(trait_all))

trait_all$bl <- terminal[rownames(trait_all),"Branch"]
#trait_all <- trait_all[!(is.na(trait_all$bl)),]
#trait_all <- unique(trait_all)
head(trait_all)
rownames(trait_all) <- trait_all$bl

gene_branch_interaction.terminal <- gene_branch_interaction[terminal$Branch,]
trait_all <- trait_all[rownames(gene_branch_interaction.terminal),]
dim(gene_branch_interaction.terminal)

trait <- trait_all
rate <- gene_branch_interaction.terminal
rate.pois.rate <- rate
dim(rate);rate[c(1:3),c(1:4)]

trait0 <- trait[trait[,trait_name] != 0.5,];dim(trait0)
trait0.test <- trait[trait[,trait_name] != 0.5,trait_name]

p0 <- mean(trait0.test)
rate00 <- as.data.frame(rate[trait[,trait_number] != 0.5,]);dim(rate00)
score0 <- abs((trait0.test - p0) %*% as.matrix(rate00))
#hist(score0)
gene0 <- order(score0,decreasing=T)
rate0 <- rate00[,gene0[1:n]]
gene.list <- gene0[1:n]
rate01 <- rate0
dim(rate01);rate01[c(1:3),c(1:3)];dim(rate01)


gene0 <- order(score0,decreasing=T)
rate0 <- rate00[,gene0[1:n]];dim(rate0)
gene.list <- gene0[1:n]
#rate00 <- as.data.frame(rate0[trait[,trait_number] != 0.5,]);dim(rate00)
#rate01 <- rate00[,gene.list]
#dim(rate01)
trait0.cv.glmnet <- cv.glmnet(x=as.matrix(rate0),trait0.test,family="binomial",nfolds = length(trait0.test))
trait0.lambda.min <- trait0.cv.glmnet$lambda.min
trait0.glmnet <- glmnet(x=as.matrix(rate0),trait0.test,lambda=trait0.lambda.min,"binomial")
trait0.coef <- trait0.glmnet$beta
trait0.coef <- as.matrix(trait0.coef)
#trait0.coef <- trait0.coef[trait0.coef[,1] != 0,]

write.table(trait0.coef, file = paste(trait_name,".",clade_name,".FDR",alpha,".n", n,".coef",".t_test.trait3.txt", sep = ""), quote = F, sep = "\t")

length(gene.list)
all_rates <- gene_branch_interaction[,gene.list]
trait0.predict <- predict(trait0.glmnet,newx=as.matrix(all_rates),type="response")
write.table(trait0.predict, file = paste(trait_name,".predict.n",n,".marine2.t_test.txt", sep = ""), quote = F, sep = "\t")

#####a0
names(trait0.glmnet)
trait0.glmnet$a0
trait0.glmnet$lambda
trait0.predict <- predict(trait0.glmnet,newx=as.matrix(all_rates),type="response")

marine_like <- log(trait0.predict/(1-trait0.predict))
write.table(marine_like, file = "linear_predictor.aquatic2.LASSO.txt", quote = F, row.names = T, col.names = F, sep = "\t")

####



###Plot prediction
####plot tree
trait_plot <- trait0.predict
tree0 <- read.tree("../mammal302.anno.BL_support.nwk")
node_support <- read.delim("../mammal.branch_support.txt")
common_name <- read.delim("../common_name.species302.txt");dim(common_name)
common_name <- unique(common_name);dim(common_name)
common_name <- common_name[common_name$Species %in% tree0$tip.label,]
dim(common_name)
rownames(common_name) <- common_name$Species
###Plot trait order
node_hide <- node_support[node_support$Support <= 50,]
dim(node_hide)
ape_order <- read.delim("../marine.mammal.Ape_Bl_march.txt")
head(ape_order)
trait_order <- trait_plot[ape_order$branch_label,]
head(trait_order)

####
rbPal <- colorRampPalette(c('#74add1','#b2182b'))

x <- data.frame(trait = trait_order)
head(x)
x0 <- (x-min(x, na.rm = T))/(max(x, na.rm = T))
#x0[is.na(x0)] <- 0.5
x0 <- as.matrix(x0);head(x0)
col0 <- rbPal(length(x0))[as.numeric(cut(x0,breaks = length(x0)))]
#col0[rownames(x0) %in% node_hide$Branch] <- "#636363"

##tips
tips <- data.frame(trait = trait[,trait_number])
rownames(tips) <- rownames(trait)
phy <- read.tree("../mammal302.anno.nwk")
tip_order <- data.frame(species = phy$tip.label, ApeOrder = c(1:length(phy$tip.label)))
rownames(tip_order) <- tip_order$species

branches <- read.delim("../mammal.branch.txt");head(branches)
terminal <- branches[branches$Species != "internal",]

rownames(terminal) <- terminal$Species
tip_order$BranchLabel <- terminal[tip_order$species,2]
rownames(tip_order) <- tip_order$BranchLabel
tip_order$traits <- tips[rownames(tip_order),1]
head(tip_order)

y <- data.frame(trait = tip_order$traits)
head(y)
y0 <- (y-min(y, na.rm = T))/(max(y, na.rm = T))
y0[is.na(y0)] <- 0.5
y0 <- as.matrix(y0);head(node_hide)

col_NA <- rgb( 0.5, 0, 1-0.5)
col1 <- rbPal(length(y0))[as.numeric(cut(y0,breaks = length(y0)))]

plot(tree0, no.margin = T, edge.color = col0, 
     direction = "r",#type = "u",
     align.tip.label = F,
     cex = 0.25, tip.color = col1)

