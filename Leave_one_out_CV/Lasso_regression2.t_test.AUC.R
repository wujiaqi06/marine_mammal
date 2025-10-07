library(glmnet)
library(ROCR)
library(minpack.lm)
library(ape)

trait_number <- 1
clade_name <- "mammal"
alpha <- 0.01


trait_all <- read.table("../MarineTrait.trait3.txt");head(trait_all);dim(trait_all)
trait_name <- colnames(trait_all)[trait_number]

gene_branch_interaction <- read.table("../aquatic2.mammal.FDR0.01.gene_branch_interaction.imputated,trait3.1315.t_test.txt")
gene_branch_interaction[c(1:8),c(1:5)]
gene_branch_interaction <- scale(gene_branch_interaction)
dim(gene_branch_interaction);gene_branch_interaction[c(1:8),c(1:5)]
n <- ncol(gene_branch_interaction);n

head(apply(gene_branch_interaction,1,sum))
head(apply(gene_branch_interaction,2,sum))


branches <- read.delim("../mammal.branch.txt");head(branches)
terminal <- branches[branches$Species != "internal",]

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

## leave-one-out by genus
leave_one_pred2 <- c()
#j <-2
rownames(terminal) <- terminal$Branch
terminal <- terminal[rownames(rate01),]
genus <- c()
for (j in 1:nrow(terminal)){
  genus <- c(genus, strsplit(terminal$Species[j],"_")[[1]][1])
}
terminal$genus <- genus
genus_uniq <- unique(genus)
trait0.test <- data.frame(trait0.test)
rownames(trait0.test) <- terminal$Branch
trait0.test$genus <- terminal$genus
for (j in genus_uniq){
  lines_test <- which(trait0.test$genus %in% j)
  leave_one.trait = trait0.test[lines_test,1]
  leave_one.rate = rate01[lines_test,]
  remain.trait = trait0.test[-lines_test,1]
  remain.rate = rate01[-lines_test ,]
  trait.remain.cv.glmnet <- cv.glmnet(x=as.matrix(remain.rate),remain.trait,family="binomial",nfolds = length(remain.trait))
  trait.remain.lambda.min <- trait.remain.cv.glmnet$lambda.min
  trait.remain.glmnet <- glmnet(x=as.matrix(remain.rate),remain.trait,lambda=trait.remain.lambda.min,"binomial")
  trait.test.coef <- trait.remain.glmnet$beta
  one.pred <- predict(trait.remain.glmnet,newx=as.matrix(leave_one.rate),type="response")
  rbind(leave_one_pred2, one.pred) -> leave_one_pred2
}

data.all.pred <- prediction(leave_one_pred2,trait0.test[,1])
data.all.perf3 <- performance(data.all.pred, measure="acc", x.measure="cutoff")
bestAccInd <- which.max(data.all.perf3@"y.values"[[1]])
leave_one_out2_accuracy <- data.all.perf3@"y.values"[[1]][bestAccInd]
mis_class.rate <- 1- leave_one_out2_accuracy
leave_one_out2_accuracy.cut <- round(leave_one_out2_accuracy, digits = 3)
mis_class.rate.cut <- round(mis_class.rate,digits=3)
cutoff <- round(data.all.perf3@"x.values"[[1]][bestAccInd], 4)
data.all.perf <- performance(data.all.pred,"tpr","fpr")
data.all.auc <- performance(data.all.pred,"auc")
data.all.auc <- unlist(slot(data.all.auc, "y.values"));data.all.auc
data.all.auc1 <- round(data.all.auc,digits = 3)
data.all.pref4 <- performance(data.all.pred , measure = "acc", x.measure = "rec")
BestRecall <- which.max(data.all.pref4@"y.values"[[1]])
data.all_recall <-  data.all.pref4@"x.values"[[1]][BestRecall]
#recall <- c()
#cbind(colnames(trait)[i],paste("RR recall ",data.all_recall)) -> recall
data.all.auc
plot(data.all.perf,
     avg='vertical',
     spread.estimate='stderror',
     lwd=3,main=paste('17343 genes of 299 mammals, AUC =', data.all.auc1),
     col='blue')
abline(0,1)


#head(recall)

leave_one.trait.all <- cbind(leave_one_pred2,trait0.test)
head(leave_one.trait.all)
plot(leave_one.trait.all)
write.table(leave_one.trait.all, file = paste(trait_name,".",clade_name,".FDR",alpha,".n", n,".LOOCV",".t_test.trait3.txt", sep = ""), quote = F,sep = "\t")

p0 <- mean(trait0.test$trait0.test)
score0 <- abs((trait0.test$trait0.test - p0) %*% as.matrix(rate00))
#hist(score0)
n
gene0 <- order(score0,decreasing=T)
rate0 <- rate00[,gene0[1:n]];dim(rate0)
gene.list <- gene0[1:n]
trait0.cv.glmnet <- cv.glmnet(x=as.matrix(rate0),trait0.test$trait0.test,family="binomial",nfolds = length(trait0.test$trait0.test))
trait0.lambda.min <- trait0.cv.glmnet$lambda.min
trait0.glmnet <- glmnet(x=as.matrix(rate0),trait0.test$trait0.test,lambda=trait0.lambda.min,"binomial")
trait0.coef <- trait0.glmnet$beta
trait0.coef <- as.matrix(trait0.coef)

write.table(trait0.coef, file = paste(trait_name,".",clade_name,".FDR",alpha,".n", n,".coef",".t_test.trait3.txt", sep = ""), quote = F, sep = "\t")

length(gene.list)
all_rates <- gene_branch_interaction[,gene.list]
trait0.predict <- predict(trait0.glmnet,newx=as.matrix(all_rates),type="response")
write.table(trait0.predict, file = paste(trait_name,".predict.n",n,".marine2.t_test.txt", sep = ""), quote = F, sep = "\t")

