library(fmsb)
library(pROC)

setwd("./out/")
# get scores
fs <- list.files("./",pattern="sscore$")
dat <- read.table(fs[1],comment.char = "",header=T)
score <- data.frame(X=paste0(dat[,1],"_",dat[,1]))
for(i in 1:length(fs)){
	dat <- read.table(fs[i],comment.char = "",header=T)
	score <- cbind(score,dat[,4])
	colnames(score)[i+1] <- gsub(".sscore","",fs[i])
}

# count snp numbers 
dat <- read.table("profile_p",header=F)
dat2 <- read.table("../prange",header=F)
pthres <- dat2[,1]
count <- rep(0,length(pthres))
for(i in 1:length(pthres)){
	count[i] <- sum(dat[,2] <= pthres[i])
}

# order count for score names
scorenames <- as.numeric(gsub("score.","",names(score)[-1]))
idx <- match(scorenames,pthres)
count <- count[idx]

########################
# association test
########################

# read phenotype
phef <- "/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/phenotype/pheno_PC"
phe <- read.table(phef,header=T)
phe$BrNum[nchar(phe$BrNum) == 5] <- sub("^Br", "Br0", phe$BrNum[nchar(phe$BrNum) == 5])

phe <- phe[phe$Dx=="Schizo" | phe$Dx=="Control",]
phe <- phe[!is.na(phe$PCArace),]
phe$id <- paste0(phe[,2],"_",phe[,2])

########### EA 
phe1 <- phe[phe$PCArace=="CAUC",]

# reorder data
id.com <- intersect(phe1$id,score$X)
phe1 <- phe1[match(id.com,phe1$id),]
score1 <- score[match(id.com,score$X),]

data <- cbind(score1[,-1],phe1)
data$Dx <- as.factor(data$Dx) #(178 scz and 343 controls)

# association test
ps <- colnames(score)[-1]
pr = data.frame(array(0, c(length(ps),4)))
rownames(pr) <- ps
colnames(pr) <- c("count","p","r2","ROC")
pr[,1] <- count
for(i in 1:length(ps)){
	#fit0 = glm(Dx ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, family=binomial, data=data)
	#fit1 = glm(Dx ~ data[,i] + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, family=binomial, data=data)
	fit0 = glm(Dx ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=binomial, data=data)
	fit1 = glm(Dx ~ data[,i] + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=binomial, data=data)
	s = summary(fit1)
	d = NagelkerkeR2(fit1)$R2 - NagelkerkeR2(fit0)$R2;
	pr[i, 2] = s$coef[2,4]
	pr[i, 3] = d * 100
	roc_obj <- roc(data$Dx, data[,i])
	pr[i, 4] <- auc(roc_obj)
}
pr <- pr[order(pr[,3] * -1),]
write.csv(pr,"assoc_EA2.csv")

########### AA
phe1 <- phe[phe$PCArace=="AA",]

# reorder data
id.com <- intersect(phe1$id,score$X)
phe1 <- phe1[match(id.com,phe1$id),]
score1 <- score[match(id.com,score$X),]

data <- cbind(score1[,-1],phe1)
data$Dx <- as.factor(data$Dx)  #(130 scz and 276 controls)

# association test
ps <- colnames(score)[-1]
pr = data.frame(array(0, c(length(ps),4)))
rownames(pr) <- ps
colnames(pr) <- c("count","p","r2","ROC")
pr[,1] <- count
for(i in 1:length(ps)){
	#fit0 = glm(Dx ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, family=binomial, data=data)
	#fit1 = glm(Dx ~ data[,i] + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, family=binomial, data=data)
	fit0 = glm(Dx ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=binomial, data=data)
	fit1 = glm(Dx ~ data[,i] + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=binomial, data=data)
	s = summary(fit1)
	d = NagelkerkeR2(fit1)$R2 - NagelkerkeR2(fit0)$R2;
	pr[i, 2] = s$coef[2,4]
	pr[i, 3] = d * 100
	roc_obj <- roc(data$Dx, data[,i])
	pr[i, 4] <- auc(roc_obj)
}
pr <- pr[order(pr[,3] * -1),]
write.csv(pr,"assoc_AA2.csv")