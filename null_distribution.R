library(plyr)
library(parallel)
library(doParallel)
library(foreach)
nodes <- read.table(paste(disease,"nodes.txt",sep = "_"),header = F)
nodes <- nodes$V1
cal_randompp <- function(candidate_gene,risk_gene_eQTL){
  iter <- 0
  dif <- 100
  num0 <- rep(0,n)
  thres <- 1e-4
  gene_name <- candidate_gene$official_name
  ind <- c(1:n)
  sample_size <- 30
  CMC <- risk_gene_eQTL[order(risk_gene_eQTL$Pvalue,decreasing = F)[1:sample_size],]
  pick <- which(is.element(candidate_gene$Name,CMC$Name),arr.ind = T)
  freq0 <- rep(0,n)
  while(dif>thres){
    iter <- iter+1
    num <- num0
    pick <- sample(ind[-pick],sample_size,replace = F,prob = candidate_gene$weight[-pick])
    num[pick] <- num[pick]+1
    if(iter>1) freq0 <- num0/((iter-1)*sample_size)
    freq <- num/(iter*sample_size)
    dif <- norm(freq-freq0,"2")
    num0 <- num
  }
  list(freq)
}
cl <- makeCluster(5)
registerDoParallel(cl)

cat("random sampling start!")
load("E:/hedan/brain diseases/eQTL/cis_and_trans_eQTL.RData")
pathi <- paste(disease,"_candidate_genes.txt",sep="")
candidate_gene <- read.delim(pathi,header = T)
risk_gene_eQTL <- merge(eQTL,candidate_gene[,c("Name","SNP")],by=c("Name","SNP"))
risk_gene_eQTL <- unique(risk_gene_eQTL)
candidate_gene <- candidate_gene[is.element(candidate_gene$official_name,nodes),]
n <- nrow(candidate_gene)
candidate_gene$weight <- 1/n
set.seed(2021)
gene <- list()
gene <- foreach(x=1:10000,.combine='c') %dopar% cal_randompp(candidate_gene,risk_gene_eQTL)
gene <- data.frame(gene)
colnames(gene) <- paste("sample",1:10000,sep = "")
gene$gene_symbol <- candidate_gene$official_name
save(gene,file = paste(disease,"_random.RData",sep = ""))

stopCluster(cl)

