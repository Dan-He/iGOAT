library(kknn)
iGOAT <- function(disease, gene_path, flanking = 1000000){
  all_gene <- read.delim(gene_path,as.is=T)
  all_gene <- all_gene[!is.na(all_gene$official_name),]
  all_gene$Name <- substr(all_gene$Name,1,15)
  GWAS_path <- paste("GWAS_file/",disease,"_gwas.txt",sep = "")
  GWAS_loci <- read.delim(GWAS_path)
  if(disease%in%c("BP","AD","PD","MG")){
    GWAS_loci <- GWAS_loci[GWAS_loci$P.VALUE<5e-8,]
    GWAS_loci <- GWAS_loci[match(unique(GWAS_loci$SNP),GWAS_loci$SNP),]
  }
  GWAS_loci$start <- unlist(sapply(GWAS_loci$Pos_hg19,function(x) max(x-250000,0)))
  GWAS_loci$end <- GWAS_loci$Pos_hg19+250000 #250kb
  if(substr(GWAS_loci$chr[1],1,3)!="chr"){GWAS_loci$chr <- paste("chr",GWAS_loci$chr,sep = "")}
  # collect candidate genes based on 2D distance
  gene1 <- list()
  for(i in 1:nrow(GWAS_loci)){
    start <- max(GWAS_loci[i,]$Pos_hg19 - flanking,0)
    end <- GWAS_loci[i,]$Pos_hg19 + flanking # flanking = 1000000
    temp <- all_gene[all_gene$chrom == GWAS_loci[i,]$chr,]
    index1 <- temp$start_hg19 < end & temp$start_hg19 > start 
    index2 <- temp$end_hg19 < end & temp$end_hg19 > start
    temp <- temp[index1|index2,]
    if(nrow(temp)>0)
    {
      temp$SNP <- GWAS_loci[i,]$SNP
      temp$SNP_chr <- GWAS_loci[i,]$chr
      temp$SNP_pos_hg19 <- GWAS_loci[i,]$Pos_hg19 
      gene1 <- rbind(gene1,temp)
    }
  } 
  
  # collect candidate genes based on 3D interaction
  brainCP <- read.table("HiC/S22_TSS_CP.txt",as.is = T,header = T)
  brainGZ <- read.table("HiC/S23_TSS_GZ.txt",as.is = T,header = T)
  gene2 <- NULL
  for(i in 1:nrow(GWAS_loci)){
    chr <- GWAS_loci$chr[i]
    cpi <- brainCP[brainCP$chr==chr,]
    gzi <- brainGZ[brainGZ$chr==chr,]
    ind1 <- cpi$interacting_bin_start<GWAS_loci$end[i] & cpi$interacting_bin_start>GWAS_loci$start[i]
    ind2 <- cpi$interacting_bin_end<GWAS_loci$end[i] & cpi$interacting_bin_end>GWAS_loci$start[i]
    cpi <- cpi[ind1|ind2,]
    ind1 <- gzi$interacting_bin_start<GWAS_loci$end[i] & gzi$interacting_bin_start>GWAS_loci$start[i]
    ind2 <- gzi$interacting_bin_end<GWAS_loci$end[i] & gzi$interacting_bin_end>GWAS_loci$start[i]
    gzi <- gzi[ind1|ind2,]
    reg_gene <- union(cpi$ENSGID_for_TSS,gzi$ENSGID_for_TSS)
    if(sum(is.element(all_gene$Name,reg_gene))>0){
      genei <- all_gene[is.element(all_gene$Name,reg_gene),]
      genei$SNP <- GWAS_loci$SNP[i]
      genei$SNP_chr <- GWAS_loci$chr[i]
      genei$SNP_pos_hg19 <- GWAS_loci$Pos_hg19[i]
      gene2 <- rbind(gene2,genei)
    }
  }
  
  candidate_gene <- NULL
  for(i in 1:nrow(GWAS_loci)){
    gene2i <- gene2[gene2$SNP==GWAS_loci$SNP[i],]
    gene1i <- gene1[gene1$SNP==GWAS_loci$SNP[i],]
    genei <- rbind(gene2i,gene1i)
    genei <- unique(genei)
    #genei <- genei[match(unique(genei$official_name),genei$official_name),]
    candidate_gene <- rbind(candidate_gene,genei)
  }
  gene <- candidate_gene
  write.table(gene,paste(disease,"candidate_genes.txt",sep = "_"),sep = "\t",quote = F,row.names = F)
  cat("collecting candidate genes ends.\n")
  # DTS
  index <- gene$strand=="+"
  tss_dist<-0
  tss_dist[index] <- abs(gene[index,]$start_hg19-gene[index,]$SNP_pos_hg19)
  tss_dist[!index] <- abs(gene[!index,]$end_hg19-gene[!index,]$SNP_pos_hg19)
  gene$dist <- tss_dist
  # DE
  if(disease!="MG"){
    De <- read.delim(paste("DE/",disease,"_geo2r.tsv",sep = ""),header = T)
    gene$de <- De$P.Value[match(gene$official_name,De$Gene.symbol)]
  }
  # EPI
  en_pro_uniq <- read.table("HiC/brain cell type-specific enhancer-promoters_aggregated.txt",header = T)
  gene$neuron_count <- en_pro_uniq$neuron_count[match(gene$Name,en_pro_uniq$Group.1)]
  gene$oligo_count <- en_pro_uniq$oligo_count[match(gene$Name,en_pro_uniq$Group.1)]
  gene$micro_count <- en_pro_uniq$micro_count[match(gene$Name,en_pro_uniq$Group.1)]
  gene$neuron_count[is.na(gene$neuron_count)] <- 0
  gene$oligo_count[is.na(gene$oligo_count)] <- 0
  gene$micro_count[is.na(gene$micro_count)] <- 0
  cat("collecting omics data ends.\n")
  # Gibbs sampling
  load("eQTL/cis_and_trans_eQTL.RData")
  load(paste(disease,"_GO+BioGRID+HIPP+DD_p11.Rdata",sep = ""))
  candidate_gene <- gene
  risk_gene_eQTL <- merge(eQTL,candidate_gene[,c("Name","SNP")],by=c("Name","SNP"))
  # risk_gene_CMC <- risk_gene_CMC[risk_gene_CMC$FDR<0.05,]
  risk_gene_eQTL <- unique(risk_gene_eQTL)
  candidate_dist <- aggregate(candidate_gene$dist,by = list(candidate_gene$official_name),FUN = min)
  gene2 <- candidate_gene[match(candidate_dist$Group.1,candidate_gene$official_name),]
  gene2$dist <- candidate_dist$x
  for(i in 1:nrow(candidate_dist)){
    snps <- candidate_gene$SNP[candidate_gene$official_name==candidate_dist$Group.1[i]]
    gene2$SNP[i] <- paste(snps,collapse = ",")
  }
  N <- nrow(gene2)
  gene2$DTS <- sapply(gene2$dist,function(x) sum(gene2$dist<=x)/N)
  gene2$neuron_P <- sapply(gene2$neuron_count,function(x) sum(gene2$neuron_count>=x)/N)
  gene2$micro_P <- sapply(gene2$micro_count,function(x) sum(gene2$micro_count>=x)/N)
  gene2$oligo_P <- sapply(gene2$oligo_count,function(x) sum(gene2$oligo_count>=x)/N)
  if(disease!="MG"){
    train <- gene2[!is.na(gene2$de),c("dist","de","neuron_count","micro_count","oligo_count")]
    test <- gene2[is.na(gene2$de),c("dist","de","neuron_count","micro_count","oligo_count")]
    knn.fit <- kknn(de~.,train = train,test = test,k = 5,kernel = "triangular") # use knn to handle NAs
    gene2[is.na(gene2$de),"de"] <- knn.fit$fitted.values
    gene2$DE <- sapply(gene2$de,function(x) sum(gene2$de<=x)/N)
    weight <- apply(gene2[,c("DTS","DE","neuron_P")],1,function(x) sum(-log(x)))
  }else{
    weight <- apply(gene2[,c("DTS","neuron_P")],1,function(x) sum(-log(x)))
  }
  gene2$weight <- weight
  nodes <- colnames(pro_p)
  write.table(nodes,file = paste(disease,"nodes.txt",sep = "_"),row.names = F,quote = F,col.names = F)
  gene2 <- gene2[is.element(gene2$official_name,nodes),]
  pro_p<-pro_p[,(is.element(nodes,unique(gene2$official_name)))]
  n <- nrow(gene2)
  iter <- 0
  dif <- 100
  num0 <- rep(0,n)
  thres <- 1e-4
  gene_name <- gene2$official_name
  ind <- c(1:n)
  set.seed(123)
  sample_size <- 30
  CMC <- risk_gene_eQTL[order(risk_gene_eQTL$Pvalue,decreasing = F)[1:sample_size],] # pick initial genes from eQTL significant genes
  pick <- which(is.element(gene2$Name,CMC$Name),arr.ind = T)
  if(length(pick)<2){
    pick <- sample(ind,sample_size,replace = F)
  }
  freq0 <- rep(0,n)
  cat("sampling starts.\n")
  while(dif > thres){
    iter <- iter+1
    num <- num0
    pp <- pro_p[,is.element(colnames(pro_p),gene_name[pick])]
    pp <- pp[match(gene_name[-pick],nodes),]
    pick <- sample(ind[-pick],sample_size,replace = F,prob = rowSums(pp)*gene2$weight[-pick])
    num[pick] <- num[pick]+1
    if(iter>1) freq0 <- num0/((iter-1)*sample_size)
    freq <- num/(iter*sample_size)
    dif <- norm(freq-freq0,"2")
    num0 <- num
  }
  cat("sampling ends.\n")
  gene2$pp <- freq
  gene2$num <- num
  gene2
}
gene2 <- iGOAT("PD","All_human_genes")#running example
