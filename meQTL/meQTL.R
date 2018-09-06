library(MatrixEQTL)
require(GEM)
setwd('/home/jcastle/VitaminUse_Study/meQTL')
DATADIR <- "../data/meQTL_norm_and_adjnorm_samples_with_geno/bc_only_loci/EUR/"
RESULTSDIR <- "normal_adjnormal_results/EUR/"
LOADDATA <- FALSE


##### USE Matrix eQTL TO SEPARATE CIS AND TRANS PAIRS #####
# Settings
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(DATADIR, "ImputationResults_all_loci_maf_geno_hwe_EURandAFR_1000G_finalAD.tab.txt", sep = "")

# Gene expression file name
expression_file_name = paste(DATADIR, "meth.txt", sep = "")

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(DATADIR, "cov.txt", sep = "")

# Output file name
output_file_name_cis = paste(RESULTSDIR, "cis_all_loci_mafcut.txt", sep='');
output_file_name_tra = paste(RESULTSDIR, "trans_all_loci_mafcut.txt", sep='');

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 10 * 0.05 / 3335526 / 1237230;
pvOutputThreshold_tra = 1 * 0.05 / 3335526 / 1237230 / 1000000000;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Load SNP and gene locations
snpspos = read.table("../data/SNPloc.txt", header = TRUE, stringsAsFactors = FALSE);
genepos = read.table("../data/methloc.txt", header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

##### SAVE RESULTS #####
save(me, file= paste(RESULTSDIR, "Matrix_eQTL_main.RDATA", sep = ""))

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

##### LOAD RESULTS #####
load( paste(RESULTSDIR, "Matrix_eQTL_main.RDATA", sep = "") )

##### RUN cis-meQTL ANALYSIS USING GEM #####
#covariate_file = "../data/meQTL_norm_and_tumor_samples_with_geno/cov.txt"
#methylation_file = "../data/meQTL_norm_and_tumor_samples_with_geno/meth.txt"
#snp_file = '../data/meQTL_norm_and_tumor_samples_with_geno/ImputationResults_EURandAFR_1000G_finalAD.tab.txt'
#Gmodel_pv = 2e-6;
#output_file = "Result_Gmodel.txt"
#GEM_Gmodel(snp_file, covariate_file, methylation_file, Gmodel_pv, output_file)

###############################################################################################################

source("https://bioconductor.org/biocLite.R")
biocLite(c("AnnotationHub", "Homo.sapiens",
           "Organism.dplyr",
           "TxDb.Hsapiens.UCSC.hg19.knownGene",
           "TxDb.Hsapiens.UCSC.hg38.knownGene",
           "BSgenome.Hsapiens.UCSC.hg19", "biomaRt",
           "TxDb.Athaliana.BioMart.plantsmart22"))
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(Homo.sapiens)
library(goseq)

geneRanges <- function(db, column="SYMBOL")
{
  g <- genes(db, columns=column)
  col <- mcols(g)[[column]]
  genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  genes
}

splitByOverlap <- function(query, subject, column="ENTREZID", ...)
{
  olaps <- findOverlaps(query, subject, maxgap = 0)
  f1 <- factor(subjectHits(olaps),
               levels=seq_len(subjectLength(olaps)))
  splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
}

gns = geneRanges(Homo.sapiens,column = 'SYMBOL')

geneRanges <- function(db, column="ENTREZID")
{
  g <- genes(db, columns=column)
  col <- mcols(g)[[column]]
  genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  genes
}

gns.Ensembl = geneRanges(Homo.sapiens, column = 'ENTREZID')
promoters_txdb <- promoters(gns)
promoters_txdb_ENTREZID <- promoters(gns.Ensembl)

##################################################################################
# MANHATTAN PLOT
##################################################################################
source("https://bioconductor.org/biocLite.R")
biocLite("seq2pathway")
library(stringr)
library(QCEWAS)

##### CIS #####
me.add.cis = me$cis$eqtls
a <- str_split_fixed(me.add.cis$gene, "-", 2)
me.add.cis$chr <- substring(a[,1],4,length(a[,1]))
me.add.cis$POS <- a[,2]
colnames(me.add.cis) <- c("snps" ,'gene',"stats" ,"P_VAL" , "fdr", 'beta',"CHR", 'POS')
me.add.cis$CHR <- as.numeric(me.add.cis$CHR)
me.add.cis$POS <- as.numeric(me.add.cis$POS)

trace("EWAS_plots",edit=TRUE)
EWAS_plots(me.add.cis,
           plot_QQ = FALSE,
           plot_Man = TRUE,
           plot_cutoff_p = 1,
           plot_QQ_bands = FALSE,
           high_quality_plots = FALSE,
           save_name = paste(RESULTSDIR, "plots/Cis", sep="")
           )

##### TRANS #####
me.add.trans = me$trans$eqtls
a <- str_split_fixed(me.add.trans$gene, "-", 2)
me.add.trans$chr <- substring(a[,1],4,length(a[,1]))
me.add.trans$POS <- a[,2]
colnames(me.add.trans) <- c("snps" ,'gene',"stats" ,"P_VAL" , "fdr", 'beta',"CHR", 'POS')
me.add.trans$CHR <- as.numeric(me.add.trans$CHR)
me.add.trans$POS <- as.numeric(me.add.trans$POS)

trace("EWAS_plots",edit=TRUE)
EWAS_plots(me.add.trans,
           plot_QQ = FALSE,
           plot_Man = TRUE,
           plot_cutoff_p = 1,
           plot_QQ_bands = FALSE,
           high_quality_plots = FALSE,
           save_name = paste(RESULTSDIR, "plots/Trans", sep = "")
           )

##################################################################################
# DESCRIPTIVE DISTRIBUTION
##################################################################################
detach("package:dplyr", unload = TRUE)
library(plyr)

##### CIS #####
me.add.cis = me$cis$eqtls
png( paste(RESULTSDIR, "plots/CisDistn_CpGs_per_SNP_pval0p05.png", sep = "") )
hist(count(me.add.cis$snps)$freq,breaks = 100,main = 'Distribution of CpGs per SNPs (cis) pval < 0.05',xlab = '# of CpGs')
dev.off()

png( paste(RESULTSDIR, "plots/CisDistn_CpGs_per_SNP_FDR0p05.png", sep =="") )
me.add.cis = me.add.cis[which(me.add.cis$FDR < 0.05),]
hist(count(me.add.cis$snps)$freq,breaks = 100,main = 'Distribution of CpGs per SNPs (cis) fdr < 0.05',xlab = '# of CpGs')
dev.off()

summary(count(me.add.cis$snps)$freq)

##### TRANS #####
me.add.trans = me$trans$eqtls
png( paste(RESULTSDIR,"plots/TransDistn_CpGs_per_SNP_pval0p05.png", sep = "") )
hist(count(me.add.trans$snps)$freq,breaks = 100,main = 'Distribution of CpGs per SNPs (trans) pval < 0.05',xlab = '# of CpGs')
dev.off()

png( paste(RESULTSDIR, "plots/TransDistn_CpGs_per_SNP_FDR0p05.png", sep = "") )
me.add.trans = me.add.trans[which(me.add.trans$FDR < 0.05),]
hist(count(me.add.trans$snps)$freq,breaks = 100,main = 'Distribution of CpGs per SNPs (trans) fdr < 0.05',xlab = '# of CpGs')
dev.off()

summary(count(me.add.trans$snps)$freq)


##################################################################################
# DISTANCE DENSITY PLOT
##################################################################################
library(latticeExtra)

##### CIS #####
me.add.cis = merge(merge(me.add.cis, snpspos, by.x = 'snps',by.y = 'snp'),genepos[,c('CpG','s1')], by.x = 'gene', by.y = 'CpG')
me.add.cis$dist = me.add.cis$pos - me.add.cis$s1
x = densityplot(me.add.cis[which(me.add.cis$dist != 0),]$dist,type = NULL,col = 'blue',lwd=5)
y = histogram(me.add.cis[which(me.add.cis$dist != 0),]$dist,main = 'Distance between CpG and SNPs (Cis)',breaks = 100, xlab = 'Distance',ylab = 'Count of SNPs', col = 'red')
png( paste(RESULTSDIR, "plots/CisDist_CpG_SNP_distance.png", sep = "") )
doubleYScale(y, x, add.ylab2 = TRUE)
dev.off()

png( paste(RESULTSDIR, "plots/Cis_beta_vs_dist.png", sep = "") )
plot(me.add.cis$dist,me.add.cis$beta)
dev.off()

png( paste(RESULTSDIR, "plots/Cis_FDR_vs_dist.png", sep = "") )
plot(me.add.cis$dist,me.add.cis$FDR)
dev.off()

me.add.cis.50 <- me.add.cis[which(abs(me.add.cis$dist) <= 50),]
x = densityplot(me.add.cis.50$dist,type = NULL)
y = histogram(me.add.cis.50$dist,main = 'Distance between CpG and SNPs (Cis within 50bp)',breaks = 50, xlab = 'Distance')
png( paste(RESULTSDIR, "plots/CisDist_CpG_SNP_distance_bp50.png", sep = "") )
doubleYScale(y, x, add.ylab2 = TRUE)
dev.off()

me.add.cis.2000 <- me.add.cis[which(abs(me.add.cis$dist) <= 2000),]
x = densityplot(me.add.cis.2000$dist,type = NULL)
y = histogram(me.add.cis.2000$dist,main = 'Distance between CpG and SNPs (Cis within 2000bp)',breaks = 100, xlab = 'Distance')
png( paste(RESULTSDIR, "plots/CisDist_CpG_SNP_distance_bp2000.png", sep = "") )
doubleYScale(y, x, add.ylab2 = TRUE)
dev.off()

##### TRANS #####
me.add.trans = merge(merge(me.add.trans, snpspos, by.x = 'snps',by.y = 'snp'),genepos[,c('CpG','s1')], by.x = 'gene', by.y = 'CpG')
me.add.trans$dist = me.add.trans$pos - me.add.trans$s1
x = densityplot(me.add.trans[which(me.add.trans$dist != 0),]$dist,type = NULL,col = 'blue',lwd=5)
y = histogram(me.add.trans[which(me.add.trans$dist != 0),]$dist,main = 'Distance between CpG and SNPs (Trans)',breaks = 100, xlab = 'Distance',ylab = 'Count of SNPs', col = 'red')
png( paste(RESULTSDIR, "plots/TransDist_CpG_SNP_distance.png", sep = "") )
doubleYScale(y, x, add.ylab2 = TRUE)
dev.off()

png( paste(RESULTSDIR, "plots/Trans_beta_vs_dist.png", sep = "") )
plot(me.add.trans$dist,me.add.trans$beta)
dev.off()

png( paste(RESULTSDIR, "plots/Trans_FDR_vs_dist.png", sep = "") )
plot(me.add.trans$dist,me.add.trans$FDR)
dev.off()

me.add.trans.50 <- me.add.trans[which(abs(me.add.trans$dist) <= 50),]
x = densityplot(me.add.trans.50$dist,type = NULL)
y = histogram(me.add.trans.50$dist,main = 'Distance between CpG and SNPs (Trans within 50bp)',breaks = 50, xlab = 'Distance')
png( paste(RESULTSDIR, "plots/TransDist_CpG_SNP_distance_bp50.png", sep = "") )
doubleYScale(y, x, add.ylab2 = TRUE)
dev.off()

me.add.trans.2000 <- me.add.trans[which(abs(me.add.trans$dist) <= 2000),]
x = densityplot(me.add.trans.2000$dist,type = NULL)
y = histogram(me.add.trans.2000$dist,main = 'Distance between CpG and SNPs (Trans within 2000bp)',breaks = 100, xlab = 'Distance')
png( paste(RESULTSDIR, "plots/TransDist_CpG_SNP_distance_bp2000.png", sep = "") )
doubleYScale(y, x, add.ylab2 = TRUE)
dev.off()


##################################################################################
# HEATMAP FOR SNPS AND CPGS 
##################################################################################
library(ggplot2)
library(dplyr) # easier data wrangling 
library(lubridate) # for easy date manipulation
library(ggExtra) # because remembering ggplot theme options is beyond me
library(tidyr) 
library(viridis)
library(reshape2)

##### CIS #####
df <- me.add.cis
df$chrnum <- as.numeric(substr(df$chr,4,length(df$chr)))
df <- df[order(df$chrnum),]
df$pos <- as.numeric(df$pos)
df$s1 <- as.numeric(df$s1)

df <- df[order(df$chrnum,df$s1,df$pos),]
df$snppos <- paste(df$chrnum,df$pos,sep = ':')
df$cpgpos <- paste(df$chrnum,df$s1,sep = ':')
df$snpchr <- df$chr
df$cpgchr <- df$chr

# Plotting starts here
df1 <- df[which(df$chrnum == 6),]
df1 <- df1[,c(8,9,4)]
df1_mat <- acast(df1, pos~s1)
df1_mat <- df1_mat[order(as.numeric(rownames(df1_mat))),order(as.numeric(colnames(df1_mat)))]
df1_melt <- melt(df1_mat)
df1_melt$Var1 <- paste('chr6:',df1_melt$Var1, sep = '')
df1_melt$Var2 <- paste('chr6:',df1_melt$Var2, sep = '')
df1_melt$Var1 <- factor(df1_melt$Var1, as.character(unique(df1_melt$Var1)))
df1_melt$Var2 <- factor(df1_melt$Var2, as.character(unique(df1_melt$Var2)))

p <-ggplot(df1_melt,aes(Var2,Var1,fill=value))+
  geom_tile(color= "white",size=0.1) + 
  scale_fill_viridis(name="value",option="C",na.value="white")
#+facet_grid(chrnumsnp ~ chrnumcpg,  scales="free", space="free",labeller = label_both)
p <-p + labs(title= NULL, x="CpG Position", y="Snp Position")
p <-p + 
  #theme(legend.position = "bottom")+
  theme(plot.title=element_text(size = 15))+
  theme(axis.text.y=element_text(size=15)) +
  theme(strip.background = element_rect(colour="black"))+
  theme(plot.title=element_text(hjust=10))+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=10))+
  theme(legend.title=element_text(size=15))+
  theme(legend.text=element_text(size=15))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  removeGrid()#ggExtra

ggsave(paste(RESULTSDIR, 'plots/Cis_chr6.png', sep = ""), device = 'png', width = 28, height = 8, unit = 'in')

##### TRANS #####
df <- me.add.trans
df$chrnum <- as.numeric(substr(df$chr,4,length(df$chr)))
df <- df[order(df$chrnum),]
df$pos <- as.numeric(df$pos)
df$s1 <- as.numeric(df$s1)

df <- df[order(df$chrnum,df$s1,df$pos),]
df$snppos <- paste(df$chrnum,df$pos,sep = ':')
df$cpgpos <- paste(df$chrnum,df$s1,sep = ':')
df$snpchr <- df$chr
df$cpgchr <- df$chr

# Plotting starts here
df1 <- df[which(df$chrnum == 6),]
df1 <- df1[,c(8,9,4)]
df1_mat <- acast(df1, pos~s1)
df1_mat <- df1_mat[order(as.numeric(rownames(df1_mat))),order(as.numeric(colnames(df1_mat)))]
df1_melt <- melt(df1_mat)
df1_melt$Var1 <- paste('chr6:',df1_melt$Var1, sep = '')
df1_melt$Var2 <- paste('chr6:',df1_melt$Var2, sep = '')

df1_melt$Var1 <- factor(df1_melt$Var1, as.character(unique(df1_melt$Var1)))
df1_melt$Var2 <- factor(df1_melt$Var2, as.character(unique(df1_melt$Var2)))

p <-ggplot(df1_melt,aes(Var2,Var1,fill=value))+
  geom_tile(color= "white",size=0.1) + 
  scale_fill_viridis(name="value",option="C",na.value="white")
#+facet_grid(chrnumsnp ~ chrnumcpg,  scales="free", space="free",labeller = label_both)
p <-p + labs(title= NULL, x="CpG Position", y="Snp Position")
p <-p + 
  #theme(legend.position = "bottom")+
  theme(plot.title=element_text(size = 15))+
  theme(axis.text.y=element_text(size=15)) +
  theme(strip.background = element_rect(colour="black"))+k
  theme(plot.title=element_text(hjust=10))+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=10))+
  theme(legend.title=element_text(size=15))+
  theme(legend.text=element_text(size=15))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  removeGrid()#ggExtra

ggsave(paste(RESULTSDIR, 'plots/Trans_chr6.png', sep = ""), device = 'png', width = 28, height = 8, unit = 'in')

##################################################################################
# GENE ANNOTATION AND PATHWAYS
##################################################################################
source("https://bioconductor.org/biocLite.R")
biocLite("seq2pathway")
library(seq2pathway)

##### FREE-UP SOME MEMORY #####
if( exists("df") ) { rm(df) }
if( exists("df1") ) { rm(df1) }
if( exists("df1_melt") ) { rm(df1_melt) }
if( exists("p") ) { rm(p) }
if( exists("genepos") ) { rm(genepos) }
if( exists("snpspos") ) { rm(snpspos) }
gc()

##### CIS #####
me.add.cis = merge(merge(me.add.cis, snpspos, by.x = 'snps',by.y = 'snp'),genepos[,c('CpG','s1')], by.x = 'gene', by.y = 'CpG')
me.add.cis$dist = me.add.cis$pos - me.add.cis$s1
me.add.cis <- me.add.cis[order(me.add.cis$FDR),]
a<-data.frame(me.add.cis$chr,me.add.cis$s1,me.add.cis$s1,me.add.cis$snps,me.add.cis$pvalue,me.add.cis$FDR)
colnames(a)<-c('chr','start','end','snps','pval','fdr')

b <- runseq2gene(as(a,'GRanges'),
                 search_radius=150000, promoter_radius=200, promoter_radius2=100,
                 genome=c("hg19"), adjacent=FALSE, SNP=FALSE,
                 PromoterStop=FALSE,NearestTwoDirection=TRUE,UTR3=FALSE)

c <- runseq2pathway(as(a,'GRanges'),
                    search_radius=150000, promoter_radius=200, promoter_radius2=100,
                    genome=c('hg19'), adjacent=FALSE, SNP= FALSE,
                    PromoterStop=FALSE, NearestTwoDirection=TRUE,UTR3=FALSE,
                    DataBase=c("GOterm"), FAIMETest=FALSE, FisherTest=TRUE,
                    collapsemethod=c("Average"),
                    alpha=5, logCheck=FALSE, B=100, na.rm=FALSE, min_Intersect_Count=5)

d <- c$gene2pathway_result.FET$GO_BP
write.csv(d, paste(RESULTSDIR, "pathways_results/Cis_BP_pathway.csv", sep = ""))
d <- c$gene2pathway_result.FET$GO_CC
write.csv(d, paste(RESULTSDIR, "pathways_results/Cis_CC_pathway.csv", sep = ""))
d <- c$gene2pathway_result.FET$GO_MF
write.csv(d, paste(RESULTSDIR,"pathways_results/Cis_MF_pathway.csv", sep = ""))

a$promoter <- paste((splitByOverlap(promoters_txdb, as(a,'GRanges'), "SYMBOL")),collapse="/")
a$promoter_ENTREZID <- paste((splitByOverlap(promoters_txdb_ENTREZID, as(a,'GRanges'), "ENTREZID")),collapse="/")
a$gene <- paste((splitByOverlap(gns, as(a,'GRanges'), "SYMBOL")),collapse="/")
a$gene_ENTREZID <- paste((splitByOverlap(gns.Ensembl, as(a,'GRanges'), "ENTREZID")),collapse="/")

t<-nearest(as(a,'GRanges'),gns,ignore.strand=T)
nearest <- gns[t]
a$nearest <- nearest$SYMBOL
a$nearest <- ifelse(a$gene == '' & a$promoter == '', a$nearest, '')

t<-nearest(as(a,'GRanges'),gns.Ensembl,ignore.strand=T)
nearest <- gns.Ensembl[t]
a$nearest_ENTREZID <- nearest$ENTREZID
a$nearest_ENTREZID <- ifelse(a$gene == '' & a$promoter == '', a$nearest, '')

a <- separate_rows(a, gene)
a <- separate_rows(a, promoter)
a$genelist <- ifelse(a$promoter == '', a$gene, a$promoter)
a$genelist <- ifelse(a$gene == '' & a$promoter == '', a$nearest, a$genelist)
genelist <- as.character(a[!duplicated(a$genelist),]$genelist)
#change this line
genelist_cis <- genelist
write.csv(a, paste(RESULTSDIR, 'pathways_results/Cis_genelist.csv', sep = ""))

gns_symbol = as.character(gns[!duplicated(gns$SYMBOL),]$SYMBOL)
gns_symbol = na.omit(gns_symbol)

DE.genes <- as.integer(gns_symbol %in% genelist[1:200])
names(DE.genes) <- gns_symbol
pwf <- nullp(DE.genes,"hg19","geneSymbol")

pvals <- goseq(pwf,'hg19','geneSymbol')
over_represented_adjPvalue <- p.adjust(pvals$over_represented_pvalue,method="BH")
pvals <- cbind(pvals,over_represented_adjPvalue)
add.cis.path <- pvals

##### TRANS #####
me.add.trans = merge(merge(me.add.trans, snpspos, by.x = 'snps',by.y = 'snp'),genepos[,c('CpG','s1')], by.x = 'gene', by.y = 'CpG')
me.add.trans$dist = me.add.trans$pos - me.add.trans$s1
me.add.trans <- me.add.trans[order(me.add.trans$FDR),]
a<-data.frame(me.add.trans$chr,me.add.trans$s1,me.add.trans$s1,me.add.trans$snps,me.add.trans$pvalue,me.add.trans$FDR)
if( exists("me.add.trans") ) { rm(me.add.trans) }
if( exists("b") ) { rm(b) }
if( exists("c") ) { rm(c) }
gc()

b <- runseq2gene(as(a,'GRanges'),
                 search_radius=150000, promoter_radius=200, promoter_radius2=100,
                 genome=c("hg19"), adjacent=FALSE, SNP=FALSE,
                 PromoterStop=FALSE,NearestTwoDirection=TRUE,UTR3=FALSE)

c <- runseq2pathway(as(a,'GRanges'),
                    search_radius=150000, promoter_radius=200, promoter_radius2=100,
                    genome=c('hg19'), adjacent=FALSE, SNP= FALSE,
                    PromoterStop=FALSE, NearestTwoDirection=TRUE,UTR3=FALSE,
                    DataBase=c("GOterm"), FAIMETest=FALSE, FisherTest=TRUE,
                    collapsemethod=c("Average"),
                    alpha=5, logCheck=FALSE, B=100, na.rm=FALSE, min_Intersect_Count=5)

d <- c$gene2pathway_result.FET$GO_BP
write.csv(d, paste(RESULTSDIR, "pathways_results/Trans_BP_pathway.csv", sep = ""))
d <- c$gene2pathway_result.FET$GO_CC
write.csv(d, paste(RESULTSDIR, "pathways_results/Trans_CC_pathway.csv", sep = ""))
d <- c$gene2pathway_result.FET$GO_MF
write.csv(d, paste(RESULTSDIR, "pathways_results/Trans_MF_pathway.csv", sep = ""))

a$promoter <- paste((splitByOverlap(promoters_txdb, as(a,'GRanges'), "SYMBOL")),collapse="/")
a$promoter_ENTREZID <- paste((splitByOverlap(promoters_txdb_ENTREZID, as(a,'GRanges'), "ENTREZID")),collapse="/")
a$gene <- paste((splitByOverlap(gns, as(a,'GRanges'), "SYMBOL")),collapse="/")
a$gene_ENTREZID <- paste((splitByOverlap(gns.Ensembl, as(a,'GRanges'), "ENTREZID")),collapse="/")

t<-nearest(as(a,'GRanges'),gns,ignore.strand=T)
nearest <- gns[t]
a$nearest <- nearest$SYMBOL
a$nearest <- ifelse(a$gene == '' & a$promoter == '', a$nearest, '')

t<-nearest(as(a,'GRanges'),gns.Ensembl,ignore.strand=T)
nearest <- gns.Ensembl[t]
a$nearest_ENTREZID <- nearest$ENTREZID
a$nearest_ENTREZID <- ifelse(a$gene == '' & a$promoter == '', a$nearest, '')

a <- separate_rows(a, gene)
a <- separate_rows(a, promoter)
a$genelist <- ifelse(a$promoter == '', a$gene, a$promoter)
a$genelist <- ifelse(a$gene == '' & a$promoter == '', a$nearest, a$genelist)
genelist <- as.character(a[!duplicated(a$genelist),]$genelist)
#change this line
genelist_trans <- genelist
write.csv(a, paste(RESULTSDIR, 'pathways_results/Trans_genelist.csv', sep = ""))

gns_symbol = as.character(gns[!duplicated(gns$SYMBOL),]$SYMBOL)
gns_symbol = na.omit(gns_symbol)

DE.genes <- as.integer(gns_symbol %in% genelist[1:200])
names(DE.genes) <- gns_symbol
pwf <- nullp(DE.genes,"hg19","geneSymbol")

pvals <- goseq(pwf,'hg19','geneSymbol')
over_represented_adjPvalue <- p.adjust(pvals$over_represented_pvalue,method="BH")
pvals <- cbind(pvals,over_represented_adjPvalue)
add.trans.path <- pvals




##################################################################################
# GENE ANNOTATION FOR SNP AND COMPARE TO EQTL
##################################################################################

##### CIS #####
promoters_txdb <- promoters(gns)
me.add.cis<- me.add.cis[order(me.add.cis$FDR),]
a<-data.frame(me.add.cis$chr,me.add.cis$s1,me.add.cis$s1,me.add.cis$snps)
colnames(a)<-c('chr','start','end','snps')

a$promoter <- paste((splitByOverlap(promoters_txdb, as(a,'GRanges'), "SYMBOL")),collapse="/")
a$gene <- paste((splitByOverlap(gns, as(a,'GRanges'), "SYMBOL")),collapse="/")

t<-nearest(as(a,'GRanges'),gns,ignore.strand=T)
nearest <- gns[t]
a$nearest <- nearest$SYMBOL
a$nearest <- ifelse(a$gene == '' & a$promoter == '', a$nearest, '')

a <- separate_rows(a, gene)
a <- separate_rows(a, promoter)
a$genelist <- ifelse(a$gene == '', a$promoter, a$gene)
a$genelist <- ifelse(a$gene == '' & a$promoter == '', a$nearest, a$genelist)
a <- a[which(a$genelist != 'NA'),]
genelist <- as.character(a[!duplicated(a$genelist),]$genelist)

gns_symbol = as.character(gns[!duplicated(gns$SYMBOL),]$SYMBOL)
gns_symbol = na.omit(gns_symbol)
DE.genes <- as.integer(gns_symbol %in% genelist)
names(DE.genes) <- gns_symbol
pwf <- nullp(DE.genes,"hg19","geneSymbol")

pvals <- goseq(pwf,'hg19','geneSymbol')
over_represented_adjPvalue <- p.adjust(pvals$over_represented_pvalue,method="BH")
pvals <- cbind(pvals,over_represented_adjPvalue)
add.cis.path <- pvals

#####Cis multiple SNPs to CpG #####

##### Cis > 1 #####
a <-plyr::count(me.add.cis$gene)
#change the criteria
a<-a[which(a$freq > 1),]
me.add.cis_multi <- me.add.cis[which(me.add.cis$gene %in% a$x),]
me.add.cis_multi<- me.add.cis_multi[order(me.add.cis_multi$FDR),]
a<- str_split_fixed(me.add.cis_multi$gene, fixed("-"), 2)
a<-data.frame(me.add.cis_multi$chr,me.add.cis_multi$s1,me.add.cis_multi$s1,me.add.cis_multi$snps)
colnames(a)<-c('chr','start','end','snps')

a$promoter <- paste((splitByOverlap(promoters_txdb, as(a,'GRanges'), "SYMBOL")),collapse="/")
a$gene <- paste((splitByOverlap(gns, as(a,'GRanges'), "SYMBOL")),collapse="/")

t<-nearest(as(a,'GRanges'),gns,ignore.strand=T)
nearest <- gns[t]
a$nearest <- nearest$SYMBOL
a$nearest <- ifelse(a$gene == '' & a$promoter == '', a$nearest, '')

a <- separate_rows(a, gene,sep = '/')
a <- separate_rows(a, promoter,sep = '/')
a$genelist <- ifelse(a$gene == '', a$promoter, a$gene)
a$genelist <- ifelse(a$gene == '' & a$promoter == '', a$nearest, a$genelist)
a <- a[which(a$genelist != 'NA'),]
genelist <- as.character(a[!duplicated(a$genelist),]$genelist)

gns_symbol = as.character(gns[!duplicated(gns$SYMBOL),]$SYMBOL)
gns_symbol = na.omit(gns_symbol)

DE.genes <- as.integer(gns_symbol %in% genelist)
names(DE.genes) <- gns_symbol
pwf <- nullp(DE.genes,"hg19","geneSymbol")
pvals <- goseq(pwf,'hg19','geneSymbol')
over_represented_adjPvalue <- p.adjust(pvals$over_represented_pvalue,method="BH")
pvals <- cbind(pvals,over_represented_adjPvalue)
add.cis_1.path <- pvals


#### Cis > 2 #####
a <-plyr::count(me.add.cis$gene)
#change the criteria
a<-a[which(a$freq > 2),]
me.add.cis_multi <- me.add.cis[which(me.add.cis$gene %in% a$x),]
me.add.cis_multi<- me.add.cis_multi[order(me.add.cis_multi$FDR),]
a<- str_split_fixed(me.add.cis_multi$gene, fixed("-"), 2)
a<-data.frame(me.add.cis_multi$chr,me.add.cis_multi$s1,me.add.cis_multi$s1,me.add.cis_multi$snps)
colnames(a)<-c('chr','start','end','snps')

a$promoter <- paste((splitByOverlap(promoters_txdb, as(a,'GRanges'), "SYMBOL")),collapse="/")
a$gene <- paste((splitByOverlap(gns, as(a,'GRanges'), "SYMBOL")),collapse="/")

t<-nearest(as(a,'GRanges'),gns,ignore.strand=T)
nearest <- gns[t]
a$nearest <- nearest$SYMBOL
a$nearest <- ifelse(a$gene == '' & a$promoter == '', a$nearest, '')

a <- separate_rows(a, gene,sep = '/')
a <- separate_rows(a, promoter,sep = '/')
a$genelist <- ifelse(a$gene == '', a$promoter, a$gene)
a$genelist <- ifelse(a$gene == '' & a$promoter == '', a$nearest, a$genelist)
a <- a[which(a$genelist != 'NA'),]
genelist <- as.character(a[!duplicated(a$genelist),]$genelist)

gns_symbol = as.character(gns[!duplicated(gns$SYMBOL),]$SYMBOL)
gns_symbol = na.omit(gns_symbol)

DE.genes <- as.integer(gns_symbol %in% genelist)
names(DE.genes) <- gns_symbol
pwf <- nullp(DE.genes,"hg19","geneSymbol")
pvals <- goseq(pwf,'hg19','geneSymbol')
over_represented_adjPvalue <- p.adjust(pvals$over_represented_pvalue,method="BH")
pvals <- cbind(pvals,over_represented_adjPvalue)
add.cis_2.path <- pvals