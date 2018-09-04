source("http://bioconductor.org/workflows.R")
workflowInstall("annotation")

source("https://bioconductor.org/biocLite.R")
biocLite(c("AnnotationHub", "Homo.sapiens",
           "Organism.dplyr",
           "TxDb.Hsapiens.UCSC.hg19.knownGene",
           "TxDb.Hsapiens.UCSC.hg38.knownGene",
           "BSgenome.Hsapiens.UCSC.hg19", "biomaRt",
           "TxDb.Athaliana.BioMart.plantsmart22"))
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tidyr)
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

geneRanges <- function(db, column="ENTREZID")
{
  g <- genes(db, columns=column)
  col <- mcols(g)[[column]]
  genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  genes
}

gns = geneRanges(Homo.sapiens,column = 'SYMBOL')
gns.Ensembl = geneRanges(Homo.sapiens, column = 'ENTREZID')
promoters_txdb <- promoters(gns)
promoters_txdb_ENTREZID <- promoters(gns.Ensembl)


#########################
DMP.GEM.VD_w <- read.table('Result_Emodel_EUR.txt',sep = '\t',header = T)
DMP.GEM.VD_w <- DMP.GEM.VD_w[which(substring(DMP.GEM.VD_w$cpg,4,4) != 'X'),]
DMP.GEM.VD_w <- DMP.GEM.VD_w[which(substring(DMP.GEM.VD_w$cpg,4,4) != 'Y'),]

#VD_w
library(stringr)
temp<- DMP.GEM.VD_w[1:100,]
temp <- data.frame(temp, str_split_fixed(temp$cpg, ":", 2))
colnames(temp) <- c('cpg', 'beta', 'stats', 'pvals', 'FDR', 'chr', 'pos')

a<-data.frame(temp$chr,temp$pos,temp$pos,temp$pvals)
colnames(a)<-c('chr','start','end','pvalue')

a$promoter <- paste((splitByOverlap(promoters_txdb, as(a,'GRanges'), "SYMBOL")),collapse="/")
a$gene <- paste((splitByOverlap(gns, as(a,'GRanges'), "SYMBOL")),collapse="/")

t<-nearest(as(a,'GRanges'),gns,ignore.strand=T)
nearest <- gns[t]
a$nearest <- nearest$SYMBOL
a$nearest <- ifelse(a$gene == '' & a$promoter == '', a$nearest, '')

a <- separate_rows(a, gene)
a <- separate_rows(a, promoter)
a$genelist <- ifelse(a$promoter == '', a$gene, a$promoter)
a$genelist <- ifelse(a$gene == '' & a$promoter == '', a$nearest, a$genelist)
a <- a[which(a$genelist != 'NA'),]
genelist <- as.character(a[!duplicated(a$genelist),]$genelist)
genelist_DSS_VD_w <- genelist
write.csv(a,'GeneList_DSS_VD_w.csv')


###################################################################################################
#Manhatten plot
###################################################################################################
#install.packages("QCEWAS")
library(QCEWAS)
trace("EWAS_plots",edit=TRUE)

# split cpg into chr and pos
temp <- read.table('/mnt/DATA/meQTL/Truseq/VD_w.txt',sep = '\t',header = T)
#temp <- read.table('Result_Emodel_EUR.txt',sep = '\t',header = T)
temp2 <- data.frame(temp, str_split_fixed(temp$cpg, ":", 2))
colnames(temp2) <- c('cpg', 'beta', 'stats', 'P_VAL', 'FDR', 'CHR', 'POS')
temp2$CHR <- substr(temp2$CHR,4,length(temp2$CHR))
temp2$CHR <- ifelse(temp2$CHR == 'X', 23,
                           ifelse(temp2$CHR == 'Y', 24, temp2$CHR))

temp2$CHR <- as.numeric(temp2$CHR)
temp2$POS <- as.numeric(as.character(temp2$POS))

# omit NAs
temp2 <- na.omit(temp2)

# Get rid of Y chromosomes
temp2 <- temp2[which(temp2$CHR <= 23),]

keep <- c('P_VAL', 'CHR', 'POS')
temp2 <- temp2[keep]

EWAS_plots(temp2,
           plot_QQ = FALSE,
           plot_Man = TRUE,
           plot_cutoff_p = 1,
           plot_QQ_bands = FALSE,
           high_quality_plots = FALSE,
           save_name = "VDUse_Manhattan")


#########################
# NEW PATHWAYS
#########################
source("https://bioconductor.org/biocLite.R")
biocLite("seq2pathway")

temp2 <- data.frame(temp, str_split_fixed(temp$cpg, ":", 2))
colnames(temp2) <- c('cpg', 'beta', 'stats', 'P_VAL', 'FDR', 'CHR', 'POS')
temp2$CHR <- substr(temp2$CHR,4,length(temp2$CHR))
temp2$CHR <- as.numeric(temp2$CHR)
temp2$POS <- as.numeric(as.character(temp2$POS))

a<-data.frame(temp2$CHR, temp2$POS, temp2$POS, temp2$P_VAL, temp2$FDR)
colnames(a)<-c('chr', 'start', 'end', 'pval', 'fdr')

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
d <- c$gene2pathway_result.FET$GO_CC
d <- c$gene2pathway_result.FET$GO_MF
write.csv(d,'new pathway.csv')
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
write.csv(a,'genelist_cis.csv')

gns_symbol = as.character(gns[!duplicated(gns$SYMBOL),]$SYMBOL)
gns_symbol = na.omit(gns_symbol)

DE.genes <- as.integer(gns_symbol %in% genelist[1:200])
names(DE.genes) <- gns_symbol
pwf <- nullp(DE.genes,"hg19","geneSymbol")

pvals <- goseq(pwf,'hg19','geneSymbol')
over_represented_adjPvalue <- p.adjust(pvals$over_represented_pvalue,method="BH")
pvals <- cbind(pvals,over_represented_adjPvalue)
add.cis.path <- pvals



