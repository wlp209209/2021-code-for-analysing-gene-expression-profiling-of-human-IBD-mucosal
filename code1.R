library(GEOquery)
library(limma)        
library(plyr)
library( hugene10sttranscriptcluster.db)
gset <- getGEO('GSE75214',destdir = ".")
gpl6244 <- getGEO('GPL6244',destdir = ".")
exprSet <- exprs(gset[[1]])
pdate <- pData(gset[[1]])
ids= toTable( hugene10sttranscriptclusterSYMBOL)
exprSet = exprSet[rownames(exprSet)%in%ids$probe_id,]
table(rownames(exprSet)%in%ids$probe_id) 
ids=ids[match(rownames(exprSet),ids$probe_id),]
tmp = by(exprSet,
         ids$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
exprSet=exprSet[rownames(exprSet)%in%probes,]
rownames(exprSet)=ids[which(ids$probe_id%in%a),'symbol']
clinicaldata=pdate[,c('geo_accession','characteristics_ch1.2')]
clinicaldata$group='normal'
clinicaldata[grep(pattern="active",b[,'characteristics_ch1.2']),'group']='active'
colData <- data.frame(row.names=colnames(exprSet), group_list=b$group)
group_list=b$group
design <- model.matrix(~0+factor(group_list,levels = c('active', 'normal')))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
contrast.matrix<-makeContrasts('active - normal',levels = design)
fit <- lmFit(exprSet,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
temOutput=topTable(fit2,coef = 1,n=Inf)
nrDEG=na.omit(temOutput)


library(ggrepel)
library(ggpubr)
library(ggthemes)
library(ggplot2)
kegmt<-read.gmt('~/c2.cp.kegg.v7.1.entrez.gmt')
NFKB=kegmt[which(kegmt$ont=='NF-kappa B signaling pathway'),]
upreg=rownames(nrDEG)[which((nrDEG$adj.P.Val<0.05)&(nrDEG$logFC>1))] 
NKup=upreg[upreg%in%NFKB$gene] 
nrDEG$logP=-log10(nrDEG$adj.P.Val)
nrDEG$Gene=rownames(nrDEG)
nrDEG$group='not-significant'
nrDEG$group[which((nrDEG$adj.P.Val<0.05)&(nrDEG$logFC>1))]='up-regulated'
nrDEG$group[which((nrDEG$adj.P.Val<0.05)&(nrDEG$logFC< -1))]='down-regulated'
nrDEG$label=ifelse(nrDEG$Gene%in%NKup,nrDEG$Gene,'')
ggscatter(nrDEG,x='logFC',y='logP',color = 'group',
          palette = c('green','gray','red'),size = 1)+theme_base()+
  geom_hline(yintercept = 1.30,linetype='dashed')+
  geom_vline(xintercept = -1,linetype='dashed')+
  geom_vline(xintercept = 1,linetype='dashed')+
  geom_text_repel(data = nrDEG, aes(label=label))



library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
deg <- nrDEG[abs(nrDEG$logFC)>=1&nrDEG$P.Value<0.05,]
gene <- rownames(deg)
gene.df <- bitr(gene, fromType ="SYMBOL" ,
                toType = c("ENTREZID" ),
                OrgDb = org.Hs.eg.db)
kegg <- enrichKEGG(gene       = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
kegg=subset(kegg@result,pvalue< 0.05)
kegg$GeneRatio=parse_ratio(kegg$GeneRatio)
kegg=kegg[significance,]
kegg$Description <- factor(kegg$Description, levels=unique(as.character(kegg$Description)) )
ggplot(kk2,aes(x=GeneRatio,y=Description)) + 
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_colour_gradient(low="blue",high="red")+
  labs(
    color=expression(-log[10](P.value)),
    size="Gene number")+
  theme_bw()

gene=bitr(rownames(nrDEG),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene$logFC=nrDEG[which(rownames(nrDEG)%in%gene$SYMBOL),'logFC']
geneList<-gene$logFC 
names(geneList)=gene$SYMBOL 
geneList=sort(geneList,decreasing = T)
gsea<-GSEA(geneList,TERM2GENE = kegmt,pvalueCutoff = 0.5)
gseaplot2(gsea,'NF-kappa B signaling pathway',color="red",pvalue_table = T)

