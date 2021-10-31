library(dplyr)
library(stringr)
library(edgeR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(heatmaps)
library(gplots)
library(dplyr)
library(pheatmap)
library(stringr)
library(grid)
library(DESeq2)

##read all sample raw-count file 

BT_1=read.csv("./BT-1.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BT_2=read.csv("./BT-2.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BT_4=read.csv("./BT-4.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BT_5=read.csv("./BT-5.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BT_6=read.csv("./BT-6.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BT_7=read.csv("./BT-7.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BT_8=read.csv("./BT-8.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BT_9=read.csv("./BT-9.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BT_10=read.csv("./BT-10.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BT_11=read.csv("./BT-11.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BT_13=read.csv("./BT-13.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BTN_1=read.csv("./BTN-1.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BTN_2=read.csv("./BTN-2.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BTN_4=read.csv("./BTN-4.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BTN_5=read.csv("./BTN-5.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BTN_6=read.csv("./BTN-6.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BTN_7=read.csv("./BTN-7.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BTN_8=read.csv("./BTN-8.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BTN_9=read.csv("./BTN-9.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BTN_10=read.csv("./BTN-10.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
BTN_11=read.csv("./BTN-11.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
NNL_297=read.csv("./NNL_297.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]
NNL_143=read.csv("./NNL_143.out.tab",sep="",skip = 4,header=F,row.names = 1)[1]

all.equal(rownames(BT_1),rownames(BT_5),rownames(BT_6),rownames(BT_13),rownames(BT_10),rownames(NNL_297))


expr<-cbind(BT_1=BT_1$V2,BT_2=BT_2$V2,BT_4=BT_4$V2,BT_5=BT_5$V2,BT_6=BT_6$V2,BT_7=BT_7$V2,BT_8=BT_8$V2,BT_9=BT_9$V2,BT_10=BT_10$V2,BT_11=BT_11$V2, BT_13=BT_13$V2,BTN_1=BTN_1$V2,BTN_2=BTN_2$V2,BTN_4=BTN_4$V2,BTN_5=BTN_5$V2,BTN_6=BTN_6$V2,BTN_7=BTN_7$V2,BTN_8=BTN_8$V2,BTN_9=BTN_9$V2,BTN_10=BTN_10$V2,BTN_11=BTN_11$V2,NNL_143=NNL_143$V2,NNL_297=NNL_297$V2)
rownames(expr)<-rownames(BT_1)
#dim(expr) 62446 23

#filter out the lower expression gene
expr=expr[rowSums(expr)>1,]
#dim(expr) 46276 23

#rename the normal patients
row.names(expr)<-str_sub(rownames(expr),1,15)

#add the SYMBOL ID
expr<-as.data.frame(expr)
expr$ENSEMBL<-row.names(expr)
entrezIDs <- bitr(geneID = rownames(expr),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb ="org.Hs.eg.db" ,drop = T) 
expr<-inner_join(expr,entrezIDs,by="ENSEMBL")

#for those symbols have more than one row, aggregate them and take their mean value
expr<-aggregate(x=expr[,1:(ncol(expr)-2)],by=list(expr$SYMBOL),mean,na.rm=TRUE) #2
rownames(expr)<-expr$Group.1
expr<-expr %>% dplyr::select(-Group.1)

expr<-expr %>% mutate(mean_NNL = rowMeans(dplyr::select(., starts_with("NNL")), na.rm = TRUE)) #nrow=874
expr<-expr %>% mutate(mean_BT = rowMeans(dplyr::select(., starts_with("BT_")), na.rm = TRUE)) 
expr<-expr %>% mutate(mean_BTN = rowMeans(dplyr::select(., starts_with("BTN_")), na.rm = TRUE))

expr<-expr %>% filter(mean_BT> mean_NNL,mean_BTN> mean_NNL) 
#dim(expr) 12298

BT_BTN<-dplyr::select(expr,starts_with("BT")|starts_with("BTN"))
BT_NNL<-dplyr::select(expr,starts_with("BT_")|starts_with("NNL"))
BTN_NNL<-dplyr::select(expr,starts_with("BTN")|starts_with("NNL"))
dat<-list(BT_NNL=BT_NNL,BTN_NNL=BTN_NNL,BT_BTN=BT_BTN)
dat<-list(BT_BTN_filter=as.data.frame(Heatmap_input))



###performing the deseq2 step
for (i in names(dat)){dat[[i]]
  condition<-as.factor(word(names(dat[[i]]),1,sep="[_]"))
  print(condition)
  colData <- data.frame(row.names=colnames(dat[[i]]), condition)
  dds <- DESeqDataSetFromMatrix(countData = round(dat[[i]]), colData = colData, design= ~ condition) 
  dds <- DESeq(dds)
  #normlzd_dds <- counts(dds,normalized=TRUE)
  #write.csv(normlzd_dds,"./normlzd_dds.csv")
  res <- results(dds, contrast=c("condition",as.vector(unique(condition))))
  #print(head(res))
  resOrdered<-res[order(res$pvalue),]
  print(head( resOrdered))
  DGE<-as.data.frame(resOrdered)
  logFC_cutoff <-1
  #logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
  k1 = (DGE$pvalue < 0.05)&(DGE$log2FoldChange < - logFC_cutoff)&(DGE$padj<0.05)
  k2 = (DGE$pvalue < 0.05)&(DGE$log2FoldChange > logFC_cutoff)&(DGE$padj<0.05)
  #DGE<-DGE %>% mutate(change=case_when(pvalue<0.05 & log2FoldChange> logFC_cutoff ~ "UP",pvalue<0.05 & log2FoldChange< -logFC_cutoff ~"DOWN",TRUE ~ "NOT"))
  DGE$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
  DGE=DGE[order(DGE$log2FoldChange,decreasing = T),]
  save_dir<-paste0(names(dat[i]),"/")
  dir.create(save_dir)
  write.csv(DGE,paste0(save_dir,paste(levels(condition),collapse = "-"),".csv"),quote = F) #DGE up and down genes are uncertain secreted protein,and this is a protein
  annotation_col = data.frame(TissueType = factor(condition)) 
  print(annotation_col$TissueType)
  rownames(annotation_col) = names(dat[[i]])
  #factor<-read.csv(paste0(gsub("factor",i,"C:/Users/37140/Downloads/factor"),".","tsv"),sep="")
  DGE_factor<-(DGE[rownames(DGE) %in% intersect(rownames(DGE),secretome$Gene),]) #until now, the gene are limited to which produce secreted protein
  head(DGE_factor)
  genes<-row.names(DGE_factor[!(DGE_factor$change=="NOT") & (!is.na(DGE_factor$change)),])
  head(annotation_col)
  #group1<-as.character(annotation_col$TissueType[1])
  #group2<-as.character(annotation_col$TissueType[2])
  annota_colors<-list(TissueType=c("BT"="firebrick","BTN"="pink3","NNL"="white"))
  #pheat<-pheatmap(log2(expr[rownames(expr) %in% genes,]+1)[,c(rownames(colData))],scale = "row", clustering_distance_rows = "correlation",border=FALSE,main = toupper(i),color = colorRampPalette(c("navy", "white", "firebrick3"))(100),annotation_col=annotation_col,annotation_colors=annota_colors)

}

#the function as saving the pheatmap pdf figures
save_pheatmap_pdf <- function(x, filename,width=10,height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

annota_colors<-list(TissueType=c("BT"="firebrick","BTN"="blue3"))
unfavo_exp<-log2(expr[rownames(expr) %in% factor(unfavo_sig),]+1)
rownames(unfavo_exp)<-factor(rownames(unfavo_exp),levels=unfavo_sig)
pheat1<-pheatmap( unfavo_exp[,c(rownames(colData))],scale = "row", clustering_distance_rows = "correlation",border=FALSE,main = toupper("SECRETOME"),annotation_col=annotation_col,annotation_colors=annota_colors, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), fontsize = 15)
save_pheatmap_pdf(pheat1, paste0("unfavo",".pdf"))

#using the normal patients as a control group, compared the relative abundant between tumor and paracancerous tissues 
#keep the genes in secretome,expr has been filtered 
expr.1<-as.data.frame(expr[rownames(expr) %in% unique(secretome$Gene),] )

#dim(expr.1) 571

#filter out the gene without expressing in normal patients 
expr.1<-expr.1[!(expr.1$mean_NNL==0),] #nrow=447,delete about 3 genes

expr_1 <- expr.1 %>% mutate(
  BT_1vNNL = log(BT_1/mean_NNL,2),
  BT_2vNNL = log(BT_2/mean_NNL,2),
  BT_4vNNL = log(BT_4/mean_NNL,2),
  BT_5vNNL = log(BT_5/mean_NNL,2),
  BT_6vNNL = log(BT_6/mean_NNL,2),
  BT_7vNNL = log(BT_7/mean_NNL,2),
  BT_8vNNL = log(BT_8/mean_NNL,2),
  BT_9vNNL = log(BT_9/mean_NNL,2),
  BT_10vNNL =log(BT_10/mean_NNL,2),
  BT_11vNNL = log(BT_11/mean_NNL,2),
  BT_13vNNL = log(BT_13/mean_NNL,2),
  BTN_1vNNL = log(BTN_1/mean_NNL,2),
  BTN_2vNNL = log(BTN_2/mean_NNL,2),
  BTN_4vNNL = log(BTN_4/mean_NNL,2),
  BTN_5vNNL = log(BTN_5/mean_NNL,2),
  BTN_6vNNL = log(BTN_6/mean_NNL,2),
  BTN_7vNNL = log(BTN_7/mean_NNL,2),
  BTN_8vNNL = log(BTN_8/mean_NNL,2),
  BTN_9vNNL = log(BTN_9/mean_NNL,2),
  BTN_10vNNL = log(BTN_10/mean_NNL,2),
  BTN_11vNNL = log(BTN_11/mean_NNL,2)
) 

Heatmap_input <- expr_1 %>% dplyr::select(ends_with("vNNL")) 
Heatmap_input<-as.matrix(Heatmap_input)

#use the Heatmap_input to do the methods1 gene differential analysis
gene<-rownames(BT_BTN_filter[!(BT_BTN_filter$change=="NOT"),])

Heatmap_input_1<-Heatmap_input[rownames(Heatmap_input) %in% sig,] #nrow , The label gene symbolizes the differentially expressed genes in tumor and paracancerous tissues(BT_BTN_RPA)
write.csv(Heatmap_input_1,"./Heatmap_input_1.csv",quote=F)
#Heatmap_input_1<-read.csv("Heatmap_input_1.csv",row.names = 1)




#Derive the colour palette
ColSideColors= c("firebrick","blue3")[as.factor(c(rep("BT",11),rep("BTN",10)))]


#Derive the colour palette
my_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
png("2021_heatmap.png",width=700,height=700,res = 100)
par(mfrow=c(8,8), mar=c(0,1,1,1))
heatmap.2(Heatmap_input_1,main = "",
          notecol="black",
          density.info="density",
          trace="none",
          col=my_palette,
          symbreaks=TRUE,
          dendrogram="both",
          Colv=TRUE,
          Rowv=TRUE,
          scale='none',
          ColSideColors = ColSideColors,
          keysize = 1,
          key.xtickfun = NULL,
          key.ytickfun = NULL,
          key.xlab = "Log2(FC)",
          #xlab = "Subtype",
          ylab = "",
          cexCol = 1,
          cexRow =1,
          symkey = FALSE,
          #font <- par("font.axis"),
          #scale = c("row"),
          labRow = rownames(Heatmap_input_1)
)
dev.off()
#Legend_sort <- Heatmap_input[order(Heatmap_input_1$Group),]

#Add legend to the heatmap at correct posotion next to colour key          
legend(y=1.23, x=.008, xpd=TRUE,      # location of the legend on the heatmap plot
       col = col1,  
       lty= 1,             
       lwd = 5,
       cex = .4)

dev.off()


  
  
