library(clusterProfiler)
library(DOSE)

library(org.Hs.eg.db)
library(topGO)
library(pathview)
library(ggplot2)
library(readxl)
library(enrichplot)
library(Cairo)
rm(list=ls())


setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI/BP/ODI_Neg_ALL_GO/Enrichment_GO')
gene1<-read.csv('GO_AllLists.csv',header=TRUE)
setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI/BP/ODI_Neg_ALL_GO/Enrichment_heatmap')
gene_ref<-read.csv('HeatmapSelectedGO_ref.csv',header=TRUE)

# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI/BP/ODI_Pos_ALL_GO/Enrichment_GO')
# gene1<-read.csv('GO_AllLists.csv',header=TRUE)
# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan_new/Results/ODI/BP/ODI_Pos_ALL_GO/Enrichment_heatmap')
# gene_ref<-read.csv('HeatmapSelectedGO_ref.csv',header=TRUE)


# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan/Results/NDI/BP/PLS-Endo/Enrichment_GO')
# gene1<-read.csv('GO_AllLists.csv',header=TRUE)
# setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan/Results/NDI/BP/PLS-Endo/Enrichment_heatmap')
# gene_ref<-read.csv('HeatmapSelectedGO_ref.csv',header=TRUE)


#setwd('D:/ZZZ/Manuscripts/Image_Gene/Gradient_test/Lifespan/Results/NDI/BP/Multi_gene_list_PLS_Pos_immune/Enrichment_GO')



#dotplot(gene1,showCategory=10,includeAll=TRUE)
gene1$Log.q.value.<-(-gene1$Log.q.value.)
gene1$Ratio<-gene1$X.GeneInGOAndHitList#/gene1$X.GeneInHitList
gene2<-gene1 %>% filter(Log.q.value.>0)



index_sel<-intersect(gene1$Description,gene2$Description)

index<-NULL
for (i in 1:nrow(gene2))
{
  
  if (!is.na(match(gene2$Description[i],gene_ref$Description)))
  {index<-c(index,i)}
}

gene1<-gene2[index,] 

  
for (i in 1:nrow(gene1))
{
  if (!length(setdiff(gene1$Category[i],'GO Biological Processes')))
  {gene1[i,'Description']<-paste0('GO:',gene1$Description[i])}
  else if (!length(setdiff(gene1$Category[i],'KEGG Pathway')))
  {gene1[i,'Description']<-paste0('KEGG:',gene1$Description[i])}
  else if (!length(setdiff(gene1$Category[i],'Hallmark Gene Sets')))
  {gene1[i,'Description']<-paste0('H:',gene1$Description[i])}
  
}

# for (i in 1:nrow(gene1))
# {
#   if (!length(setdiff(gene1$GeneList[i],'A_HC')))
#   {gene1[i,'GeneList']<-paste0('HC')}
#   else if (!length(setdiff(gene1$GeneList[i],'B_AD')))
#   {gene1[i,'GeneList']<-paste0('AD')}
#   else if (!length(setdiff(gene1$GeneList[i],'C_PD')))
#   {gene1[i,'GeneList']<-paste0('PD')}
#   else if (!length(setdiff(gene1$GeneList[i],'D_SVD')))
#   {gene1[i,'GeneList']<-paste0('SVD')}
#   else if (!length(setdiff(gene1$GeneList[i],'F_MS')))
#   {gene1[i,'GeneList']<-paste0('MS')}
#   
# }

#gene1$GeneList1<-NULL
# gene1$GeneList[A_HC]<-'HC'
# gene1$GeneList[B_AD]<-'AD'
# gene1$GeneList[A_PD]<-'PD'
# gene1$GeneList[A_SVD]<-'SVD'
# gene1$GeneList[A_MS]<-'MS'
lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 16, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))


#pdf(file="Popo_figure.pdf")

gene1$GeneList <-factor(gene1$GeneList,levels=c('HC','AD','PD','SVD','MS'))

p<-ggplot(gene1,aes(x=GeneList,y=reorder(Description,Log.q.value.),size=Ratio,colour=Log.q.value.))+
#p<-ggplot(gene1,aes(x=GeneList,y=Description,size=Ratio,colour=Log.q.value.))+
  geom_point(shape=16)+
  labs(X=" ",y=" ")+
  scale_colour_continuous(name="-log(pFDR)", low="blue", high="red")+
  scale_radius(range = c(1,6),name="GeneNumber")+
  
  # theme(axis.text   = element_text(size = rel(0.8)), 
  #       strip.text  = element_text(size = rel(0.8)),
  #       legend.text = element_text(size = rel(0.8)),
  #       plot.title  = element_text(size = rel(1.2)),
  #       panel.grid.minor = element_line(size = rel(0.5)))+
  # guides(color=guide_colorbar(order=1),size=guide_legend(order=2))+
  theme(text = element_text(family = "Arial",
                            size=14,colour='black'),
        axis.text = element_text(colour='black'),
        axis.text.x = element_text(colour='black'),
        axis.text.y = element_text(colour='black'),
        axis.title.x = element_text(colour='black'),
        axis.title.y = element_text(colour='black'),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12))
  
 # theme_bw(base_size = 20, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt)
#print(p)
#dev.off()
ggsave("Rplot01.pdf",p,device = cairo_pdf,width=30,height=12,units='cm')