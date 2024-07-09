library(GO.db)
library(org.Hs.eg.db)
library(fgsea)
library(tidyr)
library(dplyr)


# Pathway Download --------------------------------------------------------


###### KEGG pathways
# the bellow csv file was taken from the CTD web site
# http://ctdbase.org/downloads/;jsessionid=2FE9A5BD51C354A8D6A0FB4BACDB8998
ptw = read.csv('CTD_genes_pathways_12112023.csv',comment.char = '#',strip.white = F,header = F)
kegg<-ptw[grepl('KEGG',ptw$V4),] 
kegg2gene<-sapply(unique(kegg$V4),function(x)as.character(kegg$V2[kegg$V4==x])) 
kegg2gene_symbol_converted = sapply(kegg2gene,function(x)mapIds(org.Hs.eg.db, x,  'SYMBOL','ENTREZID'))
kegg_annot = kegg[,c(4,3)]
colnames(kegg_annot) = c('pathway','description')
kegg_annot=unique(kegg_annot)
kegg_name_to_annotation = function(x){
  return(kegg_annot$description[kegg_annot$pathway == x])
}
names(kegg2gene_symbol_converted) = sapply(names(kegg2gene_symbol_converted),function(x)kegg_name_to_annotation(x))
# # topas pathways --------------------------------------------------------

path_to_basket_annotation_file =  "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/TUPAC_SCORING_4th gen_5th gen 231206.xlsx"
basket_file = xlsx::read.xlsx(path_to_basket_annotation_file,sheetIndex = 1)


all_topas_basekts = unique(basket_file$BASKET)
all_topas_basket_list = list()
for (i in 1:length(all_topas_basekts)){
  all_topas_basket_list[[i]] = unique(basket_file$GENE.NAME[which(basket_file$BASKET %in%  all_topas_basekts[i])])
}
names(all_topas_basket_list) = all_topas_basekts



jaccrd_similarity_of_two_lists<-function(list1,list2){
  #retunrs n*m from a list with gene sets
  jac<-matrix(NA, nrow = length(list1), ncol = length(list2))
  for (i in 1:length(list1)){
    for (j in 1:length(list2)){
      jac[i,j]<-length(intersect(list1[[i]],list2[[j]]))/length(union(list1[[i]],list2[[j]]))
    } #end for j
  }#end for i
  rownames(jac)<-names(list1)
  colnames(jac)<-names(list2)
  return(jac)
}#end function

jac = jaccrd_similarity_of_two_lists(kegg2gene_symbol_converted,all_topas_basket_list)
hist(as.vector(jac),xlab = 'Jaccard Similarity Distribution')

names(kegg2gene_symbol_converted)[which(apply(jac, 1, max,na.rm=T) %in% max(jac))]
#[1] "Cell cycle - G1/S transition" "TGF-beta signaling"          
 names(all_topas_basket_list)[which(apply(jac, 2, max,na.rm=T) %in% max(jac))]
#[1] "CDK4_6 activity" "TGFb" 