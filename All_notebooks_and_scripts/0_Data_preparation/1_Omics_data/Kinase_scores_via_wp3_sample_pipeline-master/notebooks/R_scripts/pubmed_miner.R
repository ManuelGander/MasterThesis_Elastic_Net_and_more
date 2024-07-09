# this script can search for a specific biomaker asoociated with a disease in litreature
# first search for the term like Chordoma in Pubmed and then select save pubmed and use
# as the input for this script


library(pubmed.mineR)
library(magrittr)
library('org.Hs.eg.db') 
library(xlsx)

abstracts <- readabs('pubmed-Chordoma-set.txt')
pmids<-abstracts@PMID
res<-list()
for (i in 1:length(pmids)){res[[i]]<-pubtator_function(pmids[i]);print(i)}
res[which(sapply(res, length)==1)]<-c()
res_genes<-lapply(res, function(x)x$Genes) # all genes reported in the papers
res_genes<-unlist(res_genes)


entrez<-res_genes %>% strsplit('>') %>% lapply('[[', 2) %>% unlist() %>% strsplit(';') %>% unlist()
gene_freq<-as.data.frame(table(entrez))


# mappping the gene Entrez to the genes symbols

gene_freq$symbol<-mapIds(org.Hs.eg.db, as.character(gene_freq$entrez),  'SYMBOL','ENTREZID')

# saving as Rdata and excel
save(gene_freq,res,file = 'chordoma_genes_frequency_14NOV2023_result.RData')
xlsx::write.xlsx(gene_freq,'publications_frequency_chordoma_14NOV2023.xlsx',sheetName = 'Frequency')
