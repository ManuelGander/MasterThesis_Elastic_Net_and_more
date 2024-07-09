library(ggvenn)
library(magrittr)

filePath = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.11.28_AhS_all_mixed_cohort/phospho_measures_z.tsv'

z_scors = read.delim(filePath)
TALL_patients = c('T022-30615-T1',  #BCR-ABL1 
                  'T022-31857-T1',  #BCR-ABL1
                  'T022-32063-T1',  #BCR-ABL1
                 #'T022-30426-T1',  #ABL1 fusions
                  'T022-29480-T1',  #IgH-CRLF2
                  'T022-29511-T1',  #IgH-CRLF2
                  'T022-29992-T1',  #ETV-RUNX1
                  'T022-30063-T1')  #ETV-RUNX1

TALL_patients = paste('zscore_',gsub(TALL_patients,pattern = '-',replacement = '.'),sep = '')

group_BCR_ABL1 = TALL_patients[1:3]
group_IgH_CRLF2 = TALL_patients[4:5]
group_ETV_RUNX1 = TALL_patients[6:7]

TALL_z_scores = z_scors[,TALL_patients]
rownames(TALL_z_scores) = paste(z_scors$Gene.names,z_scors$Modified.sequence,sep = '_')


# patientswise: list of the genes in one patient with respect to the others
colnames(TALL_z_scores)[colnames(TALL_z_scores) %in% group_BCR_ABL1]='group_BCR_ABL1'
colnames(TALL_z_scores)[colnames(TALL_z_scores) %in% group_IgH_CRLF2]='group_IgH_CRLF2'
colnames(TALL_z_scores)[colnames(TALL_z_scores) %in% group_ETV_RUNX1]='group_ETV_RUNX1'

# ANOVA test across all genes

BCR_ABL1_proteins = TALL_z_scores[,colnames(TALL_z_scores) == 'group_BCR_ABL1'] %>% na.omit() %>% rownames()
IgH_CRLF2_proteins = TALL_z_scores[,colnames(TALL_z_scores) == 'group_IgH_CRLF2'] %>% na.omit() %>% rownames()
ETV_RUNX1_proteins = TALL_z_scores[,colnames(TALL_z_scores) == 'group_ETV_RUNX1'] %>% na.omit() %>% rownames()

unique_protein_in_BCR_ABL1 = setdiff(BCR_ABL1_proteins,unlist(list(IgH_CRLF2_proteins,ETV_RUNX1_proteins)))
unique_protein_in_IgH_CRLF2_proteins = setdiff(IgH_CRLF2_proteins,unlist(list(BCR_ABL1_proteins,ETV_RUNX1_proteins)))

set.seed(20190708)

x <- list(
  BCR_ABL1_Phosphopeptides = BCR_ABL1_proteins, 
  IgH_CRLF2_Phosphopeptides = IgH_CRLF2_proteins, 
  ETV_RUNX1_Phosphopeptides= ETV_RUNX1_proteins
  
)



ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)



TALL_z_scores = na.omit(TALL_z_scores) # no imputations and removing the genes with NA values 
p_values=rep(NA,nrow(TALL_z_scores))
for (i in 1:nrow(TALL_z_scores)){
  
  z = TALL_z_scores[i,]  #z_scores
  rr = as.data.frame(t(z))
  rr$names = names(z)    # group names
  colnames(rr) = c('values','ind')
  res.aov<-aov(values ~ ind,data= rr)
  print(rownames(TALL_z_scores)[i])
  p_values[i] = summary(res.aov)[[1]][["Pr(>F)"]][1]
  
}

names(p_values) = rownames(TALL_z_scores)
p_values = as.data.frame(p_values)
p_values$ProteinNames = rownames(p_values)
# p_adjusted values
p_values$p.adjusted_fdr = p.adjust(p_values$p_values, method ='fdr', n = length(p_values$p_values))

int = p_values[p_values$p_values < .01 & p_values$p.adjusted_fdr <0.3 ,]
sig_proteins = TALL_z_scores[intersect(rownames(int),rownames(TALL_z_scores)),]


sig_proteins$mean_BCR_ABL1 = apply(sig_proteins,1,function(x)mean(x[1],x[2],x[3]))
sig_proteins$mean_IgH_CRLF2 = apply(sig_proteins,1,function(x)mean(x[4],x[5]))
sig_proteins$mean_RUNK1 = apply(sig_proteins,1,function(x)mean(x[6],x[7]))

box_data = sig_proteins[,grep(colnames(sig_proteins),pattern = 'mean')]
box_data$Protein = rownames(box_data)
tt = tidyr::gather(box_data,key='group',value='z_score',-Protein)

library(plotly)
fig <- plot_ly(
  data = tt,
  y = ~z_score,
  x = ~Protein,
  type = "box"
)

fig



write.csv(p_values,'/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.11.28_AhS_all_mixed_cohort/ANOVA_DEGs_pp.csv')
