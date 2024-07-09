library(tidyr)
report_dir = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.08.18_CJ_pancancer_celllines_pdx"
#basket_annotation_path = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/TUPAC_SCORING_4th gen_230628.xlsx"
BASKET_SCORES_FILE = 'basket_scores_4th_gen_zscored.tsv'
Tupac_classification_file_path = '/home/amir/Desktop/Annika_files/Tupac/TUPAC classification.xlsx'
####


basket_df = read.csv(paste(report_dir,BASKET_SCORES_FILE,sep = '/'),sep = '\t',check.names = F)
tupac_classification_df = xlsx::read.xlsx(Tupac_classification_file_path,sheetIndex = 1)
tupac_classificaiton_types = unique(tupac_classification_df$classification)
#basket_annotation_df = xlsx::read.xlsx(basket_annotation_path,sheetIndex = 1)



patient = c('H021-QSF764-T1-Q1')
df = basket_df[basket_df$Sample == patient,intersect(colnames(basket_df),
                                                     tupac_classification_df$TUPAC)]

df = t(df)
colnames(df) = 'y'
df = as.data.frame(df)

df$x = rownames(df)
df  = merge(df,tupac_classification_df,by.x = 'x',by.y = 'TUPAC')

df[is.na(df)] = 0
df[df<0] = 0

df$rank=1:nrow(df)
df = df[order(df$classification),]
categories = unique(df$classification)
labels = df$x
# lolipop Plot
ggplot(df, aes(x=factor(x,levels = labels), y=y)) +
 geom_segment( aes(x=factor(x,levels = labels), xend=x, y=0, yend=y,color=factor(classification,levels = categories))) +
 geom_point( aes(color=factor(classification,levels = categories), size=y)) +
  theme_light() +
  theme(
   # legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab(patient) +
  ylab("RTK Baskets")+
  coord_polar()
