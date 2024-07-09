# for all proteins it makes  the swarm plots as one multi page pdf files
library(ggplot2)
library(ggbeeswarm)
library(gridExtra)
library(patchwork)

df = read.csv('/media/kusterlab/internal_projects/active/CLINSPECT_WP4/MTB_pipeline_FP/23.04.21_CJ_WP3_pipeline/preprocessed_fp.csv')
file_to_save = 'glioma_FP_across_all_patients.pdf'


get_sub_df = function(df,protein_name){
  # @ df: a data frame with rownames as the protein names and the columns are the patients
  # @ protein_name: the protein_name to extract
  sub_df = as.data.frame(t(df[protein_name,]))
  colnames(sub_df) = 'intensity'
  sub_df$x = protein_name
  return(sub_df)
  
}

# swarm plot function in R
draw_swarm <- function(sub_df,protein_name,min_df,max_df){
  # @ sub_df: data_frame with two columns x: protein_names 
  p <- ggplot(sub_df, aes(
    x = x,
    y = intensity,
    color = intensity
  )) +
    geom_quasirandom() +
    scale_y_continuous( limits = c(min_df, max_df)) +
    theme(axis.text.x = element_text(angle = 45)) +
    xlab('Distribution of all patients') + ylab('log10 (MAXLFQ_intensities)') + theme_bw()
  return(p)
}



df = df[,-grep(colnames(df),pattern = 'Identification')]

rownames(df) = df$Gene.names
df = df[,!colnames(df) %in% 'Gene.names']
df$stdev = apply(df, 1,sd,na.rm=T)
df$counts = apply(df,1,function(x)sum(!is.na(x)))
df = df[order(df$counts,df$stdev,decreasing = T),]



df = df[,!colnames(df) %in% c('stdev','counts')]

min_df =  min(df,na.rm = T)
max_df =  max(df,na.rm = T)



pdf(file_to_save, onefile = TRUE)
for (i in seq(1,nrow(df),2)) {
 print(i)
 j=i+1 
  protein_name = rownames(df)[i]
  sub_df = get_sub_df(df,protein_name)
  p1 = draw_swarm(sub_df,protein_name,min_df,max_df)
  
  protein_name = rownames(df)[j]
  sub_df = get_sub_df(df,protein_name)
  p2 = draw_swarm(sub_df,protein_name,min_df,max_df)
  #print(grid.arrange(p1, p2, nrow = 1,ncol=2,))
  combined <- p1 + p2 & theme(legend.position = "bottom")
  print(combined)
}
dev.off()
