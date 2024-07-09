report_directory = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.03.21_CJ_batch91'
protein_abundance_file = paste(report_directory,'full_proteome_measures_z.tsv',sep = '/')
fpkms_file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/FPKM/fpkms_221219mitINFORM.csv'
fpkm_log10_abundance = "/media/kusterlab/internal_projects/active/TOPAS/WP31/FPKM/log10_fpkm_informa_added.RData"

# correlation_analysis across all proteins----------------------------------------------------

protein_abundance = read.csv(protein_abundance_file,
                             sep = '\t',
                             check.names = F)

rownames(protein_abundance) = protein_abundance$`Gene names`

fpkms = read.csv(fpkms_file,
                 check.names = F)

rownames(fpkms) = fpkms[, 1]

colnames(fpkms) = gsub(colnames(fpkms), pattern = ' Z-score', replacement = '')
colnames(protein_abundance) = gsub(colnames(protein_abundance), pattern = 'zscore_', replacement = '')

intersected_patients = intersect(colnames(fpkms),colnames(protein_abundance))
intersected_proteins = intersect(rownames(fpkms),rownames(protein_abundance))

fpkms = fpkms[intersected_proteins,intersected_patients]
protein_abundance = protein_abundance[intersected_proteins,intersected_patients]


final_df = as.data.frame(list(rownames(protein_abundance)))
colnames(final_df) = 'protein'
final_df$pearson_correlation = NA
final_df$spearman_correlation = NA
final_df$number_patients = NA
final_df$pearson_p_value = NA
final_df$spearman_p_value = NA

for (i in 1:nrow(final_df)){
  print(i)
  protein_name = final_df$protein[i]
  df = as.data.frame(list(t(protein_abundance[protein_name, ]), t(fpkms[protein_name,])))
  df = na.omit(df)
  if (nrow(df) > 2){
    cor_test = cor.test(df[, 1], df[, 2], method ='pearson')
    final_df$pearson_p_value[i] = cor_test$p.value
    cor_test = cor.test(df[, 1], df[, 2], method ='spearman')
    final_df$spearman_p_value[i] = cor_test$p.value
  }
  final_df$pearson_correlation[i] = cor(df, method = 'pearson')[1,2]
  final_df$spearman_correlation[i] = cor(df, method = 'spearman')[1,2]
  final_df$number_patients[i] = nrow(df)
}
final_df = na.omit(final_df)


# Correlation analysis ----------------------------------------------------

# Density plot for the significant none significant correlations ----------------------------------------------------------------
library(plotly)
method = 'spearman'

correlation_type = paste(method,'correlation',sep = '_')
correlation_p_value = paste(method,'p_value',sep = '_')

data = final_df
data = data[data$number_patients > 30 ,]


#data = data[!is.na(data$basket_annotation),]  # proteins of interest

data$significant = 'not_significant'
data$significant[which(data[[correlation_p_value]] < 0.05)] = 'significant'

section1 <- data[which(data$significant == "significant" & 
                         data[correlation_type] > 0 ),]

section2 <- data[which(data$significant != "significant" & 
                         data[correlation_type] > 0 ),]

section3 <- data[which(data[correlation_type] <= 0),]

dx1 = table(cut(section1[[correlation_type]],breaks = seq.int(from = 0, to = 1, by = 0.01)))
dx1_range = seq(.005,.995,.01)
df1 = as.data.frame(list(dx1,dx1_range))

dx2 = table(cut(section2[[correlation_type]], breaks = seq.int(from = 0, to = 1, by = 0.01)))
df2 = as.data.frame(list(dx2,dx1_range))
# negative correlation
dx3 = table(cut(section3[[correlation_type]], breaks = seq.int(from = -1, to = 0, by = 0.01)))
dx3_range = seq(-0.995,0,.01)
df3 = as.data.frame(list(dx3,dx3_range))

colnames(df1) = colnames(df2) = colnames(df3) = c('range','y','x')


clrs <- c('#adcbe3','#4b86b4','#2a4d69','#e78ac3')
fig <- plot_ly(x = ~df1$x,
               y = ~df1$y,
               type = 'scatter',
               mode = 'lines',
               name = 'Significant (p_value < 0.05)',
               line = list(color = clrs[1], width = 1),
               fillcolor = clrs[1],
               fill = 'tozeroy')

fig <- fig %>% add_trace(x = ~df2$x,
                         y = ~df2$y,
                         #text = ~df2$range,
                         name = 'not significant',
                         line = list(color = clrs[2], width = 1),
                         fill = 'tozeroy',
                         fillcolor = clrs[2])

fig <- fig %>% add_trace(x = ~df3$x,
                         y = ~df3$y, 
                         name = 'negative correlation',
                         fill = 'tozeroy',
                         line = list(color = clrs[3], width = 1),
                         fillcolor = clrs[3])

fig <- fig %>% layout(xaxis = list(zeroline = F,
                                   # zerolinecolor = 'white',
                                   #zerolinewidth = 0,
                                   title = correlation_type),
                      yaxis = list(
                        title = '# Proteins'))
fig


