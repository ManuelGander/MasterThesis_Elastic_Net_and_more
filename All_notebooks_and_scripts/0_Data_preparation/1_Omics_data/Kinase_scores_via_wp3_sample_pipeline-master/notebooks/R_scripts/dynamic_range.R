report_directory = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.12.29_AhS_batch73_74_imp'
protein_abundance_file = paste(report_directory,'full_proteome_measures_z.tsv',sep = '/')
fpkms_file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/FPKM/fpkms_221219.csv'
normalized_intensities_file = paste(report_directory,'preprocessed_fp.csv',sep = '/')
sample_annotation_file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Searches/patient_annotation_mixed_cohort_221229.csv'
basket_annotation_file = "/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/Baskets_4th gen_230103.xlsx"
fpkm_log10_abundance = "/media/kusterlab/internal_projects/active/TOPAS/WP31/FPKM/log10_fpkm.RData"
ibaq_table = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/MT/2022.01.10_iBAQ_tables/ibaq_table.tsv'
# correlation_analysis ----------------------------------------------------

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


# protein_per_batch_counts ----------------------------------------------

get_proteins_per_batch = function(batch_number,normalized_intensities,sample_annotation,num_peptides_data){
  patients_name = sample_annotation$Sample.name[sample_annotation$Batch.Name == batch_number]
  intensity_per_patients = normalized_intensities[,which(colnames(normalized_intensities) %in% patients_name)]
  intensity_per_patients = na.omit(intensity_per_patients)
  primarylist_proteins = rownames(intensity_per_patients)
  num_pep_batch = as.data.frame(num_peptides_data[primarylist_proteins,which(colnames(num_peptides_data) %in% patients_name)])
  num_pep_batch$min_peps = as.vector(apply(num_pep_batch,1,min,na.rm =T))
  return(rownames(num_pep_batch)[which(num_pep_batch$min_peps > 1)])
}

normalized_intensities = read.csv(normalized_intensities_file,check.names = F)
rownames(normalized_intensities) = normalized_intensities$`Gene names`

num_peptides_data = normalized_intensities[,grep(colnames(normalized_intensities),pattern = 'metadata')]
rownames(num_peptides_data) = normalized_intensities$`Gene names`

num_peptides_data = apply(num_peptides_data,2,function(x)gsub(x, pattern = 'num_peptides=|;',replacement = ''))
num_peptides_data = apply(num_peptides_data,2,function(x)gsub(x, pattern = 'detected in batch',replacement = ''))
colnames(num_peptides_data) = gsub(colnames(num_peptides_data), pattern = "Identification metadata ",replacement = '')


sample_annotation = read.csv(sample_annotation_file)

Gene_names_uniq = unique(normalized_intensities$`Gene names`)
batch_uniq = unique(sample_annotation$Batch.Name)

batch_protein_mat = matrix(0,nrow = length(batch_uniq) ,ncol = length(Gene_names_uniq))
rownames(batch_protein_mat) = batch_uniq
colnames(batch_protein_mat) = Gene_names_uniq

for (i in 1:nrow(batch_protein_mat)){
  print(i)
  batch_no = rownames(batch_protein_mat)[i]
  list_protein = get_proteins_per_batch(batch_no,normalized_intensities,sample_annotation,num_peptides_data)
  batch_protein_mat[i,which(colnames(batch_protein_mat) %in% list_protein)] = 1
}
batch_protein_df = as.data.frame(t(batch_protein_mat))

batch_protein_df$percent_overlapping = (rowSums(batch_protein_df)/ncol(batch_protein_df))*100
batch_protein_final_df = as.data.frame(batch_protein_df[,'percent_overlapping'])
colnames(batch_protein_final_df) = 'overlapping_percent'
batch_protein_final_df$protein = rownames(batch_protein_df)

final_df_overlapping = merge(final_df,batch_protein_final_df)


# adding basket_annotations ------------------------------------------------------

basket_annotation =xlsx::read.xlsx(basket_annotation_file,sheetIndex = 1)
basket_annotation = basket_annotation[,c("BASKET", "GENE.NAME" )]
basket_annotation = unique(basket_annotation)
final_df_overlapping$basket_annotation = NA
for (j in 1: nrow(basket_annotation)){
  final_df_overlapping$basket_annotation[which(final_df_overlapping$protein %in% basket_annotation$GENE.NAME[j])] = basket_annotation$BASKET[j]
  
}

xlsx::write.xlsx(final_df_overlapping,paste(report_directory,'Correlation_with_fpkm_values.xls',sep = '/'))
# Density plot for the significant none significant correlations ----------------------------------------------------------------
library(plotly)
method = 'spearman'

correlation_type = paste(method,'correlation',sep = '_')
correlation_p_value = paste(method,'p_value',sep = '_')


data = final_df_overlapping[final_df_overlapping$overlapping_percent > 90,]
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




# ibaq vs fpkm density plot--------------------------------------------------------------------

protein_ibaqs = read.csv(ibaq_table,sep = '\t') # proteins and their corresponding in silico No.peptides

proteins_of_interest = final_df_overlapping$protein[!is.na(final_df_overlapping$basket_annotation)]

#intersected_proteins = proteins_of_interest #to just do for proteins of interest (comment this line to do it for all proteins)

intersected_proteins = intersected_proteins[!sapply(intersected_proteins, function(x)sapply(strsplit(x,';'),length)) > 1]


intensities_df = normalized_intensities[rownames(normalized_intensities) %in% intersected_proteins,
                                        colnames(normalized_intensities) %in% intersected_patients]



IBAQ_values = intensities_df
IBAQ_values[,] = NA


# calculation of the IBAQ values
for (protein_idx in 1:length(rownames(intensities_df))){
  protein_name = rownames(intensities_df)[protein_idx]
  peptide_num = protein_ibaqs$num.iBAQ.peptides[protein_ibaqs$Gene.name == protein_name]
  print(paste(protein_name,peptide_num,sep = '_'))
  IBAQ_values[protein_idx,] = log10(10^(intensities_df[protein_idx,]) / peptide_num)
  
}
protein_ibaq = IBAQ_values

load(fpkm_log10_abundance)
fpkm_log10 = rna_df_protein_mapped[rownames(rna_df_protein_mapped) %in% intersected_proteins,
                                   colnames(rna_df_protein_mapped) %in% intersected_patients,]
density_protein <- density(as.vector(protein_ibaq[!is.na(protein_ibaq)]))
density_fpkm <- density(as.vector(fpkm_log10[!is.na(fpkm_log10)]))

fig <- plot_ly(x = ~density_fpkm$x,
               y = ~density_fpkm$y,
               type = 'scatter',
               mode = 'lines',
               name = 'mRNA[FPKM]',
               line = list(color = clrs[1], width = 1),
               fillcolor = clrs[1],
               fill = 'tozeroy')

fig <- fig %>% add_trace(x = ~density_protein$x,
                         y = ~density_protein$y,
                         name = 'Protein [iBAQ]',
                         line = list(color = clrs[2], width = 1),
                         fill = 'tozeroy',
                         fillcolor = clrs[2])


fig <- fig %>% layout(xaxis = list(zeroline = F,
                                   # zerolinecolor = 'white',
                                   # zerolinewidth = 0,
                                   title = 'Log10[abundance]'),
                      yaxis = list(
                        title = 'Density'))
fig



