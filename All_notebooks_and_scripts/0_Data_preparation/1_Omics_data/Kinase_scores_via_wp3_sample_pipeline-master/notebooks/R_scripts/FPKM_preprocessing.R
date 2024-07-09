# 1- new Master  ----------------------------------------------------------

# converting the Heidelberg FPKM data for the portal use
# load('MatchedProteinRNA_221216.RData')
# the path to all FPKMs
fpkm_dir = "/media/kusterlab/internal_projects/active/TOPAS/WP31/FPKM"

process_fpkm_data = function(new_fpkm_data,proteomics_patients_df){
  colnames(new_fpkm_data) = gsub(colnames(new_fpkm_data),pattern = 'tumor0',replacement = 'T')
  colnames(new_fpkm_data) = gsub(colnames(new_fpkm_data),pattern = 'tumor',replacement = 'T1')
  colnames(new_fpkm_data) = gsub(colnames(new_fpkm_data),pattern = 'metastasis0',replacement = 'M')
  colnames(new_fpkm_data) = gsub(colnames(new_fpkm_data),pattern = 'metastasis',replacement = 'M1')
  colnames(new_fpkm_data) = gsub(colnames(new_fpkm_data),pattern = '_',replacement = '-')
  colnames(new_fpkm_data) = gsub(colnames(new_fpkm_data),pattern = '-undefined-neoplasia\\d{2}',replacement = '')
  colnames(new_fpkm_data) = gsub(colnames(new_fpkm_data),pattern = '-undefined-neoplasia',replacement = '')
  transctip_ids = colnames(new_fpkm_data)[! colnames(new_fpkm_data) %in% c("EnsemblGeneID" ,"Gene")]
  annotation_df = as.data.frame(transctip_ids)
  colnames(annotation_df) = 'transcriptom'
  annotation_df$proteom = NA
  for (i in 1:nrow(annotation_df)){
    if(any(grep(x=proteomics_patients_df$Sample.name,pattern = annotation_df$transcriptom[i]))){
      annotation_df$proteom[i] = proteomics_patients_df$Sample.name[grep(x=proteomics_patients_df$Sample.name,pattern = annotation_df$transcriptom[i])[1]]
    }
  }
  annotation_df = na.omit(annotation_df)
  RNA_ENS_GENE_MAPPER = new_fpkm_data[,colnames(new_fpkm_data) %in% c("EnsemblGeneID" , "Gene" ), ]
  rownames(new_fpkm_data) = new_fpkm_data$EnsemblGeneID
  new_fpkm_data = new_fpkm_data[,!colnames(new_fpkm_data) %in% c("EnsemblGeneID" ,"Gene" )]
  new_fpkm_data = new_fpkm_data[,colnames(new_fpkm_data) %in% annotation_df$transcriptom]
  for (i in 1:ncol(new_fpkm_data))colnames(new_fpkm_data)[i] = annotation_df$proteom[annotation_df$transcriptom == colnames(new_fpkm_data)[i]]
  return(list(processed_fpkm = new_fpkm_data,RNA_ENS_GENE_MAPPER = RNA_ENS_GENE_MAPPER,list_genes = RNA_ENS_GENE_MAPPER$EnsemblGeneID))
}



proteomics_patients_df = xlsx::read.xlsx("/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/Metadata_PAN CANCER_230621.xlsx",
                                         sheetIndex = 1)
# second group of tha patients from MASTER
fpkm_data_1 = read.csv(paste(fpkm_dir,"FPKM_for_missing_PIDs_2023-01-16.csv",
                             sep = '/'),
                       check.names = F)

fpkm_data_2 = read.csv(paste(fpkm_dir,"FPKM_for_missing_PIDs_2023-07-19.csv",
                             sep = '/'),
                       check.names = F)


fpkm_data_3 = read.csv(paste(fpkm_dir,"FPKM_for_missing_PIDs_2023-07-25.csv",
                             sep = '/'),
                       check.names = F)
master_fpkms = lapply(list(fpkm_data_1,fpkm_data_2,fpkm_data_3), function(x)process_fpkm_data(x,proteomics_patients_df))
list_genes = lapply(master_fpkms,function(x)x$list_genes)
common_genes = Reduce('intersect',list_genes)
ordered_fpkms = lapply(master_fpkms,function(x)x$processed_fpkm[common_genes,])
new_fpkm_data = do.call('cbind',ordered_fpkms)
RNA_ENS_GENE_MAPPER = master_fpkms[[1]]$RNA_ENS_GENE_MAPPER
# 2- others ---------------------------------------------------------------

# converting the Heidelberg FPKM data for the portal use
# load('MatchedProteinRNA_221216.RData')
# the path to all FPKMs
library(magrittr)
library(dplyr)
library(tidyr)
patientRNA2protein = function(x){
  # @ :parma x: is the patient id 
  ans = coldata$Proteomics.SampleID[which(coldata$RNA.SampleID %in% x)]
  print(ans[1]) # the first patient will be taken replicates ignored
  return(as.character(ans[1]))
}




# 1st group of the MASTER patients from the Heidelburg data
load(paste(fpkm_dir,"MatchedProteinRNA_221219.RData",sep = '/'))
rna_df = as.data.frame(assay_rna)
coldata = coldata[,c(8,9)]
ids = sapply(colnames(rna_df),function(x)patientRNA2protein(x))
colnames(rna_df) = ids



# INFORM patients from the kids center
INFORM_data = read.csv(paste(fpkm_dir,'INFORM_specIDs_FPKM_values_230323.tsv',sep = '/'),
                       sep = '\t',stringsAsFactors = F,check.names = F)

mapped_INFORM_DF = merge(INFORM_data,RNA_ENS_GENE_MAPPER,by.x='gene_id',by.y='EnsemblGeneID')



int_proteins = intersect(rownames(rna_df),rownames(new_fpkm_data))
new_fpkm_data = new_fpkm_data[int_proteins,]
rna_df = rna_df[int_proteins,]
rna_df = cbind(rna_df[,setdiff(colnames(rna_df),colnames(new_fpkm_data))],
               new_fpkm_data[,setdiff(colnames(new_fpkm_data),colnames(rna_df))])



rna_df$protein_ens = rownames(rna_df)
rowdata = rowdata[,c(3,4)]

# mapping gene ids to protein names
rna_df_protein_mapped = merge(rna_df,rowdata,by.x='protein_ens',by.y = 'RNA.EnsemblGeneID' )
rna_df_protein_mapped = rna_df_protein_mapped[!duplicated(rna_df_protein_mapped$Protein),]
rownames(rna_df_protein_mapped) = rna_df_protein_mapped$Protein
rna_df_protein_mapped = rna_df_protein_mapped[,!colnames(rna_df_protein_mapped) %in% c("protein_ens",'Protein')]

rna_df_protein_mapped$Gene_names = rownames(rna_df_protein_mapped)
rna_df_protein_mapped = rna_df_protein_mapped %>%  mutate(Gene_names=strsplit(as.character(Gene_names),';',fixed =T)) %>%  unnest(Gene_names)
rna_df_protein_mapped = as.data.frame(rna_df_protein_mapped)
rownames(rna_df_protein_mapped) = rna_df_protein_mapped$Gene_names
rna_df_protein_mapped = rna_df_protein_mapped[,!colnames(rna_df_protein_mapped) %in% 'Gene_names']

# adding INFORM patients
intersecting_genes = intersect(rownames(rna_df_protein_mapped),mapped_INFORM_DF$Gene)
mapped_INFORM_DF = mapped_INFORM_DF[mapped_INFORM_DF$Gene %in% intersecting_genes,]
mapped_INFORM_DF = mapped_INFORM_DF[!duplicated(mapped_INFORM_DF$Gene),]
rownames(mapped_INFORM_DF) = mapped_INFORM_DF$Gene
mapped_INFORM_DF = mapped_INFORM_DF %>% select(matches('I0'))

mapped_INFORM_DF = mapped_INFORM_DF[intersecting_genes,]
rna_df_protein_mapped = rna_df_protein_mapped[intersecting_genes,]
if(all(rownames(rna_df_protein_mapped)==rownames(mapped_INFORM_DF))){
  rna_df_protein_mapped = cbind(rna_df_protein_mapped,mapped_INFORM_DF)
}


# # data for the data 12.12.2023 ------------------------------------------
meta_data_path = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/METADATA_PAN CANCER_Batch159_AS.xlsx'
file_path = paste0(fpkm_dir,'/FPKM_for_missing_PIDs_2023-12-12.csv')
df = read.csv(file_path)
meta_Data = xlsx::read.xlsx(meta_data_path,sheetIndex = 1)
new_Df = df[,!colnames(df) %in% c("EnsemblGeneID"  , "Gene" )]
colnames(new_Df) = colnames(new_Df) %>% strsplit('_') %>% sapply('[[',1) %>% gsub(pattern = ('H021.'),replacement = '')

meta_Data$newsample = meta_Data$Sample.name %>% strsplit('-') %>% sapply('[[',2)
id_to_samplename = function(x){
  return(meta_Data$Sample.name[which(meta_Data$newsample==x)][1])
}

colnames(new_Df) = sapply(colnames(new_Df),function(x)id_to_samplename(x))
new_Df$Gene = df$Gene
new_Df = new_Df[!duplicated(new_Df$Gene),]
rownames(new_Df) = new_Df$Gene

intersected_genes = intersect(rownames(rna_df_protein_mapped),new_Df$Gene)
new_Df = new_Df[,!colnames(new_Df) %in% 'Gene']
new_Df = new_Df[intersected_genes,]
rna_df_protein_mapped = rna_df_protein_mapped[intersected_genes,]
duplicated_patients = intersect(colnames(rna_df_protein_mapped),colnames(new_Df))
new_Df = new_Df[,!colnames(new_Df) %in%  duplicated_patients]


rna_df_protein_mapped = cbind(rna_df_protein_mapped,new_Df)



# log transformations
rna_df_protein_mapped = log10(rna_df_protein_mapped)

save(rna_df_protein_mapped,file = paste(fpkm_dir,"log10_fpkm_informa_added.RData",sep = '/'))
# calculating the z_score for each gene across all patients
z_rna_df_protein_mapped =  rna_df_protein_mapped
sds = rep(0,nrow(rna_df_protein_mapped))
for (j in 1:nrow(rna_df_protein_mapped)){
  
  sample = as.numeric(rna_df_protein_mapped[j,]) # a gene across all patients
  sample[is.infinite(sample)] = NA
  mean_sample = median(sample,na.rm = T)
  sd_sample = sd(sample,na.rm = T)
  sds[j] = sd_sample
  print(paste(j,sd_sample,mean_sample))
  z_sample = (sample - mean_sample) /sd_sample
  z_rna_df_protein_mapped[j,] = z_sample
}
z_rna_df_protein_mapped = z_rna_df_protein_mapped[!sds ==0,]
z_rna_df_protein_mapped = z_rna_df_protein_mapped[!is.na(sds),]
colnames(z_rna_df_protein_mapped) = paste(colnames(z_rna_df_protein_mapped), 'Z-score',sep = ' ')
z_rna_df_protein_mapped$genes = rownames(z_rna_df_protein_mapped)

rna_df_protein_mapped = rna_df_protein_mapped[which(rownames(rna_df_protein_mapped) %in% rownames(z_rna_df_protein_mapped)),]
colnames(rna_df_protein_mapped) = paste(colnames(rna_df_protein_mapped), 'Z-score',sep = ' ')
rna_df_protein_mapped$genes = rownames(rna_df_protein_mapped)

# saving the files for the MTB
write.csv(z_rna_df_protein_mapped,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/FPKM/fpkms_230725mitINFORM.csv')
write.csv(rna_df_protein_mapped,file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/FPKM/not_z_scored_fpkms_230725mitINFORM.csv')
