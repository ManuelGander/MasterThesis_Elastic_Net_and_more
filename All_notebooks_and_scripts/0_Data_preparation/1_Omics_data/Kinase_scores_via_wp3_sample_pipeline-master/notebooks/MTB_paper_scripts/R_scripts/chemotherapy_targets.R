# source: https://db.idrblab.net/ttd/full-data-download           P1-06-Target_disease.txt       P2-01-TTD_uniprot_all.txt
# goal: find the number of the drugtargetable proteins in the MASTER/INFORM identified proteins
#       find the muber of identified kinases in the MASTER/INFORM patients
#       
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
# wp3 output dir
output_dir = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2024.01.13_CJ_pancancer_166_new_PSP'

# db: drug targetable proteins source: https://db.idrblab.net/ttd/full-data-download
df = read.csv('/home/amir/topas/Annika/cancer_drug_targets/P1-06-Target_disease.txt',
              sep = '\t',row.names = NULL,header = F)

colnames(df) = c('target','indication','targetID','disease')

# normalized intensities of the wp3 pipeline output
intensitities = paste0(output_dir,'/preprocessed_fp.csv')
int_df = read.csv(intensitities,check.names = F)

all_psites_path = paste0(output_dir,'/preprocessed_pp.csv')
all_psites_df = read.csv(all_psites_path) 

uniprot_mapping = read.csv('/home/amir/topas/Annika/cancer_drug_targets/P2-01-TTD_uniprot_all.txt',sep = '\t',row.names = NULL,header = F)
colnames(uniprot_mapping) = c('target','indication','uniprot_ID')
# adding kinases
kinases_df = xlsx::read.xlsx('Kincat_Hsap.08.02_FB.xls',sheetIndex = 1)
mapped_gene_to_kinase = read.csv('Kinase_to_Genes.txt',comment.char = '#',sep = '\t')
kinhub_list = read.csv('kinhub.csv')
##############
uniprot_mapping = uniprot_mapping[uniprot_mapping$indication == 'UNIPROID',]
df = df[df$target != '',]
cancer_related_df = df[grep(x= df$disease,
                            pattern = 'malignan|cancer|neoplasia|carcinoma|Melanoma|lymphoma|neoplasm|sarcoma|tumor|leukaemia|tumour',
                            ignore.case = T),] 


final_df = merge(cancer_related_df,uniprot_mapping,by='target')
final_df$uniprot_ID = gsub(final_df$uniprot_ID,pattern = '_HUMAN',replacement = '')


# to compare the db list with our list

int_df = tidyr::gather(int_df,key = 'patients',value = 'intensities',-`Gene names`)

# un-netsting the protein groups by ;
int_df = int_df %>%   
  mutate(`Gene names`=strsplit(as.character(`Gene names`),';',fixed =T)) %>%
  unnest(`Gene names`)

# counting the number of the patients per protein
int_df = int_df[-grep(int_df$patients,pattern = 'metadata'),]
int_df = na.omit(int_df)
int_df <- int_df %>% 
  group_by(`Gene names`) %>% 
  mutate(count = n())


### DRUG TARGETABLE PROTOEINS
int_df$targetable = 'not targetable'
all_cancer_targsts = final_df$uniprot_ID %>% 
  strsplit(';') %>%
  unlist() %>% strsplit('-') %>% 
  unlist() %>% gsub(pattern = ' ',replacement = '') %>%
  unique()

int_df = int_df[,colnames(int_df) %in% c("Gene names",  "count", "targetable" )]
int_df = unique(int_df)

int_df$targetable[int_df$`Gene names` %in% all_cancer_targsts] = 'Targetable'

### KINASE 
### KINASE 
# we use the union of all databases reporting the kinases
kinses = xlsx::read.xlsx('Full Kinase List_human kinome.xlsx',sheetIndex = 1)
all_kinases_names = as.character(kinses$HGNC.Name[!is.na(kinses$HGNC.Name)])

#all_kinases_names = Reduce('union',list(kinases_df$Name,
#                                        mapped_gene_to_kinase$Kinase, # from flo list
#                                        mapped_gene_to_kinase$Gene,  
#                                        kinases_df$Name,
#                                        kinhub_list$Manning.Name,     # from the kin hub
#                                        kinhub_list$HGNC.Name,
#                                        kinhub_list$Kinase.Name,
#                                        kinhub_list$xName))

# we use the union of all databases reporting the kinases

int_df$kinase = 'not_kinase'
int_df$kinase[int_df$`Gene names` %in% all_kinases_names] = 'kinase'
print(paste0(length(unique(int_df$`Gene names`[int_df$kinase == 'kinase'])),' out of ',length(unique(kinases_df$Name)), ' kinases are found'))
print(paste0(length(unique(int_df$`Gene names`[int_df$targetable == 'Targetable'])),' out of ',length(unique(all_cancer_targsts)), ' drug targetable proteins are found'))



all_psites_df = all_psites_df[,-grep(colnames(all_psites_df),pattern = 'meta')]
all_psites_df = all_psites_df %>%
  mutate(count = rowSums(!is.na(select(., -c("Gene.names" ,"Modified.sequence",  "Proteins")))))
all_psites_df = all_psites_df[which(all_psites_df$count != 0),]

kinase_scored_peptides_path = paste0(output_dir,'/kinase_results/scored_peptides.tsv')
kinase_scores_df = read.csv(kinase_scored_peptides_path,sep = '\t')



peptides_df = all_psites_df[,colnames(all_psites_df) %in% c('Gene.names','Modified.sequence')]
peptides_df <- peptides_df %>%  group_by(Gene.names) %>%  mutate(count = n())
peptides_df = unique(peptides_df[,colnames(peptides_df) %in% c("Gene.names","count" )])

kinase_substrate_df = kinase_scores_df[,colnames(kinase_scores_df) %in% c('PSP.Kinases','Modified.sequence')]
kinase_substrate_df <- kinase_substrate_df %>%  group_by(PSP.Kinases) %>%  mutate(count = n())
kinase_substrate_df = unique(kinase_substrate_df[,colnames(kinase_substrate_df) %in% c("PSP.Kinases","count" )])


colnames(kinase_substrate_df) = c("Gene names",'number_identified_substrates')
colnames(peptides_df) = c("Gene names",'number_identified_peptides')

int_df = merge(int_df,kinase_substrate_df,all.x = T)
int_df = merge(int_df,peptides_df,all.x = T)

RTKs_list = c('ABL1', 'ABL2', 'ALK', 'AXL', 'BLK', 'CSK', 'DDR1', 'DDR2', 'EGFR', 'EPHA1', 'EPHA10', 'EPHA2', 'EPHA3', 'EPHA4', 'EPHA5', 'EPHA7', 'EPHA8', 'EPHB1', 'EPHB2', 'EPHB3', 'EPHB4', 'ERBB2', 'ERBB3', 'ERBB4', 'FAK', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FGR', 'FLT1', 'FLT4', 'FYN', 'HCK', 'IGF1R', 'JAK1', 'JAK2', 'JAK3', 'KDR', 'KIT', 'LCK', 'LYN', 'MERTK', 'MET', 'NTRK1', 'NTRK2', 'NTRK3', 'PDGFRA', 'PDGFRB', 'RET', 'SRC', 'YES1')

int_df$is_RTK = FALSE
int_df$is_RTK[int_df$`Gene names` %in% RTKs_list] = T


subset = int_df[int_df$kinase=='kinase',colnames(int_df) %in% c('Gene names','count')]
write.csv(subset,file='occurece_kinases.csv')
colnames(subset) = c('Protein','number_patients')




## plotting
ggplot(subset, aes(x = reorder(Protein,-number_patients), y = number_patients)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#geom_text(aes(label = number_patients), vjust = 1.5, colour = "white")


#barplot(subset$number_patients, 
#        main = "Occurence",
#        xlab = "Protein", 
#        ylab = "Num_patients", 
#        names.arg = subset$Protein)


#write.csv(int_df,file = 'data_of_identification.csv')

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

found_targets_length = length(unique(int_df$`Gene names`[int_df$targetable == 'Targetable']))
all_targetes_length = length(unique(all_cancer_targsts))
Prop = c((all_targetes_length-found_targets_length),found_targets_length)
pie(Prop,labels = c('Undetected Targetable proteins','Detected Targetable proteins'))


found_kinasea_length = length(unique(int_df$`Gene names`[int_df$kinase == 'kinase']))
all_kinases_length = length(unique(kinases_df$Name))
Prop = c((all_kinases_length-found_kinasea_length),found_kinasea_length)
pie(Prop,labels = c('Undetected Kinases','Detected Kinases'))


all_id_peptided = unique(all_psites_df$Modified.sequence)
all_id_peptides_kinase = unique(all_psites_df$Modified.sequence[all_psites_df$Gene.names %in% unique(kinases_df$Name) ])
all_id_peptides_kinase_substrated = unique(kinase_scores_df$Modified.sequence[kinase_scores_df$PSP.Kinases %in% unique(kinases_df$Name) ])

Prop = c((length(all_id_peptided)-length(all_id_peptides_kinase)),length(all_id_peptides_kinase))
pie(Prop,labels = c('Other','Peptides of Kinases'))


Prop = c((length(all_id_peptided)-length(all_id_peptides_kinase_substrated)),length(all_id_peptides_kinase_substrated))
pie(Prop,labels = c('Other','Peptides of Kinases Substrates'))


