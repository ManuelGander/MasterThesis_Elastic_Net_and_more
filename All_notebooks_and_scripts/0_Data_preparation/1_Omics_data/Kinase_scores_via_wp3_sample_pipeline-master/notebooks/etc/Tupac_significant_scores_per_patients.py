import pandas as pd


tupac_scores = pd.read_csv('/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.12.07_CJ_PAPER_final/basket_scores_4th_gen_zscored.tsv',sep='\t')
expresssion_z_scores = pd.read_csv('/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.12.07_CJ_PAPER_final/full_proteome_measures_z.tsv',sep='\t')
meta_data_df = pd.read_excel("/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_MTBs_Evaluation/METADATA_PAPER_Dez23_Batch158.xlsx")


RTK_tupacs = [
'ABL1',
'ABL2',
'ALK',
'AXL',
'BLK',
'CSK',
'DDR1',
'DDR2',
'EGFR',
'EPHA1',
'EPHA10',
'EPHA2',
'EPHA3',
'EPHA4',
'EPHA5',
'EPHA7',
'EPHA8',
'EPHB1',
'EPHB2',
'EPHB3',
'EPHB4',
'ERBB2',
'ERBB3',
'ERBB4',
'FAK',
'FGFR1',
'FGFR2',
'FGFR3',
'FGFR4',
'FGR',
'FLT1',
'FLT4',
'FYN',
'HCK',
'IGF1R',
'JAK1',
'JAK2',
'JAK3',
'KDR',
'KIT',
'LCK',
'LYN',
'MERTK',
'MET',
'NTRK1',
'NTRK2',
'NTRK3',
'PDGFRA',
'PDGFRB',
'RET',
'SRC',
'YES1'
]

down_stream_signaing = ['AKT',
'AURK',
'CDK4_6 activity',
'ERK',
'HIPPO',
'JNK',
'MEK',
'MTOR',
'MTORC1',
'MTORC2',
'NFKB',
'NOTCH',
'P38',
'PI3K',
'PLK1',
'Proximal RTK signaling',
'RAF',
'RAS',
'SHH',
'TGFb',
'WNT',
'cellcycle'
]



tumor_antigen = [
'PTK7',
'MAGEA10',
'TACSTD2',
'CEACAM7',
'LRRC15',
'CEACAM1',
'PRAME',
'PROM1',
'GPC3',
'CDH17',
'FAP',
'NECTIN2',
'PMEL',
'F3',
'FOLH1;FOLH1B',
'MAGEA4',
'MSLN',
'SLC39A6',
'FOLR1',
'CTAG1A',
'NECTIN4',
'SLC18A1',
'PSCA',
'CLDN6',
'SLC18A2',
'SSTR2',
'CLDN18-2',
'CTAG2',
'MAGEA3;MAGEA6',
'CD274'
]



# 
tupac_scores.set_index('Sample',inplace=True)
significant_RTK_df = pd.DataFrame(tupac_scores[tupac_scores[RTK_tupacs] > 2].count(axis=1),columns=['significant (R)TK TUPAC scores (z-score >2)'])
significant_down_stream_signaling_df = pd.DataFrame(tupac_scores[tupac_scores[down_stream_signaing] > 2].count(axis=1),columns =['significant downstram signaling TUPAC scores (z-score >2)'])
all_significant_tupacs = pd.DataFrame(tupac_scores[tupac_scores > 2].count(axis=1),columns=['significant TUPAC scores (z-score >2'])

significant_tupacs = pd.concat([significant_RTK_df,significant_down_stream_signaling_df,all_significant_tupacs],axis=1)
significant_tupacs['Sample_name'] = significant_tupacs.index
expression = expresssion_z_scores.set_index('Gene names').T
significant_tumor_antigens = pd.DataFrame(expression[expression[tumor_antigen] > 2].count(axis=1),columns= ['significant Tumorantigen z-scores (z-score >2)'])
significant_tumor_antigens['Sample_name'] = [x.split('_')[1] for x in significant_tumor_antigens.index.tolist() ]
tupac_df = significant_tumor_antigens.merge(significant_tupacs,on='Sample_name')
res = []
for i in tumor_antigen:
    print(i)
    temp = pd.DataFrame(expression[expression[[i]] > 2].count(axis=1),columns= [i])
    res.append(temp)
final_tumor_antigen_df = pd.concat(res,axis=1)
final_tumor_antigen_df['Sample_name'] = [x.split('_')[1] for x in final_tumor_antigen_df.index.tolist() ]
tupac_df = tupac_df.merge(final_tumor_antigen_df,on='Sample_name')
[x for x in tupac_df['Sample_name'].tolist() if x not in meta_data_df['Sample name'].str.strip().tolist() ]
meta_data_df['Sample name'] = meta_data_df['Sample name'].str.strip()
final_df = meta_data_df[['Sample name','Tumor cell content','code_oncotree']].merge(tupac_df,right_on='Sample_name',left_on='Sample name')
final_df.to_excel('/home/amir/Desktop/significant_tupac_counts.xlsx')
