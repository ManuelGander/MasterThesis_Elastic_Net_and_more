import pandas as pd
import seaborn as sns
import math
import seaborn as sns
import plotly.express as px

regulation_mapping = {
    'up'  : 'Lysate',
    'down': 'KryoTissue'
}

def pfunction(x):
    try:

        return math.log10(float(1/x)) 
    except:
        return 0
path_to_file = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2023.06.22_AhS_PAPER_COHORT/FP_signatures_Lysate.xlsx'
signature_df = pd.read_excel(path_to_file)
signature_df['delta'] = signature_df['means_group1'] - signature_df['means_group2']


signature_df['type'] = signature_df.up_down.map(regulation_mapping)
signature_df['log(1/p_value)'] = signature_df.p_values.apply(lambda  x:pfunction(x))

# plotly interactive plot
fig = px.scatter(signature_df,
                 x="delta",
                 y="log(1/p_value)",
                 color="type",
                 hover_name="names")
fig.add_hline(y=20)
fig.show()
fig.write_html('/home/amir/Desktop/protocol_effect.html')

# seaborn plot
sns.scatterplot(data=signature_df, x="delta", y="log(1/p_value)", hue="type")
print(len(signature_df))
signature_df_final = signature_df[signature_df['log(1/p_value)'] < 20]
print(len(signature_df_final))
final_list = signature_df_final.names
final_list.to_csv('/home/amir/Desktop/Annika_files/final_list_after_protocol_removal.txt',index=False,header=False)
signature_df_removed = signature_df[signature_df['log(1/p_value)'] > 20]
signature_df_removed = signature_df_removed.names
signature_df_removed.to_csv('/home/amir/Desktop/Annika_files/removed_proteins.txt',index=False,header=False)