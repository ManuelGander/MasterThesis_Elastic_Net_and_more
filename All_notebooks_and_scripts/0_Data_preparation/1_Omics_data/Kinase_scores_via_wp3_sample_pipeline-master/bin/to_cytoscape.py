import os.path
import os

import pandas as pd
import numpy as np


location = '/home/cjensen/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.03.11/cytoscape'

print(len(os.listdir(location)))
files = os.listdir(location)
print(files)

for i, file in enumerate(files):
    if '_R1' in file or '_R2' in file or '_R3' in file:
        sample_id = file.split('_')[0] + '_' + file.split('_')[1]
    else:
        sample_id = file.split('_')[0]
    print(sample_id)

    # fp first
    # pairwise check and create list
    if (i % 2) == 0:
        df1 = pd.read_csv(location + '/' + files[i])
        df2 = pd.read_csv(location + '/' + files[i+1])
        in_both = list(set(df1['Protein IDs'].values.tolist()).intersection(set(df2['Gene names'].values.tolist())))
        df1 = df1.rename(columns={"Protein IDs": "Gene names"})
        combined_data = [df1, df2]
        combined_data = pd.concat(combined_data, axis=0, ignore_index=True)
        # from combined data remove those in the intersection
        combined_data.loc[combined_data['Gene names'].isin(in_both), 'mapped'] = 'FPPP'

        # check for duplicates
        combined_data = combined_data.drop_duplicates()

        # we are overwriting now

        # save
        combined_data.to_csv(location + '/' + sample_id + '_for_cytoscape.txt', index=False)
