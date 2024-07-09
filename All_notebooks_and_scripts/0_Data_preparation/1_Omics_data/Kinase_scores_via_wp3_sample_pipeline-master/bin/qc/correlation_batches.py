import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

results_folder = '/media/kusterlab/internal_projects/active/TOPAS/WP31/Playground/Retrospective_study/2022.11.14_CJ_no_kick_of_failed_samples/'

df = pd.read_csv(os.path.join(results_folder, 'preprocessed_fp_with_ref.csv'))
df = df.set_index('Gene names')
# df = df.filter(regex=r'^[H|I]\d')

print(df.head())

# batch 59 sample
#I007-126-1T1-P1
#I034-036-2T1-P1
#I026-012-1M1-P1
#I070-029-1M1-P1
#I070-030-0T2-P1
#I043-011-2M1-P1
#I018-026-1M1-P1
#I034-033-1T1-P1   X

# batch 16
# H021-GDTFYK-M2   X
# H021-FD8YQK-M1
# H021-RYWQDW-T4
# H021-RY13PM-M1
# H021-CEFU8J-T1
# H021-4RWKD3-M3
# H021-X2UVYS-T1
# H021-ATE46U-T1

# batch 20
# H021-YHPFMB-T2
# H021-GDTFYK-M2-R2
# H021-MXRQ8Q-M1-R2
# H021-2TVADY-M1
# H021-LGESLH-T2
# H021-BRW4ZQ-T1
# H021-RYWQDW-T4-R2
# H021-4RWKD3-M3-R2

# batch 21
# I007-007-80690    X
# I007-031-108742-R2
# I062-008-87120    X
# I052-003-93894    X
# I008-006-99146    X
# I137-004-100626   X
# I034-018-101046   X
# I070-014-103148   X

# batch 22
# I002-010-106444   X
# I036-012-107936   X
# I007-029-108166   X
# I024-020-127332   X
# I034-044-186620   X
# I078-024-223164
# I034-052-222996
# I049-015-223128   X

batch_16 = ['H021-GDTFYK-M2', 'H021-FD8YQK-M1', 'H021-RYWQDW-T4', 'H021-RY13PM-M1', 'H021-CEFU8J-T1', 'H021-4RWKD3-M3',
            'H021-X2UVYS-T1', 'H021-ATE46U-T1']
batch_20 = ['H021-YHPFMB-T2', 'H021-GDTFYK-M2-R2', 'H021-MXRQ8Q-M1-R2', 'H021-2TVADY-M1', 'H021-LGESLH-T2', 'H021-BRW4ZQ-T1',
            'H021-RYWQDW-T4-R2', 'H021-4RWKD3-M3-R2']
batch_21 = ['I007-007-80690', 'I007-031-108742-R2', 'I062-008-87120', 'I052-003-93894', 'I008-006-99146', 'I137-004-100626',
            'I034-018-101046', 'I070-014-103148']
batch_22 = ['I002-010-106444', 'I036-012-107936', 'I007-029-108166', 'I024-020-127332', 'I034-044-186620', 'I078-024-223164',
            'I034-052-222996', 'I049-015-223128']
batch_59 = ['I007-126-1T1-P1', 'I034-036-2T1-P1', 'I026-012-1M1-P1', 'I070-029-1M1-P1', 'I070-030-0T2-P1', 'I043-011-2M1-P1',
            'I018-026-1M1-P1', 'I034-033-1T1-P1']
batch_45 = ['H021-YTWLBM-T3-E1', 'H021-G6YVGM-T2-E2', 'H021-G6YVGM-T14-E2', 'H021-9DVBZG-M1-E2', 'H021-EAP7LV-M1-E1', 
            'H021-DXJNNJ-M1-E1', 'H021-DXJNNJ-M2-E1', 'H021-WFU8BQ-M5-R2']
batch_20_ref = ['Reporter intensity corrected 9 Sarcoma_Batch20', 'Reporter intensity corrected 10 Sarcoma_Batch20', 'Reporter intensity corrected 11 Sarcoma_Batch20']

batches = [batch_16, batch_20, batch_21, batch_22, batch_59, batch_45, batch_20_ref]
batch_no = [16, 20, 21, 22, 59, 45, 20.5]
stds = []
# add extra loop for comparison
for i, batch in enumerate(batches):
    df_sub = df[[c for c in df.columns if c in batch]]
    y = df_sub.notna().sum()
    print(df_sub.std())
    stds.append(df_sub.std())
    print(y)
    fig = plt.figure()
    y.plot.bar()
    plt.savefig(os.path.join(results_folder, f'outlier_batch{batch_no[i]}.png'))
    plt.close()

print(stds)
stds = [std.tolist() for std in stds]
stds = [std for sub_std in stds for std in sub_std]
print(stds)
fig = plt.figure()
stds = pd.Series(stds)
stds.plot.bar()
plt.savefig(os.path.join(results_folder, f'outlier_batches_std.png'))
# plot stds


exit()

df_batch20 = df[[c for c in df.columns if c in batch_20]]
df_batch21 = df[[c for c in df.columns if c in batch_21]]
df_batch22 = df[[c for c in df.columns if c in batch_22]]
df_batch59 = df[[c for c in df.columns if c in batch_59]]
df_batch16 = df[[c for c in df.columns if c in batch_16]]
df_batch45 = df[[c for c in df.columns if c in batch_45]]

print(df.shape)
print(df_batch20.shape)

# remove only nan rows
df_batch20 = df_batch20.dropna(how='all', axis=0)
print(df_batch20.shape)

df_batch21 = df_batch21.dropna(how='all', axis=0)
print(df_batch21.shape)

df_batch22 = df_batch22.dropna(how='all', axis=0)
print(df_batch22.shape)

df_batch59 = df_batch59.dropna(how='all', axis=0)
print(df_batch59.shape)

df_batch16 = df_batch16.dropna(how='all', axis=0)
print(df_batch16.shape)

df_batch45 = df_batch45.dropna(how='all', axis=0)
print(df_batch45.shape)
#
fig = plt.figure()
# get overlap venn diagram gene names
b20 = set(df_batch20.index.tolist())
b59 = set(df_batch59.index.tolist())
b20_b59 = list(b20.intersection(b59))
b20_not_b59 = b20.difference(b59)
b59_not_b20 = b59.difference(b20)

print(len(b20_not_b59))
print(len(b59_not_b20))

# Use the venn2 function
venn2(subsets=(len(b20_not_b59), len(b59_not_b20), len(b20_b59)), set_labels=('Group A', 'Group B'))
plt.savefig(os.path.join(results_folder, f'venn_outlier_batch20.png'))
plt.close()

fig = plt.figure()
# get overlap venn diagram gene names
b45 = set(df_batch45.index.tolist())
b59 = set(df_batch59.index.tolist())
b45_b59 = list(b45.intersection(b59))
b45_not_b59 = b45.difference(b59)
b59_not_b45 = b59.difference(b45)

print(len(b45_not_b59))
print(len(b59_not_b45))
# Use the venn2 function
venn2(subsets=(len(b45_not_b59), len(b59_not_b45), len(b45_b59)), set_labels=('Group A', 'Group B'))
plt.savefig(os.path.join(results_folder, f'venn_outlier_batch45.png'))
plt.close()

# # count refs
# print(df.isna().sum())
#
# # correlation
# df_batch21.columns = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
# df_batch45.columns = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
# print(df_batch21.corrwith(df_batch45, axis=0))
