import numpy as np
import pandas as pd
import csv
from tqdm.auto import tqdm
from scipy.stats import rankdata
import Levenshtein
import matplotlib.pyplot as plt
import seaborn as sns

# Wild type sequence provided in the "Dataset Description":
wt = 'VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQRVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGTNAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKALGSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK'

# Read testing set sequences and pH:
test_df = pd.read_csv("input/novozymes-enzyme-stability-prediction/test.csv")


# Add mutation information to testing set:
result = []
for _, row in test_df.iterrows():
    ops = Levenshtein.editops(wt, row['protein_sequence'])
    #print("n")
    #print(ops)
    assert len(ops) <= 1
    if len(ops) > 0 and ops[0][0] == 'replace':
        idx = ops[0][1]
        result.append(['SUB', idx + 1, wt[idx], row['protein_sequence'][idx]])
    elif len(ops) == 0:
        result.append(['WT', 0, '', ''])
    elif ops[0][0] == 'insert':
        assert False, "Ups"
    elif ops[0][0] == 'delete':
        idx = ops[0][1]
        result.append(['DEL', idx + 1, wt[idx], '_'])
    else:
        assert False, "Ups"

test_df = pd.concat([test_df, pd.DataFrame(data=result, columns=['type', 'resid', 'wt', 'mut'])], axis=1)

# Read VMD/NAMD output:
namd = pd.read_csv('input/RMSD_data/novozymes-md2/residue_rmsd_sasa_last.dat', sep='\t', header=None, names=['resid','rmsd','sasa0','sasaf'])


# Add VMD/NAMD results to the testing set:
test_df = pd.merge(
    test_df,
    namd[['resid','rmsd']],
    left_on='resid',
    right_on='resid',
    how='left'
)

test_df.loc[test_df['type']=='WT','rmsd'] = test_df['rmsd'].dropna().max()
# test_df.loc[test_df['type']=='WT','sasaf'] = test_df['sasaf'].dropna().max()

test_df['rmsd_rank'] = rankdata(test_df['rmsd'])

print(test_df.head())
print(test_df.info())

#Create submission file
rmsd_submission = test_df[['seq_id', 'rmsd_rank']]
rmsd_submission.rename(columns={'rmsd_rank':'tm'},inplace = True)
rmsd_submission.to_csv('output/rmsd_submission.csv', index=False)


print('done')