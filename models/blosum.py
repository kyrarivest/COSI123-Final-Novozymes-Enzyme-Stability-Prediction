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



###########BLOSUM#############
def blosum_apply(row):
    if row['type'] == 'SUB':
        return blosum.loc[row['wt'], row['mut']]
    elif row['type'] == 'DEL':
        return -10
    elif row['type'] == 'WT':
        return 0
    else:
        assert False, "Ups"

blosum = pd.read_csv("input/blosum_data/BLOSUM100.txt", sep='\s+', comment='#')
test_df['blosum'] = test_df.apply(blosum_apply, axis=1)
test_df['blosum_rank'] = rankdata(test_df['blosum'])


print(test_df.head())
print(test_df.info())

#Create submission file
blosum_submission = test_df[['seq_id', 'blosum_rank']]
blosum_submission.rename(columns={'blosum_rank':'tm'},inplace = True)
blosum_submission.to_csv('output/blosum_submission.csv', index=False)

print('done')