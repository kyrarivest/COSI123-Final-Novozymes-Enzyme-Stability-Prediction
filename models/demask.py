import pandas as pd
pd.set_option('display.max_columns', 500)
import matplotlib.pyplot as plt
from scipy.stats import rankdata
import numpy as np
import pyrosettacolabsetup
import builtins
import Levenshtein
from biopandas.pdb import PandasPdb
import demask

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

testDF = pd.concat([test_df, pd.DataFrame(data=result, columns=['type', 'resid', 'wt', 'mut'])], axis=1)

demask_submission = pd.read_csv('input/novozymes-enzyme-stability-prediction/sample_submission.csv')


# Run DeMaSk on https://demask.princeton.edu/query/ using the wildtype sequence as parameter
demask = pd.read_csv('input/demask_data/demask_prediction.txt', sep='\t', usecols=[0,1,2,3], names=['resid','wt','mut','demask'], skiprows=1)

# Merge the DeMaSk output with the test data
testDF = pd.merge(
    testDF.set_index(['wt','resid','mut']),
    demask.set_index(['wt','resid','mut']),
    left_index=True,
    right_index=True,
    how='left'
).reset_index()

print(testDF)

# Set the DeMaSk score for wildtype and deletion sequences to 0
testDF.loc[testDF['type']=='WT','demask'] = 0
testDF.loc[testDF['type']=='DEL','demask'] = testDF['demask'].dropna().min()

# Compute the rank of each DeMaSk score
testDF['demask_rank'] = rankdata(testDF['demask'])
demask_submission['tm'] = testDF['demask_rank']
demask_submission.to_csv('output/demask_submission.csv', index=False)

print('done')