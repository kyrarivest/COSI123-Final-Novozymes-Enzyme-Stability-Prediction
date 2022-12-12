import pandas as pd
pd.set_option('display.max_columns', 500)
import matplotlib.pyplot as plt
from scipy.stats import rankdata
import numpy as np
import pyrosettacolabsetup
import builtins
from biopandas.pdb import PandasPdb
import demask

testDF = pd.read_csv('input/ensemble/EDA-before-tms.csv')
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
demask_submission.to_csv('demask/demask_submission.csv', index=False)
testDF.to_csv('demask/demask_testDF.csv', index=False)