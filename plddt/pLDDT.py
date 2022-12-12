import pandas as pd
pd.set_option('display.max_columns', 500)
import matplotlib.pyplot as plt
from scipy.stats import rankdata
import numpy as np
import pyrosettacolabsetup
import builtins
from biopandas.pdb import PandasPdb

testDF = pd.read_csv('input/ensemble/EDA-before-tms.csv')
plddt_submission = pd.read_csv('input/novozymes-enzyme-stability-prediction/sample_submission.csv')
# pLDDT scores
# Read AlphaFold2 result for wild type sequence:
plddt  = pd.read_csv('input/novozymes-enzyme-stability-prediction/wildtype_structure_prediction_af2.pdb', sep='\s+', header=None)[[0,5,10]].rename(columns={0:'atom', 5:'resid', 10:'plddt'}).query('atom=="ATOM"').drop_duplicates()

# Add B factor to the testing set:
testDF = pd.merge(
    testDF,
    plddt,
    left_on='resid',
    right_on='resid',
    how='left'
)
testDF['plddt_rank'] = rankdata(-1*testDF['plddt'])
plddt_submission['tm'] = testDF['plddt_rank']
plddt_submission.to_csv('plddt/plddt_submission.csv', index=False)


from biopandas.pdb import PandasPdb


# Wild type result:
# differential pLDDT scores
plddtdiff_submission = pd.read_csv('input/novozymes-enzyme-stability-prediction/sample_submission.csv')

atom_df0 = PandasPdb().read_pdb('input/diff_pLDDT_data/nesp-kvigly-test-mutation-pdbs/WT_unrelaxed_rank_1_model_3.pdb')
atom_df0 = atom_df0.df['ATOM']
wildtype = atom_df0.groupby('residue_number').b_factor.agg('first').values

diffs = []

for index,row in testDF.iterrows():
    aa1 = row.wt
    aa2 = row.mut
    pos = row.resid
    d = 0
    try:
        atom_df = PandasPdb().read_pdb(f'input/diff_pLDDT_data/nesp-kvigly-test-mutation-pdbs/{aa1}{pos}{aa2}_unrelaxed_rank_1_model_3.pdb')
        atom_df = atom_df.df['ATOM']
        mut = atom_df.groupby('residue_number').b_factor.agg('first').values
        d = mut[pos-1] - wildtype[pos-1]
        diffs.append(d)
    except:
        diffs.append(np.nan)

testDF['plddtdiff'] = diffs
testDF['plddtdiff_rank'] = rankdata(testDF['plddtdiff'])
plddtdiff_submission['tm'] = testDF['plddtdiff_rank']
plddtdiff_submission.to_csv('plddt/plddtdiff_submission.csv', index=False)
testDF.to_csv('plddt/plddt_testDF.csv', index=False)



