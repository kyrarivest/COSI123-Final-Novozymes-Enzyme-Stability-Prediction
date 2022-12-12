import os
import pandas as pd
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import pyrosettacolabsetup;
import builtins
import numpy as np
import pandas as pd
from scipy.stats import rankdata

curr_dir = os.getcwd()

relaxed_path = 'input/rosetta_data/wildtypeA'
files = os.listdir(relaxed_path)
relaxed_pdbs = [file for file in files if file.endswith(".pdb")]
mutation_to_pdb = {}
mutation_to_pdb['mutation'] = [file.split('_')[1] for file in relaxed_pdbs]
mutation_to_pdb['path'] = [file for file in relaxed_pdbs]
mutation_to_pdb_df = pd.DataFrame(mutation_to_pdb)


"""
#pyrosettacolabsetup.install_pyrosetta(cache_wheel_on_google_drive=False)
import pyrosetta
pyrosetta.init()
from pyrosetta.teaching import *

scores = []
sfxn = get_score_function(True)
for i in range(len(mutation_to_pdb_df)):
    if not i%100: print(i)
    pose = pyrosetta.pose_from_pdb(os.path.join(relaxed_path, mutation_to_pdb_df.iloc[i]['path']))
    scores.append(sfxn(pose))

mutation_to_pdb_df['scores'] = scores
scores = mutation_to_pdb_df

print("SCORES: ")
print(scores.head())

scores.to_csv('rosetta_scores.csv', index=False)
print('done')
"""

scores = pd.read_csv('input/rosetta_data/rosetta_scores.csv')

def find_mut(row):
    mut = row.mutant_seq
    seq = row.sequence
    same = True
    for i,(x,y) in enumerate(zip(seq,mut)):
        if x!=y: 
            same = False
            break
    row['WT'] = seq[i]
    row['position'] = i+1
    if not same:
        if len(mut) < len(seq):
            row['MUT'] = 'X'
        else:
            row['MUT'] = mut[i]
    else: 
        row['position'] = -1
        row['MUT'] = 'X'
    row['mutation_key'] = row['WT']+str(row['position']) + row['MUT']
    return row

novo_test = pd.read_csv("input/novozymes-enzyme-stability-prediction/test.csv")
novo_test = novo_test.rename({'protein_sequence': 'mutant_seq', 'seq_id': 'source_df_id'}, axis = 1)
novo_test['sequence'] = 'VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQRVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGTNAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKALGSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK'
novo_test = novo_test.apply(find_mut,axis=1)
novo_test = novo_test.join(scores.set_index('mutation'), on='mutation_key')
novo_test['scores'] = -novo_test['scores']
novo_test.loc[novo_test['scores'].isna(), 'scores'] = novo_test.loc[~novo_test['scores'].isna()].quantile(q=0.25)['scores']
novo_test['scores_rank'] = rankdata(novo_test['scores'])

print("NOVO TEST")
print(novo_test.head())

submission_rosetta_scores = novo_test[['source_df_id','scores_rank']]
submission_rosetta_scores = submission_rosetta_scores.rename({'source_df_id': 'seq_id', 'scores_rank': 'tm'}, axis = 1)
submission_rosetta_scores.to_csv(curr_dir + '/output/rosetta_submission.csv', index=False)

print("submission file: ")
print(submission_rosetta_scores.head())

