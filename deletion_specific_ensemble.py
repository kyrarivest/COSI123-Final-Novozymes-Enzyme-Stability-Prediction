import numpy as np
import pandas as pd
import csv
import os
from tqdm.auto import tqdm
from scipy.stats import rankdata
import Levenshtein
import matplotlib.pyplot as plt
import seaborn as sns

two_colors = sns.xkcd_palette(['red', 'bright blue'])


# Plot rank distributions
def plot_rank_dist(name, ax, show_del=False):

    sns.kdeplot(
        data=test_df.query('type=="SUB"'),
        x='{}_rank'.format(name),
        bw_adjust=0.3,
        lw=3,
        label='SUB',
        ax=ax,
        color='k'
    )

    ax.vlines(
        test_df.query('type=="DEL"')['{}_rank'.format(name)],
        ax.get_ylim()[0],
        ax.get_ylim()[1],
        lw=5,
        label='DEL',
        color=two_colors[0]
    )

    ax.vlines(
        test_df.query('type=="WT"')['{}_rank'.format(name)],
        ax.get_ylim()[0],
        ax.get_ylim()[1],
        lw=5,
        label='WT',
        color=two_colors[1]
    )

    if show_del:
        sns.kdeplot(
            data=test_df.query('type=="DEL"'),
            x='{}_rank'.format(name),
            bw_adjust=0.3,
            lw=3,
            label='DEL',
            ax=ax,
            color=two_colors[0]
        )

        ax.vlines(
            test_df.query('type=="DEL"')['{}_rank'.format(name)],
            ax.get_ylim()[0],
            ax.get_ylim()[1],
            lw=5,
            label='DEL',
            color=two_colors[0]
        )

    ax.set_xlim(-50,2550)
    ax.set_title('{} rank distribution'.format(name), fontsize=20)
    ax.set_xlabel('{}_rank'.format(name), fontsize=20)
    ax.set_ylabel('Density', fontsize=20)

    ax.tick_params(labelsize=12)
    ax.legend(loc=1)

    return ax



curr_dir = os.getcwd()

# Wild type sequence provided in the "Dataset Description":
wt = 'VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQRVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGTNAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKALGSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK'

# Read testing set sequences and pH:
test_df = pd.read_csv(curr_dir + "/input/novozymes-enzyme-stability-prediction/test.csv")


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

#print(test_df.head())
#print(test_df.info())


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

blosum = pd.read_csv(curr_dir + '/input/blosum_data/BLOSUM100.txt', sep='\s+', comment='#')
test_df['blosum'] = test_df.apply(blosum_apply, axis=1)
test_df['blosum_rank'] = rankdata(test_df['blosum'])

#print(blosum.head())
#print(blosum.info())
print(test_df.head())
print(test_df.info())




###########PLDDT#############
# Read AlphaFold2 result for wild type sequence:
plddt = (
    pd.read_csv(curr_dir + '/input/novozymes-enzyme-stability-prediction/wildtype_structure_prediction_af2.pdb', sep='\s+', header=None)[[0,5,10]]
    .rename(columns={0:'atom', 5:'resid', 10:'plddt'})
    .query('atom=="ATOM"')
    .drop_duplicates()
)


# Add B factor to the testing set:
test_df = pd.merge(
    test_df,
    plddt,
    left_on='resid',
    right_on='resid',
    how='left'
)



test_df['plddt_rank'] = rankdata(-1*test_df['plddt'])

########### DIFFERNTIAL PLDDT #############
plddtdiff = []

# Wild type result:
wt_plddt = (
    pd.read_csv(curr_dir + '/input/diff_pLDDT_data/nesp-kvigly-test-mutation-pdbs/WT_unrelaxed_rank_1_model_3.pdb', sep='\s+')
    .loc['ATOM'].reset_index()
    .loc[:, ['level_4', 'MODEL']].drop_duplicates()
    .rename(columns={'level_4':'resid', 'MODEL':'plddt'})
    .astype({'resid':int})
    .set_index('resid')
)

# Add difference in pLDDTto the testing set:>
for _,row in tqdm(test_df.iterrows(), total=test_df.shape[0]):
    file_path = curr_dir + '/input/diff_plddt_data/nesp-kvigly-test-mutation-pdbs/{}{}{}_unrelaxed_rank_1_model_3.pdb'.format(row['wt'], row['resid'], row['mut'])
    if os.path.exists(file_path):
        tdf = (
            pd.read_csv(file_path, sep='\s+')
            .loc['ATOM'].reset_index()
            .loc[:, ['level_4', 'MODEL']].drop_duplicates()
            .rename(columns={'level_4':'resid', 'MODEL':'plddt'})
            .astype({'resid':int})
            .set_index('resid')
        )
        plddtdiff.append((tdf.loc[row['resid']] - wt_plddt.loc[row['resid']]).values[0])
    else:
        plddtdiff.append(np.nan)

test_df['plddtdiff'] = plddtdiff
test_df['plddtdiff_rank'] = rankdata(test_df['plddtdiff'])



########### DEEP DDG #############
print('########### DEEP DDG #############')

#had to reproduce model (I dont think the out file is exactly the same, but I got the ddg_rank somehow)
ddg = pd.read_csv(curr_dir + '/output/ddgout.csv')
ddg_rank = ddg['ddg_rank']
test_df.merge(ddg_rank)

#Check to make sure the order of rows is the same:
"""
ddg_rank = ddg['seq_id']
blosum_rank = test_df['seq_id']

for i in range(len(blosum_rank)):
    if blosum_rank[i] != ddg_rank[i]:
        print("!!! " +str(i))
        break

print('done')
"""

########### DEMASK #############
print('########### DEMASK #############')

demask = pd.read_csv(curr_dir + '/input/demask_data/demaskout.txt', sep='\t', usecols=[0,1,2,3], names=['resid','wt','mut','demask'], skiprows=1)
print(demask)

# Add DeMask output to the testing set:
test_df = pd.merge(
    test_df.set_index(['wt','resid','mut']),
    demask.set_index(['wt','resid','mut']),
    left_index=True,
    right_index=True,
    how='left'
).reset_index()

test_df.loc[test_df['type']=='WT','demask'] = 0
test_df.loc[test_df['type']=='DEL','demask'] = test_df['demask'].dropna().min()


test_df['demask_rank'] = rankdata(test_df['demask'])


########### RMSD #############
print('########### RMSD #############')

# Read VMD/NAMD output:
namd = pd.read_csv(curr_dir + '/input/RMSD_data/novozymes-md2/residue_rmsd_sasa_last.dat', sep='\t', header=None, names=['resid','rmsd','sasa0','sasaf'])

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


########### SASA #############
print('########### SASA #############')

# Read VMD/NAMD output:
namd = pd.read_csv(curr_dir + '/input/SASA_data/novozymes-md/residue_rmsd_sasa_last.dat', sep='\t', header=None, names=['resid','rmsd','sasa0','sasaf'])

# Add VMD/NAMD results to the testing set:
test_df = pd.merge(
    test_df,
    namd[['resid','sasaf']],
    left_on='resid',
    right_on='resid',
    how='left'
)

# test_df.loc[test_df['type']=='WT','rmsd'] = test_df['rmsd'].dropna().max()
test_df.loc[test_df['type']=='WT','sasaf'] = test_df['sasaf'].dropna().max()
test_df['sasaf_rank'] = rankdata(test_df['sasaf'])


########### ROSETTA #############
print('########### ROSETTA #############')




