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


print(test_df.head())
print(test_df.info())

#Create submission file
blosum_submission = test_df[['seq_id', 'blosum_rank']]
blosum_submission.rename(columns={'blosum_rank':'tm'},inplace = True)
blosum_submission.to_csv(curr_dir + '/output/blosum_submission.csv', index=False)

print('done')