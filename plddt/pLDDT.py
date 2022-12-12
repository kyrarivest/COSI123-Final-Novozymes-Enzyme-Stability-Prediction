import pandas as pd
pd.set_option('display.max_columns', 500)
import matplotlib.pyplot as plt
from scipy.stats import rankdata
import numpy as np
import pyrosettacolabsetup
import builtins
from biopandas.pdb import PandasPdb


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

testDF = pd.concat([test_df, pd.DataFrame(data=result, columns=['type', 'resid', 'wt', 'mut'])], axis=1)

#testDF = pd.read_csv('input/ensemble/EDA-before-tms.csv')
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



