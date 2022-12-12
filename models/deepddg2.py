import numpy as np 
import pandas as pd 
import random
import requests
from biopandas.pdb import PandasPdb
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from scipy.stats import rankdata
from tqdm.notebook import tqdm
tqdm.pandas()

# Files management
import os

# Ploting
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns

# Libraris for Data process
#import blosum as bl
import Levenshtein
from Levenshtein import distance as levenshtein_distance
from Bio.SubsMat import MatrixInfo

#helper function
def read_list_from_file(list_file):
    with open(list_file) as f:
        lines  = f.readlines()
    return lines


test_df = pd.read_csv("input/ddg_data/novoesp-deepddg-server-predictions-sub-only/test.more.csv")
ddg = read_list_from_file("input/ddg_data/novoesp-deepddg-server-predictions-sub-only/wildtype_structure_prediction_af2.deepddg.ddg.txt")

header = ddg[0]
data = [s.split() for s in ddg[1:]]

df = pd.DataFrame(data, columns = ['chain', 'WT', 'ResID', 'Mut', 'ddG'])
df.ddG = df.ddG.astype(np.float32)
df.ResID = df.ResID.astype(int)  
df.loc[:,'location'] = df.ResID -1  #change to 0-indexing

test_df.location = test_df.location.fillna(-1)
test_df.location = test_df.location.astype(int)

#generate submission csv 
if 1:
    df.loc[:,'mut_string'] = df.WT+df.location.astype(str)+df.Mut
    test_df.loc[:,'mut_string'] =  test_df.wild_type+test_df.location.astype(str)+test_df.mutation
    
    test_df = test_df.merge(df[['ddG','mut_string']], on='mut_string',how='left')
    submit_df = pd.DataFrame({
        'seq_id': test_df.seq_id.values,
        'tm': test_df.ddG.values,
    })
    submit_df.tm = submit_df.tm.fillna(0)
    submit_df.to_csv('output/deepddg2_submission.csv', index=False) #lb0.335



print('done')
