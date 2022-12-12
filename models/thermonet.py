import math
import multiprocessing
import os
import sys

import Levenshtein
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.model_selection
import tensorflow as tf
from keras import layers, callbacks
from keras import models
from keras import optimizers
from keras.saving.save import load_model
from sklearn.model_selection import GroupKFold
from tqdm import tqdm

from plotly.offline import init_notebook_mode
init_notebook_mode(connected=True)

import plotly.express as px

MULTIPROCESSING = False
BOXSIZE = 16
VOXELSIZE = 1
EPOCHS = 200
N_FOLDS = 10
GROUP_KFOLD = False  # True CV: 39.3, False CV: 39.9
ROSETTA = False
CLEAN_OUTPUT=True
MODELS_PATH = 'models'
DEBUG=False


if DEBUG:
    EPOCHS = 10
    N_FOLDS = 3
    CLEAN_OUTPUT=False


PARAMS = {
    'conv_layer_sizes': (16, 24, 32),
    'dense_layer_size': 24,
    'dropout_rate': 0.5,
    'learning_rate': 0.001,
    'batch_size': 8,
    'scheduler_patience': 10,
    'scheduler_factor': math.sqrt(0.1),
    'early_stopping_patience': 20,
}

try:
    import kaggle_secrets
    
    print('Running in Kaggle')
    WILDTYPE_PDB = '../input/novozymes-enzyme-stability-prediction/wildtype_structure_prediction_af2.pdb'
    PDB_PATH = '../input/thermonet-wildtype-relaxed'
    TRAIN_FEATURES_PATH = '../input/thermonet-features/Q3214.npy'
    TRAIN_TARGETS_PATH = 'data/datasets/Q3214_direct.csv'
    TEST_CSV='../input/novozymes-enzyme-stability-prediction/test.csv'
    TEST_FEATURES_PATH = '../input/thermonet-features/nesp_features.npy'
except Exception as ex:
    print('Running locally')
    WILDTYPE_PDB = 'nesp/thermonet/wildtypeA.pdb'
    PDB_PATH = 'nesp/thermonet/'
    TRAIN_FEATURES_PATH = 'data/datasets/Q3214.npy'
    TRAIN_TARGETS_PATH = 'data/datasets/Q3214_direct.csv'
    TEST_FEATURES_PATH = 'nesp/thermonet/nesp_features.npy'
    TEST_CSV='nesp/test.csv'

def gen_mutations(name, df,
                  wild="VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQ""RVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGT""NAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKAL""GSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK"):
    result = []
    for _, r in df.iterrows():
        ops = Levenshtein.editops(wild, r.protein_sequence)
        assert len(ops) <= 1
        if len(ops) > 0 and ops[0][0] == 'replace':
            idx = ops[0][1]
            result.append([ops[0][0], idx + 1, wild[idx], r.protein_sequence[idx]])
        elif len(ops) == 0:
            result.append(['same', 0, '', ''])
        elif ops[0][0] == 'insert':
            assert False, "Ups"
        elif ops[0][0] == 'delete':
            idx = ops[0][1]
            result.append(['delete', idx + 1, wild[idx], '-'])
        else:
            assert False, "Ups"

    df = pd.concat([df, pd.DataFrame(data=result, columns=['op', 'idx', 'wild', 'mutant'])], axis=1)
    df['mut'] = df[['wild', 'idx', 'mutant']].astype(str).apply(lambda v: ''.join(v), axis=1)
    df['name'] = name
    return df


df_test = gen_mutations('wildtypeA', pd.read_csv(TEST_CSV))
df_train = pd.read_csv(f'{THERMONET_PATH}/{TRAIN_TARGETS_PATH}')

def gen_features(pdb_chain, pos, wt, mt):
    from utils import pdb_utils
    
    pdb_dir = os.path.abspath(PDB_PATH)
    wt_pdb_path = os.path.join(pdb_dir, pdb_chain, pdb_chain + '_relaxed.pdb')
    features_wt = pdb_utils.compute_voxel_features(pos, wt_pdb_path, boxsize=BOXSIZE,
                                                   voxelsize=VOXELSIZE)
    features_wt = np.delete(features_wt, obj=6, axis=0)  # removing 0 'metallic' layer

    mt_pdb_path = os.path.join(pdb_dir, pdb_chain, pdb_chain + '_' + wt + str(pos) + mt + '_relaxed.pdb')
    features_mt = pdb_utils.compute_voxel_features(pos, mt_pdb_path, boxsize=BOXSIZE,
                                                   voxelsize=VOXELSIZE)
    features_mt = np.delete(features_mt, obj=6, axis=0)  # removing 0 'metallic' layer
    features_combined = np.concatenate((features_wt, features_mt), axis=0)
    return features_combined


def thermonet_features(df):
    install_htmd()
    os.environ['HTMD_NONINTERACTIVE'] = '1'
    if MULTIPROCESSING:
        with multiprocessing.Pool() as pool:
            thermonet_features = pool.starmap(gen_features,
                                              [[r['name'], r.idx, r.wild, r.mutant] for _, r in df.iterrows()])
    else:
        thermonet_features = [gen_features(r['name'], r.idx, r.wild, r.mutant) for _, r in
                              tqdm(df.iterrows(), total=len(df))]
    thermonet_features = np.array(thermonet_features)
    return thermonet_features


if not os.path.exists(TRAIN_FEATURES_PATH):    
    np.save(thermonet_features(df_train.rename(columns={'pdb_id': 'name', 'position': 'idx', 'wild_type': 'wild'})), 'train_features.npy')
    TRAIN_FEATURES_PATH = 'train_features.npy'
X = np.load(TRAIN_FEATURES_PATH)
X = np.moveaxis(X, 1, -1)
X.shape

from plotly.offline import init_notebook_mode
init_notebook_mode(connected=True)

def plot_voxels():
    for i in [123, 124, 125, 126]:
        df = pd.DataFrame([(x, y1, z) for x in range(16) for y1 in range(16) for z in range(16)], columns=['x', 'y', 'z'])
        df['occupancy1'] = X[i, :, :, :, 6].flatten() > 0.9
        df['occupancy2'] = X[i, :, :, :, 13].flatten() > 0.9
        df.loc[df.occupancy1 | df.occupancy2, 'color'] = 'blue'
        df.loc[~df.occupancy1 & df.occupancy2, 'color'] = 'red'
        df.loc[df.occupancy1 & ~df.occupancy2, 'color'] = 'green'
        ddg = df_train.ddg[i]
        fig = px.scatter_3d(df.dropna(), x='x', y='y', z='z', color='color', title=f"Train idx:{i}; ddg={ddg}")
        fig.show()

#TRAIN
pdb_ids = df_train.pdb_id
y = df_train.ddg

def gen_model(params):
    def build_model(params):
        conv_layer_sizes = params['conv_layer_sizes']
        dense_layer_size = params['dense_layer_size']
        dropout_rate = params['dropout_rate']
        model = models.Sequential()
        model.add(layers.Conv3D(filters=conv_layer_sizes[0], kernel_size=(3, 3, 3), input_shape=(16, 16, 16, 14)))
        model.add(layers.Activation(activation='relu'))

        for ls in conv_layer_sizes[1:]:
            model.add(layers.Conv3D(filters=ls, kernel_size=(3, 3, 3)))
            model.add(layers.Activation(activation='relu'))

        model.add(layers.MaxPooling3D(pool_size=(2, 2, 2)))
        model.add(layers.Flatten())

        model.add(layers.Dropout(rate=dropout_rate))
        model.add(layers.Dense(units=dense_layer_size, activation='relu'))
        model.add(layers.Dropout(rate=dropout_rate))
        model.add(layers.Dense(units=1))
        return model

    model = build_model(params)
    model.compile(loss='mse', optimizer=tf.keras.optimizers.Adam(
        learning_rate=params['learning_rate'],
        beta_1=0.9,
        beta_2=0.999,
        amsgrad=False
    ), metrics=['mae'])
    return model


def train_model(X_train, y_train, X_val, y_val, params, path):
    model = gen_model(params)
    scheduler = tf.keras.callbacks.ReduceLROnPlateau(
        monitor="val_loss",
        factor=params['scheduler_factor'],
        patience=params['scheduler_patience'],
        verbose=0,
        mode="auto",
        min_delta=0.0001,
        cooldown=0,
        min_lr=1e-5,
    )

    checkpoint = callbacks.ModelCheckpoint(path, monitor='val_loss', verbose=1, save_best_only=True, mode='min')
    early_stopping = callbacks.EarlyStopping(monitor='val_loss', patience=params['early_stopping_patience'])
    result = model.fit(X_train, y_train, validation_data=(X_val, y_val),
                       epochs=EPOCHS, batch_size=params['batch_size'], verbose=1, callbacks=[scheduler, checkpoint, early_stopping])
    return load_model(path), result


!mkdir -p models
kfold = GroupKFold(N_FOLDS)
thermonet_models = []
val_losses = []
if GROUP_KFOLD:
    groups = pdb_ids
else:
    groups = range(len(pdb_ids))

for fold, (train_idx, val_idx) in enumerate(tqdm(kfold.split(X, y, groups=groups), total=N_FOLDS, desc="Folds")):
    X_train = X[train_idx]
    y_train = y[train_idx]
    X_val = X[val_idx]
    y_val = y[val_idx]
    path = f'{MODELS_PATH}/model{fold}.h5'
    model, result = train_model(X_train, y_train, X_val, y_val, PARAMS, path)
    thermonet_models.append(model)
    val_losses.append(result.history['val_loss'])

import glob
thermonet_models = [load_model(f) for f in glob.glob(f'{MODELS_PATH}/model*.h5')]

test_features = np.load(TEST_FEATURES_PATH)
test_features = np.moveaxis(test_features, 1, -1)


test_ddg = np.stack([model.predict(test_features) for model in thermonet_models])
test_ddg = np.mean(test_ddg, axis=0).flatten()

df_test.loc[df_test.op == 'replace', 'ddg'] = -test_ddg
df_test.loc[df_test['op'] == "delete", 'ddg'] = df_test[df_test["op"]=="replace"]["ddg"].quantile(q=0.25)
df_test.loc[df_test['op'] == "same", 'ddg'] = 0.
df_test.rename(columns={'ddg': 'tm'})[['seq_id', 'tm']].to_csv('submission.csv', index=False)