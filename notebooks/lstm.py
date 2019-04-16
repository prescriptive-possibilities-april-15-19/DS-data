import pandas as pd
import numpy as np
import pickle

from keras.models import Model
from keras.layers import Input, Dense, Embedding, Concatenate
from keras.layers import LSTM

seed = 42
np.random.seed(seed)
num_samples = 20
lstm_size = 128
dropout = 0.2
rec_dropout = 0.2
epochs = 20
batch_size = 128
pad_length = 1000
pad_length_2 = 5000


def sml_vec(pad_length, smiles):
    vec = np.zeros((pad_length, 128))
    ind = [ord(x) for x in smiles]
    for i in range(len(ind)):
        vec[i, ind[i]] = 1.0
    return vec


with open("../models/tfidf.pickle", "rb") as fp:
    tfidf = pickle.load(fp)
    seq_vec_length = len(tfidf.vocabulary_)
with open("../models/PUC.pickle", "rb") as fp:
    PUC_models = pickle.load(fp)
for k in PUC_models:
    PUC_models[k] = pickle.loads(PUC_models[k])
print(seq_vec_length)


df_ligands = pd.read_csv("../data/ligands.csv", index_col="id", usecols=["id", "SMILES"])
df_sequences = pd.read_csv("../data/sequences.csv", index_col=0)
df_bnd = pd.read_csv("../data/lig2seq.csv")

PUC_models_l = [x for x in PUC_models]
df_bnd = df_bnd.loc[~df_bnd["lig"].isin(PUC_models_l)].sample(n=1000, random_state=seed)

X_test2_a = np.array([x for x in df_bnd["lig"].apply(lambda x: sml_vec(pad_length, df_ligands.loc[df_ligands.index==x,"SMILES"].values[0])).values])
X_test2_b = np.array([x for x in df_bnd["seq"].apply(lambda x: sml_vec(pad_length_2, df_sequences.loc[df_sequences.index==x,"sequence"].values[0])).values])
print(X_test2_a.shape)

df_sequences["seq_len"] = df_sequences["sequence"].apply(lambda x: len(x))
df_sequences = df_sequences.loc[df_sequences["seq_len"]<=pad_length_2]



df_total = pd.DataFrame()
for i, k in enumerate(PUC_models):
    print(i)
    df_temp = pd.DataFrame()
    bc_model = PUC_models[k]
    smiles = df_ligands.loc[i, "SMILES"]
    if len(smiles) < pad_length:
        seq_samp = df_sequences.sample(n=num_samples, random_state=seed+i)
        df_temp["seq_id"] = seq_samp.index.values
        df_temp["seq_sequence"] = seq_samp["sequence"].values
        df_temp["seq_vec"] = df_temp["seq_sequence"].apply(lambda x: sml_vec(pad_length_2, x))
#         df_temp["seq_vec"] = df_temp["seq_vec"].apply(lambda x: np.repeat(x.reshape((1,-1)), pad_length, axis=0))
        df_temp["lig_id"] = i
        df_temp["lig_SMILES"] = smiles
        smlvec = sml_vec(pad_length, smiles)
        df_temp["lig_vec"] = 0
        df_temp["lig_vec"] = df_temp["lig_vec"].apply(lambda x: smlvec)
        df_temp["pred_binding"] =df_temp["seq_sequence"].apply(lambda x: bc_model.predict([tfidf.transform([x])[0].toarray()[0]]))
        df_total = pd.concat([df_total, df_temp], ignore_index=True)

# df_total["t_vec"] = 0
# for i in df_total.index:
#     lig_vec = df_total.loc[i,"lig_vec"]
#     seq_vec = df_total.loc[i,"seq_vec"]
#     df_total.loc[i, ["t_vec"]] = np.concatenate((lig_vec, seq_vec), axis=1)


indices = df_total.index.values

df_test = df_total.loc[df_total.index.isin(indices[:num_samples*len(PUC_models)//5])]
df_fit = df_total.loc[df_total.index.isin(indices[num_samples*len(PUC_models)//5:])]
df_train = df_fit.sample(frac=0.8, random_state=seed)
df_val  = df_fit.loc[~df_fit.index.isin(df_train.index)]

X_train_a = np.array([x for x in df_train["lig_vec"].values])
X_train_b = np.array([x for x in df_train["seq_vec"].values])
X_val_a = np.array([x for x in df_val["lig_vec"].values])
X_val_b = np.array([x for x in df_val["seq_vec"].values])
X_test_a = np.array([x for x in df_test["lig_vec"].values])
X_test_b = np.array([x for x in df_test["seq_vec"].values])
y_train = df_train["pred_binding"].values
y_val = df_val["pred_binding"].values
y_test = df_test["pred_binding"].values


inp_a = Input(shape=(pad_length, 128,))
x_a = LSTM(lstm_size, dropout=dropout, recurrent_dropout=rec_dropout)(inp_a)
inp_b = Input(shape=(pad_length_2, 128,))
x_b = LSTM(lstm_size, dropout=dropout, recurrent_dropout=rec_dropout)(inp_a)
x = Concatenate()([x_a,x_b])
x = Dense(256, activation='relu')(x)
x = Dense(1, activation='sigmoid')(x)
model = Model(inputs=[inp_a,inp_b], outputs=[x])
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
model.fit([X_train_a,X_train_b], y_train, batch_size=batch_size, epochs=epochs, validation_data=([X_val_a,X_val_b], y_val), verbose=1)


with open("../models/lstm.pickle", "wb") as fp:
    pickle.dump(model, fp)


score, acc = model.evaluate([X_test_a,X_test_b], y_test, batch_size=batch_size)
print('Test score:', score)
print('Test accuracy:', acc)


score = model.predict([X_test2_a,X_test2_b]).mean()
print('Test score:', score)

print(model.predict([X_test2_a,X_test2_b]))
