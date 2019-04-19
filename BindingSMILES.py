''' reads in each of the 3ish writes from the BindingDB.py file

plus reads in another from url

''' 

import pandas as pd
import random
import json


d = {}
# get Smiles
with open("3820719629006742284.txt", "r") as fp:
    ''' download from pubchem, get smiles

   https://pubchem.ncbi.nlm.nih.gov/pc_fetch/pc_fetch.cgi

   click BROWSE button at the above link and upload "BndDB_ligand_CID.txt" that was written in BindingDB.py 
        BndDB_ligand_CID is a list of chempub ids

    then select SMILES under the FORMAT button

    result will be the txt you want
    ''' 
    line = fp.readline()
    vals = line.strip().split("\t")
    d[int(vals[0])] = vals[1]
    while line:
        line = fp.readline()
        if line.strip() != "":
            vals = line.strip().split("\t")
            d[int(vals[0])] = vals[1]

with open("BndDB_ligand2.json", "r") as fp:
    ''' this file was written in BindingDB.py '''
    df = pd.read_json(fp)

df = df.loc[df["PubChem CID"].isin(d)].drop_duplicates(subset="PubChem CID").dropna()
df["PubChem CID"] = df["PubChem CID"].astype(int)
df["SMILES"] = df["PubChem CID"].apply(lambda x: d.get(x))
df.index = df["id"].astype(int)
df = df.drop(columns=["id"])
print(df.columns)
print(df.shape)

with open("BndDB_protein.json", "r") as fp:
    ''' this file was written in BindingDB.py ''' 
    prot_dict = json.load(fp)

with open("BndDB_seq.json", "r") as fp:
    ''' this file was written in BindingDB.py
    ''' 
    df_seq = pd.read_json(fp)

df_seq = df_seq.drop(columns="index").loc[df_seq["sequence"]!=""].drop_duplicates(subset="sequence").reset_index().drop(columns="index")

cid_length = len(df["PubChem CID"].unique())
cid_list = list(df["PubChem CID"].unique())
random.shuffle(cid_list)
total = 0
with open("lig2seq.csv", "w") as fp:
    ''' This can use pandas 
 
    '''
    fp.write("lig,seq\n")
    for i, ind in enumerate(cid_list):
        if i%10 == 0:
            print(i, cid_length-1, total)
        pchemcid = df.loc[df["PubChem CID"]==ind]
        if len(pchemcid) > 0:
            sequences = pchemcid["sequences"].values[0]
            index = pchemcid.index[0]
            seq_ids = []
            for seq in sequences:
                seq_ids += list(df_seq.loc[df_seq.sequence.apply(lambda x: seq in x)].index.values)
                if len(seq_ids) > 1000:
                    break
            seq_ids = list(set(seq_ids))
            total += len(seq_ids)
            if len(seq_ids) > 0:
                fp.write("\n".join(["{},{}".format(index,seq_id) for seq_id in seq_ids])+"\n")
        if total > 1000000:
            break


# writes the CSVs!
df = df.drop(columns=["PubChem CID", "sequences", "pdb_ids"])
df.to_csv("ligands.csv")
df_seq.to_csv("sequences.csv")
