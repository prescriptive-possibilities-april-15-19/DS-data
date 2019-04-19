import json

stopnext = False
mols = []
count = 0

d = {}
d_key = "Structural Model"
d[d_key] = ""
data = []

i = 0
j = 1

with open('BindingDB_All_2D.sdf', 'r', encoding="utf8") as file:
    line = file.readline()
    count = int(line.strip())
    while line:
        line = file.readline()
        if line[:4] == "$$$$":
            n_d = {"id": i, "pdb_ids": []}
            for k in d:
                if d[k].strip() != "":
                    if k in ["PubChem CID", "BindingDB MonomerDB"]:
                        n_d[k] = d[k].strip()
                    elif k == "BindingDB Target Chain Sequence":
                        n_d["sequences"] = d[k].strip().split(" ")
            data.append(n_d)
            stopnext = False
            d = {}
            d_key = "Structural Model"
            d[d_key] = ""
            i += 1
            if i % 1000 == 0:
                print(i, j, count)
        elif line[:3] == "> <":
            if line[3:37] == "Number of Protein Chains in Target":
                d_key = "Number of Protein Chains in Target"
            else:
                d_key = line[3:line.find(">", 3)].strip()
            d[d_key] = ""
        else:
            d[d_key] += line
        j += 1

with open("BndDB_ligand2.json", "w") as fp:
    fp.write(json.dumps(data))

used_monomerid = []
for i, x in enumerate(data):
    if "PubChem CID" in x:
        used_monomerid.append(x["PubChem CID"])
    if i%1000 == 0:
        print(i)
used_monomerid = sorted(list(set(used_monomerid)))

with open("BndDB_ligand_CID.txt", "w") as fp:
    fp.write(",".join(used_monomerid))

data_pdb = {}
data_seq = []
with open("pdb_seqres.txt", "r") as fp:
    line = fp.readline()
    v1 = ""
    v2 = ""
    val = ""
    skip = False
    i = 0
    if line[0] == ">":
        if line[1:].split(" ")[1] != "mol:protein":
            skip = True
        else:
            skip = False
            v = line[1:].split(" ")[0]
            v1 = v.split("_")[0]
            v2 = v.split("_")[-1]
    elif not skip:
        val += line.strip()
    while line:
        line = fp.readline()
        if line != "":
            if line[0] == ">":
                if v1 not in data_pdb:
                    data_pdb[v1] = {}
                val_ind = 0
                if val not in data_seq:
                    val_ind = len(data_seq)
                    data_seq.append(val)
                else:
                    val_ind = data_seq.index(val)
                data_pdb[v1][v2] = val_ind
                val = ""
                if line[1:].split(" ")[1] != "mol:protein":
                    skip = True
                else:
                    skip = False
                    v = line[1:].split(" ")[0]
                    v1 = v.split("_")[0]
                    v2 = v.split("_")[-1]
            elif not skip:
                val += line.strip()
            i += 1
            if i % 1000 == 0:
                print(i)

with open("BndDB_protein.json", "w") as fp:
    fp.write(json.dumps(data_pdb))
with open("BndDB_seq.json", "w") as fp:
    fp.write(json.dumps([{"index": i, "sequence": x} for i, x in enumerate(data_seq)]))
