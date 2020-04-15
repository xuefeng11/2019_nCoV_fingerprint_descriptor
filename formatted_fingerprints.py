from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors
import pandas as pd
import numpy as np
from rdkit.Chem import MACCSkeys



import pybel
import pandas as pd


#input smile file should use tab as seperator and have the followng format
#smile  drug_id
#OC(=O)C(=O)C	DB-17
_drug_db_path="../ena+db/adrp_nsp10_nsp15_plpro_cov_dump.csv"
df = pd.read_csv(_drug_db_path, sep=',',dtype={"drug_name":str,"smiles":str,"dock_score":float,"smiles_dbase":str,"receptor":str})

PLPro_pocket6_df = df[df["receptor"] == "PLPro_pocket6"]

print("raws to process:",PLPro_pocket6_df.shape)

_canonical_smile_arr=[]
_drug_id_arr=[]
_dock_score_arr=[]
_smile_arr=[]
_ecfp2_arr=[]
_ecfp4_arr=[]
_ecfp6_arr=[]
_ecfp2_fingers_2048_arr=[]
_ecfp4_fingers_2048_arr=[]
_ecfp6_fingers_2048_arr=[]


_maccs_key_arr=[]

nbits = 1024

equal_counter=0
ecfp_1024_arr=[]
ecfp_2048_arr=[]
counter=0
for index, row in PLPro_pocket6_df.iterrows():
    try:
        canonical_smile=pybel.readstring("smi", row["smiles"]).write("can").strip()

    except Exception as e:
        print("warning fail in processing(skip) smile:",row["smiles"],", error msg:",e)
        continue

    try:
        maccs_key=np.array(MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(canonical_smile))).tolist()

    except Exception as e:
        print("warning fail in processing(skip) maccs:", row["smiles"], ", error msg:", e)
        continue

    try:
        m1 = Chem.MolFromSmiles(canonical_smile)
        # ecfp2

        ecfp2 = AllChem.GetMorganFingerprintAsBitVect(m1, 1, nBits=nbits)
        ecfp2_fingers = np.array(ecfp2).tolist()
        ecfp2_fingers_2048 = np.array(AllChem.GetMorganFingerprintAsBitVect(m1, 1, nBits=2048)).tolist()

        # ecfp4
        ecfp4 = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=nbits)
        ecfp4_fingers = np.array(ecfp4).tolist()
        ecfp4_fingers_2048 = np.array(AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=2048)).tolist()

        # ecfp6
        ecfp6 = AllChem.GetMorganFingerprintAsBitVect(m1, 3, nBits=nbits)
        ecfp6_fingers = np.array(ecfp6).tolist()
        ecfp6_fingers_2048 = np.array(AllChem.GetMorganFingerprintAsBitVect(m1, 3, nBits=2048)).tolist()


        np.savetxt("ecfp4_1024.txt",ecfp4_fingers,fmt="%s")
        temp1=np.sum(ecfp4_fingers)
        np.savetxt("ecfp4_2048.txt",ecfp4_fingers_2048,fmt="%s")
        temp2=np.sum(ecfp4_fingers_2048)

    except Exception as e:
        print("warning fail in processing(skip) edfp:", row["smiles"], ", error msg:", e)
        continue

    _canonical_smile_arr.append(canonical_smile)
    _drug_id_arr.append(row["drug_name"])
    _smile_arr.append(row["smiles"])
    _dock_score_arr.append(row["dock_score"])
    _maccs_key_arr.append(maccs_key)
    _ecfp2_arr.append(ecfp2_fingers)
    _ecfp4_arr.append(ecfp4_fingers)
    _ecfp6_arr.append(ecfp6_fingers)
    _ecfp2_fingers_2048_arr.append(ecfp2_fingers_2048)
    _ecfp4_fingers_2048_arr.append(ecfp4_fingers_2048)
    _ecfp6_fingers_2048_arr.append(ecfp6_fingers_2048)

    counter=counter+1

    if temp1==temp2:
        equal_counter = equal_counter + 1

    ecfp_1024_arr.append(temp1)
    ecfp_2048_arr.append(temp2)


#    if counter >2:
#        break

    if counter%1000==0:
        print("processing counter",counter)

print("counter:",counter)
print("ecfp_1024_mean:",np.mean(ecfp_1024_arr))
print("ecfp_2048_mean:",np.mean(ecfp_2048_arr))




data_format = {"canonical_smile":_canonical_smile_arr,"ID": _drug_id_arr,"smile":_smile_arr,"dock_score":_dock_score_arr,"maccs_key":_maccs_key_arr
               ,"ecfp2":_ecfp2_arr,"ecfp4":_ecfp4_arr,"ecfp6":_ecfp6_arr,"ecfp2_2048":_ecfp2_fingers_2048_arr,"ecfp4_2048":_ecfp4_fingers_2048_arr,"ecfp6_2048":_ecfp6_fingers_2048_arr}

df_out = pd.DataFrame(data_format)
df_out.to_csv("ena+db_fingerprints.csv", sep='\t', index=False)


nbits = 1024
canonical_smile=pybel.readstring("smi", row["smiles"]).write("can").strip()
m1 = Chem.MolFromSmiles(canonical_smile)
ecfp4 = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=nbits)
ecfp4_fingers = np.array(ecfp4).tolist()
