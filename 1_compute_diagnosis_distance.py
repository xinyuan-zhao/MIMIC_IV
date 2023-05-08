# Import necessary libraries
import pandas as pd
import numpy as np
import stad
import igraph
import usedist
import multiprocessing

# Function: Diagnosis similarity
def new_similarity(v1, v2):
    indices = np.where((~np.isnan(v1)) & (~np.isnan(v2)))
    values = np.c_[v1[indices], v2[indices]]
    return np.sum(np.log(1 + 1/np.max(values, axis=1)))

# Function: Standardize distance between 0 and 1
def r01(x):
    return (x - np.min(x, axis=0, keepdims=True)) / (np.max(x, axis=0, keepdims=True) - np.min(x, axis=0, keepdims=True))

# Read ICD9 Codes
# ----------------------------------------------------------------------------------
# Set the working directory
import os
os.chdir("icd_list")
# ----------------------------------------------------------------------------------
icd = pd.read_csv("icd_diagnostic_categories.csv", dtype={'ICD_DIAGNOSTIC_CATEGORY': str})
icd = icd.rename(columns={'ICD_DIAGNOSTIC_CATEGORY': 'ICD', 'SEQUENCE_NUMBER': 'SEQ_NUM'})

## DIAGNOSES: Exploring the ICD9 for patients with patients
# ----------------------------------------------------------------------------------
### Read DIAGNOSES_ICD: List of patients and diagnoses MIMIC-III database
diagnoses = "Use a valid function to read the file e.g. pd.read_csv('...')"
# ---------------------------------------------------------------------------------

# Function: Encapsulate computation process in a function to iterate over multiple ICD codes
def distance_icd(icd):
    # Filter values
    diag = diagnoses.loc[diagnoses['ICD9_CODE'] == icd]

    # Remove duplicates
    diag = diag.groupby('HADM_ID').apply(lambda x: x[x['SEQ_NUM']==x['SEQ_NUM'].min()]).reset_index(drop=True)
    diag = diag.drop('SEQ_NUM', axis=1)

    # Spread diagnosis
    diag_s = diag.groupby(['HADM_ID', 'ICD'])['SEQ_NUM'].min().reset_index()
    diag_s = pd.pivot_table(diag_s, values='SEQ_NUM', index='HADM_ID', columns='ICD', fill_value=np.nan)

    # Compute distance matrix
    diag_s_values = diag_s.values
    new_sim = usedist.dist_make(diag_s_values, new_similarity, "New similarity (custom)")
    distance_matrix =  1 - r01(new_sim)

    # Save numpy array as npy file
    np.save("a_"+icd+"_distance.npy", distance_matrix)

    return("Finished " + icd)

# ----- Parallelize the computation -------
if __name__ == '__main__':
    no_cores = multiprocessing.cpu_count() - 1
    pool = multiprocessing.Pool(processes=no_cores)
    results = [pool.apply_async(distance_icd, args=(icd,)) for icd in icd['ICD'].astype(str)]
    output = [p.get() for p in results]
    pool.close()
    pool.join()
