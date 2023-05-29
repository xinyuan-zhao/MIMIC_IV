import math
import pandas as pd

def similarity(patient1, patient2):
    '''
    Calculate the similarity between the diagnosis codes of two patient encounters. 
    Weights contribution of ICD code similarity by line order.

    See Alcaide et. al.
    '''
    similarity = 0

    for ix1, dx1 in enumerate(patient1):
        if dx1 not in patient2:
            continue

        ix2 = patient2.index(dx1)
        similarity += math.log( 1 + 1/max(ix1+1, ix2+1) )

    return similarity
def convert_similarity_to_distance(possible_edges):
    '''
    The similarity measure S needs to be converted into a distance measure D = 1 − Snormalized where Snormalized = S/max(S).
    :return: distance
    '''
    # DataFrame named 'possible_edges' with columns 'hadm1', 'hadm2', 'similarity'

    # Calculate Snormalized

    print(possible_edges['similarity'].max())
    possible_edges['s_normalized'] = possible_edges['similarity'] / possible_edges['similarity'].max()

    possible_edges['distance'] = 1 - possible_edges['s_normalized']

    result = possible_edges[['hadm1', 'hadm2', 'distance']]
    return result


def make_list_of_possible_edges(data, condition):
    '''
    input data is DIAGNOSES_ICD.csv from MIMIC-III
    '''
    
    from tqdm import tqdm # pip3 install tqdm

    # filter to patients with condition of interest
    data.HADM_ID = data.HADM_ID.astype(str)
    hadm_list = list(data[data.ICD9_CODE == condition]['HADM_ID'])
    data = data[data.HADM_ID.isin(hadm_list)]

    # create a dictionary of patients and their diagnoses
    codes = { hadm: list(data[data.HADM_ID == hadm].sort_values(by='SEQ_NUM')['ICD9_CODE'])
              for hadm 
              in data.HADM_ID.unique() }

    # list of possible edges is similarity of all unique pairs of patients
    possible_edges = []
    for i,hadm1 in tqdm(enumerate(hadm_list), total=len(hadm_list)):
        for j,hadm2 in enumerate(hadm_list):
            if i < j: # only need to calculate similarity once for each pair
                possible_edges.append( [hadm1, hadm2, similarity(codes[hadm1], codes[hadm2])] )

    return pd.DataFrame(possible_edges, columns=['hadm1', 'hadm2', 'similarity'])

def find_parent(node, parent):
    if parent[node] != node:
        parent[node] = find_parent(parent[node], parent)
    return parent[node]

def union(u, v, parent):
    parent[find_parent(v, parent)] = find_parent(u, parent)

def create_minimum_spanning_tree(possible_edges_distance):
    # Drop rows with missing or inconsistent values
    possible_edges_distance = possible_edges_distance.dropna()
    possible_edges_distance = possible_edges_distance.astype(int)

    sorted_edges = possible_edges_distance.sort_values(by='distance')
    parent = {}
    unique_nodes = set(possible_edges_distance['hadm1'].unique()) | set(possible_edges_distance['hadm2'].unique())

    for node in unique_nodes:
        parent[node] = node

    minimum_spanning_tree = []

    for _, row in sorted_edges.iterrows():
        u, v, distance = row['hadm1'], row['hadm2'], row['distance']
        if u in parent and v in parent and find_parent(u, parent) != find_parent(v, parent):
            minimum_spanning_tree.append((u, v, distance))
            union(u, v, parent)

    return minimum_spanning_tree



if __name__ == "__main__":
    # test functions using example from Alcaide et. al. 
    # Patient A: '115057'
    # Patient B: '117154'
    # similarity: 0.56

    # test similarity function
    patientA = ['99662', '99591', '5990', '4019']
    patientB = ['4329', '43491', '99702', '99591', '5990', '4019']
    test1 = round(similarity(patientA, patientB), 2)

    # test make_list_of_possible_edges function
    df = pd.read_csv('./DIAGNOSES_ICD.csv')
    possible_edges = make_list_of_possible_edges(df, condition='99591') # 99591 is sepsis
    # print(possible_edges)
    test2 = round(possible_edges[(possible_edges.hadm1 == '115057') & (possible_edges.hadm2 == '117154') |
                                 (possible_edges.hadm1 == '117154') & (possible_edges.hadm2 == '115057')].
                  similarity.values[0], 2)

    possible_edges_distance = convert_similarity_to_distance(possible_edges)
    print(convert_similarity_to_distance(possible_edges))

    # create minimum spanning tree
    minimum_spanning_tree = create_minimum_spanning_tree(possible_edges_distance)
    for edge in minimum_spanning_tree:
        print(edge)


    print(f'Test 1: {test1} | passed: {test1==0.56}' )
    print(f'Test 2: {test2} | passed: {test2==0.56}' )
