import math
import csv

def read_csv(file):
    # need 'hadm_id', 'sequence_number'，'icd_diagnostic_category'
    # ['subject_id', 'hadm_id', 'sequence_number', 'level', 'icd_diagnostic_category']

    # ['2', '163353', '1', '2', 'V30']
    # file = 'file.csv'
    # Open the CSV file

    with open(file, 'r') as file:

        # Create a reader object
        csv_reader = csv.reader(file)
        # Counter variable to keep track of number of rows printed
        row_count = 0

        # Loop through each row in the CSV file
        for row in csv_reader:
            # Print each row
            print(row)
            # Increment the row count
            row_count += 1

            # Exit the loop once you've printed the first 10 rows
            if row_count == 10:
                break

def distance(patient1, patient2):
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

def make_distance_matrix(patients):
    '''
    codes in order by SEQ_NUM
    input: [ {HADM_ID: [code1, code2, ... ,coden]}, ... ]
    output: [ [HADM_ID1, HADM_ID2, distance], ...]
    '''

    # matrix = [patient1, patient2, distance]
    print(patients)
    matrix = []
    for patient1 in patients:
        for patient2 in patients:
            # print(list(patient1.keys())[0], list(patient2.keys())[0])
            # print(list(patient1.values())[0], list(patient2.values())[0])
            d = distance(list(patient1.values())[0], list(patient2.values())[0])
            # TODO: add this to the matrix
            matrix.append([list(patient1.keys())[0], list(patient2.keys())[0], d])
    print(matrix)


def make_list_of_possible_edges():
    # Sorting the list by distance

    return None



if __name__ == "__main__":
    # Example from Alcaide et. al. Similarity should be 0.56.
    patient1 = ['99662', '99591', '5990', '4019']
    patient2 = ['4329', '43491', '99702', '99591', '5990', '4019']

    print( round(distance(patient1, patient2), 2) )

    patient3 = {16: ['99662', '99591', '5990', '4019']}

    patient4 = {17: ['4329', '43491', '99702', '99591', '5990', '4019']}
    make_distance_matrix([patient3, patient4])
    read_csv('icd_diagnostic_categories.csv')