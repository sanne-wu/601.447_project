import pandas as pd
def read_distance_csv(filePath):
    """Reads a CSV file containing the distance matrix."""
    distanceDf = pd.read_csv(filePath, index_col=0, na_values=['-'])
    distanceDf = distanceDf.astype(float)
    return distanceDf

temp1 = read_distance_csv('test/input/minhash_out_CDS.csv')
temp2 = read_distance_csv('temp/input/distance_matrix_cds.csv')
print(temp1)
print(temp2)