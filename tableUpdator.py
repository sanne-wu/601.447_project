import pandas as pd
import numpy as np

# Load the CSV file
file_path = '/users/harry/desktop/Computational Genetics/Final project/code/distance_matrix_cds.csv'  # Replace with the path to your CSV file
df = pd.read_csv(file_path, index_col=0)

# Convert to a numpy array for manipulation
matrix = df.to_numpy()

# Fill the lower triangle with the upper triangle values
for i in range(matrix.shape[0]):
    for j in range(i + 1, matrix.shape[1]):
        matrix[j, i] = matrix[i, j]

# Set the diagonal to 0
np.fill_diagonal(matrix, 0)

# Convert back to DataFrame
completed_df = pd.DataFrame(matrix, index=df.index, columns=df.columns)

# Sort the matrix alphabetically by row and column names
sorted_indices = sorted(completed_df.index)
sorted_df = completed_df.loc[sorted_indices, sorted_indices]

# Save the sorted and completed matrix to a new CSV file
sorted_df.to_csv('distance_matrix_cds.csv')

print("The sorted and completed matrix has been saved to 'sorted_completed_distance_matrix.csv'")
