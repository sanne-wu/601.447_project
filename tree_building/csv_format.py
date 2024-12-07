import pandas as pd
import numpy as np
import sys

# Load the CSV file

def main(fn_in, fn_out):
    # Read TSV file with tab separator explicitly specified
    df = pd.read_csv(fn_in, index_col=0, sep='\t')
    
    # Clean up index and column names by taking just the filename part
    df.index = df.index.map(lambda x: x.split('\t')[0])
    df.columns = df.columns.map(lambda x: x.split('\t')[0])
    
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
    sorted_df.to_csv(fn_out)

    print(f"The sorted and completed matrix has been saved to {fn_out}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
