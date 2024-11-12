import pandas as pd
import numpy as np
from skbio import TreeNode, DistanceMatrix
from skbio.tree import nj

def readSimilarityCsv(filePath):
    # Reads a CSV file containing the similarity matrix.
    similarityDf = pd.read_csv(filePath, index_col=0)
    return similarityDf

def similarityToDistance(similarityDf):
    # Ensure all similarities are between 0 and 1
    similarityDf = similarityDf.clip(lower=0, upper=1)
    epsilon = 1e-10
    similarityDf += epsilon
    distanceDf = -np.log(similarityDf)
    np.fill_diagonal(distanceDf.values, 0)
    return distanceDf

def constructTree(distanceDf):
    ids = distanceDf.index.tolist()
    distanceArray = distanceDf.values
    dm = DistanceMatrix(distanceArray, ids)
    tree = nj(dm)
    return tree

def main(similarityCsvPath, outputTreePath):
    distanceDf = readSimilarityCsv(similarityCsvPath)
    # print("Similarity matrix loaded.")
    # distanceDf = similarityToDistance(similarityDf)
    print("Similarity matrix converted to distance matrix.")
    tree = constructTree(distanceDf)
    print("Phylogenetic tree constructed using Neighbor-Joining method.")
    tree.write(outputTreePath)
    print(f"Phylogenetic tree saved to {outputTreePath}.")
    print("\nConstructed Phylogenetic Tree:")
    print(tree.ascii_art())

if __name__ == "__main__":
    similarityCsvPath = "/users/harry/desktop/Computational Genetics/Final project/code/dashing_geno_distance_matrix.csv"  # Replace with your CSV file path
    outputTreePath = "/users/harry/desktop/Computational Genetics/Final project/phylogenetic_tree.nwk"  # Output Newick file path
    main(similarityCsvPath, outputTreePath)
