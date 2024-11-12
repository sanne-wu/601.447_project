import pandas as pd
import numpy as np
from skbio import DistanceMatrix
from skbio.tree import nj

def readSimilarityCsvs(filePaths):
    similarityDfs = []
    for filePath in filePaths:
        similarityDf = pd.read_csv(filePath, index_col=0)
        similarityDfs.append(similarityDf)
    return similarityDfs

def similarityToDistance(similarityDfs):
    distanceDfs = []
    for similarityDf in similarityDfs:
        # Ensure all similarities are between 0 and 1
        similarityDf = similarityDf.clip(lower=0, upper=1)
        epsilon = 1e-10
        similarityDf += epsilon
        distanceDf = -np.log(similarityDf)

        # Set the diagonal to zero for a valid distance matrix
        np.fill_diagonal(distanceDf.values, 0)
        distanceDfs.append(distanceDf)
    return distanceDfs

def computeConsensusDistance(distanceDfs):
    # Stack the distance matrices along a new axis
    stackedDistances = np.stack([df.values for df in distanceDfs], axis=2)
    consensusDistances = np.mean(stackedDistances, axis=2)
    ids = distanceDfs[0].index.tolist()
    consensusDf = pd.DataFrame(consensusDistances, index=ids, columns=ids)
    return consensusDf

def constructTree(distanceDf):
    ids = distanceDf.index.tolist()
    distanceArray = distanceDf.values
    dm = DistanceMatrix(distanceArray, ids)
    tree = nj(dm)
    return tree

def main(similarityCsvPaths, outputTreePath):

    distanceDfs = readSimilarityCsvs(similarityCsvPaths)
    # print(f"{len(similarityDfs)} similarity matrices loaded.")

    # distanceDfs = similarityToDistance(similarityDfs)
    # print("Similarity matrices converted to distance matrices.")

    consensusDistanceDf = computeConsensusDistance(distanceDfs)
    print("Consensus distance matrix computed by averaging.")

    tree = constructTree(consensusDistanceDf)
    print("Phylogenetic tree constructed using Neighbor-Joining method.")

    with open(outputTreePath, 'w') as f:
        f.write(tree.write("/users/harry/desktop/Computational Genetics/Final project/phylogenetic_tree.nwk") + ";")
    print(f"Phylogenetic tree saved to {outputTreePath}.")

    print("\nConstructed Phylogenetic Tree:")
    print(tree.ascii_art())

if __name__ == "__main__":
    similarityCsvPaths = [
        "/users/harry/desktop/Computational Genetics/Final project/code/dashing_geno_distance_matrix.csv",
        "/users/harry/desktop/Computational Genetics/Final project/code/dashing_distance_matrix_cds.csv",
    ]
    outputTreePath = "/users/harry/desktop/Computational Genetics/Final project/phylogenetic_tree.nwk"
    main(similarityCsvPaths, outputTreePath)