import pandas as pd
import numpy as np
from skbio import DistanceMatrix
from skbio.tree import nj
from dendropy import Tree, TaxonNamespace
from dendropy.calculate import treecompare
import glob
import matplotlib.pyplot as plt
from io import StringIO
from Bio import Phylo
import io
import dendropy

def readSimilarityCsv(filePath):
    # Reads a CSV file containing the similarity matrix.
    similarityDf = pd.read_csv(filePath, index_col=0)
    return similarityDf

def similarityToDistance(similarityDf):
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

def skbioToDendropyTree(skbioTree):
    # Use StringIO to capture the Newick output as a string
    newickIo = io.StringIO()
    skbioTree.write(newickIo, format='newick')
    newickStr = newickIo.getvalue()
    newickIo.close()
    
    # Convert Newick string to a DendroPy Tree
    dendropyTree = dendropy.Tree.get(data=newickStr, schema="newick")
    return dendropyTree

def calculateBootstrapSupport(originalTree, bootstrapTrees):

    if not bootstrapTrees:
        raise ValueError("The list of bootstrap trees is empty. Please provide bootstrap trees for support calculation.")
    
    for node in originalTree.postorderNodeIter():
        node.bootstrapSupport = 0

    originalTree.encodeBipartitions()
    splitCounts = {bipartition: 0 for bipartition in originalTree.bipartitionEdgeMap}

    for bt in bootstrapTrees:
        bt.encodeBipartitions()
        for bipartition in bt.bipartitionEdgeMap:
            if bipartition in splitCounts:
                splitCounts[bipartition] += 1

    numBootstrapTrees = len(bootstrapTrees)
    for edge in originalTree.postorderEdgeIter():
        if edge.headNode is not None and not edge.headNode.isLeaf():
            bipartition = edge.bipartition
            support = (splitCounts.get(bipartition, 0) / numBootstrapTrees) * 100
            edge.headNode.label = f"{support:.1f}%"

    return originalTree

def visualizeTree(dendropyTree, outputImagePath):
    newickStr = dendropyTree.asString(schema='newick')

    from Bio import Phylo
    handle = StringIO(newickStr)
    phyloTree = Phylo.read(handle, 'newick')

    # Draw the tree using matplotlib
    fig = plt.figure(figsize=(12, 8))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(phyloTree, do_show=False, axes=axes)

    # Save the figure to a file
    plt.savefig(outputImagePath, format='png', dpi=300)
    plt.close(fig)
    print(f"Annotated phylogenetic tree visualization saved to {outputImagePath}.")

def main(originalSimilarityCsv, bootstrapSimilarityDir, outputTreePath, outputImagePath):
    similarityDf = readSimilarityCsv(originalSimilarityCsv)
    print("Original similarity matrix loaded.")

    distanceDf = similarityToDistance(similarityDf)
    print("Original similarity matrix converted to distance matrix.")

    originalSkbioTree = constructTree(distanceDf)
    print("Original phylogenetic tree constructed using Neighbor-Joining method.")

    originalTree = skbioToDendropyTree(originalSkbioTree)

    bootstrapTrees = []
    bootstrapFiles = glob.glob(f"{bootstrapSimilarityDir}/*.csv")
    print(f"Found {len(bootstrapFiles)} bootstrap similarity matrices.")

    for i, filePath in enumerate(bootstrapFiles):
        simDf = readSimilarityCsv(filePath)
        distDf = similarityToDistance(simDf)
        skbioTree = constructTree(distDf)
        dendropyTree = skbioToDendropyTree(skbioTree)
        bootstrapTrees.append(dendropyTree)
        print(f"Processed bootstrap tree {i+1}/{len(bootstrapFiles)}")

    annotatedTree = calculateBootstrapSupport(originalTree, bootstrapTrees)
    print("Bootstrap support values calculated and annotated on the original tree.")
    annotatedTree.write(path=outputTreePath, schema='newick')
    print(f"Annotated phylogenetic tree saved to {outputTreePath}.")
    visualizeTree(annotatedTree, outputImagePath)

if __name__ == "__main__":
    originalSimilarityCsv = "/users/harry/desktop/Computational Gentics/Final project/similarity_matrix.csv"  # Path to the original similarity matrix
    bootstrapSimilarityDir = "/users/harry/desktop/Computational Gentics/Final project/bootstrap_matrices"    # Directory containing bootstrap similarity matrices
    outputTreePath = "/users/harry/desktop/Computational Gentics/Final project/annotated_phylogenetic_tree.nwk"  # Output Newick file path
    outputImagePath = "/users/harry/desktop/Computational Gentics/Final project/annotated_phylogenetic_tree.png" # Output image file path
    main(originalSimilarityCsv, bootstrapSimilarityDir, outputTreePath, outputImagePath)
