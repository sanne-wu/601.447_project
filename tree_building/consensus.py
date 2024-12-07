import argparse
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from Bio import Phylo
from Bio.Phylo.Newick import Tree, Clade
from collections import Counter
import glob
from NN import NN
from tree_building.UPGMA import UPGMA
import logging
from HRA import HRA
from HRC import HRC

def setupLogging(logFilePath):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(logFilePath),
            logging.StreamHandler()
        ]
    )

def saveToPhylip(distances, labels, filePath):
    n = len(labels)
    with open(filePath, 'w') as f:
        f.write(f"{n}\n")
        for i, label in enumerate(labels):
            formattedLabel = label[:10].ljust(10)
            f.write(f"{formattedLabel} ")
            f.write(" ".join(f"{distances[i, j]:.5f}" for j in range(n)))
            f.write("\n")

def loadDistanceMatrix(filePath):
    df = pd.read_csv(filePath, index_col=0)
    df.replace("-", np.nan, inplace=True)
    distanceMatrix = df.to_numpy(dtype=float)
    iUpper = np.triu_indices_from(distanceMatrix, 1)
    distanceMatrix[(iUpper[1], iUpper[0])] = distanceMatrix[iUpper]
    return distanceMatrix, df.columns.tolist()

def linkageToPhylo(linkageMatrix, labels):
    treeRoot, _ = to_tree(linkageMatrix, rd=True)
    def buildClade(node):
        if node.is_leaf():
            return Clade(name=labels[node.id])
        else:
            clade = Clade()
            clade.clades.append(buildClade(node.left))
            clade.clades.append(buildClade(node.right))
            clade.branch_length = node.dist
            return clade
    phyloTree = Tree(root=buildClade(treeRoot))
    return phyloTree

def upgmaTree(distanceMatrix, labels, title, ax):
    condensedDist = distanceMatrix[np.triu_indices_from(distanceMatrix, 1)]
    linkageMatrix = linkage(condensedDist, method='average')
    dendrogram(linkageMatrix, labels=labels, orientation='left', leaf_rotation=0, leaf_font_size=10, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Distance")
    ax.set_ylabel("Clusters")
    tree = linkageToPhylo(linkageMatrix, labels)
    return tree

def getSplits(tree):
    splits = []
    clades = list(tree.find_clades(order='level'))
    for clade in clades:
        if clade.is_terminal():
            continue
        taxaNames = frozenset(taxon.name for taxon in clade.get_terminals())
        splits.append(taxaNames)
    return splits

def buildConsensusTree(allTaxa, majoritySplits):
    consensusTree = Tree()
    consensusTree.root.clades = [Clade(name=taxon) for taxon in sorted(allTaxa)]
    changed = True
    while changed:
        changed = False
        for split in list(majoritySplits):
            for clade in consensusTree.get_nonterminals(order='postorder'):
                cladeTaxa = {taxon.name for taxon in clade.get_terminals()}
                if split < cladeTaxa and split != cladeTaxa:
                    taxaInSplit = []
                    taxaNotInSplit = []
                    for subclade in clade.clades:
                        subcladeTaxa = {taxon.name for taxon in subclade.get_terminals()}
                        if subcladeTaxa & split:
                            taxaInSplit.append(subclade)
                        else:
                            taxaNotInSplit.append(subclade)
                    if taxaInSplit and taxaNotInSplit:
                        newClade = Clade()
                        newClade.clades = taxaInSplit
                        clade.clades = [newClade] + taxaNotInSplit
                        changed = True
                        majoritySplits.remove(split)
                        break
    return consensusTree

def buildConsensus(generalFolder):
    consensusFolder = os.path.join(generalFolder, "result", "Consensus")
    os.makedirs(consensusFolder, exist_ok=True)
    nnFolder = os.path.join(generalFolder, "result", "NN")
    upgmaFolder = os.path.join(generalFolder, "result", "UMPGA")
    hrcFolder = os.path.join(generalFolder, "result", "HRC")
    hraFolder = os.path.join(generalFolder, "result", "HRA")
    treeFiles = glob.glob(os.path.join(nnFolder, "*.nex")) + glob.glob(os.path.join(upgmaFolder, "*.nex")) + glob.glob(os.path.join(hrcFolder, "*.nex")) + glob.glob(os.path.join(hraFolder, "*.nex"))
    categories = ['cds', 'whole_genome','gene']
    trees_by_category = {category: [] for category in categories}
    for filePath in treeFiles:
        if not os.path.exists(filePath):
            logging.warning(f"File not found: {filePath}")
            continue
        treesInFile = list(Phylo.parse(filePath, 'nexus'))
        if not treesInFile:
            logging.warning(f"No trees found in file: {filePath}")
            continue
        for idx, tree in enumerate(treesInFile):
            if idx < 3:
                category = categories[idx]
                trees_by_category[category].append(tree)
    for category in categories:
        logging.info(f"Processing consensus for category: {category}")
        splitCounts = Counter()
        allTaxa = set()
        for tree in trees_by_category[category]:
            allTaxa.update(taxon.name for taxon in tree.get_terminals())
            splits = getSplits(tree)
            splitCounts.update(splits)
        majorityThreshold = len(trees_by_category[category]) / 2
        majoritySplits = {split for split, count in splitCounts.items() if count > majorityThreshold}
        logging.info(f"Number of majority splits for {category}: {len(majoritySplits)}")
        if majoritySplits:
            consensusTree = buildConsensusTree(allTaxa, majoritySplits)
            outputFileNex = os.path.join(consensusFolder, f"Consensus_Majority_Rule_Tree_{category}.nex")
            outputFileNwk = os.path.join(consensusFolder, f"Consensus_Majority_Rule_Tree_{category}.nwk")
            Phylo.write(consensusTree, outputFileNex, 'nexus')
            Phylo.write(consensusTree, outputFileNwk, 'newick')
            plt.figure(figsize=(10, 8))
            Phylo.draw(consensusTree, do_show=False)
            plt.title(f'Majority-Rule Consensus Tree - {category}')
            consensusPlotPath = os.path.join(consensusFolder, f"Majority_Rule_Consensus_Tree_{category}.png")
            plt.savefig(consensusPlotPath, dpi=300)
            plt.close()
            logging.info(f"Consensus tree saved to {outputFileNex}")
            logging.info(f"Consensus tree saved to {outputFileNwk}")
            logging.info(f"Consensus tree plot saved at: {consensusPlotPath}")
        else:
            logging.warning(f"No majority splits found to build a consensus tree for {category}.")

def main():
    parser = argparse.ArgumentParser(description="Process distance matrices to generate UPGMA trees and consensus trees.")
    parser.add_argument('generalFolder', type=str, help='Path to the general folder containing input and result directories.')
    parser.add_argument('--log', type=str, default='phylo_pipeline.log', help='Path to the log file.')
    args = parser.parse_args()
    generalFolder = args.generalFolder
    logFilePath = os.path.join(generalFolder, args.log)
    setupLogging(logFilePath)
    if not os.path.isdir(generalFolder):
        logging.error(f"The specified general folder does not exist: {generalFolder}")
        return
    logging.info(f"Starting phylogenetic analysis in folder: {generalFolder}")
    NN(generalFolder)
    UMPGA(generalFolder)
    minEvolution(generalFolder)
    FM(generalFolder)
    HRC(generalFolder)
    HRA(generalFolder)
    buildConsensus(generalFolder)
    logging.info("NN, UPGMA clustering, and consensus tree generation completed successfully.")

if __name__ == "__main__":
    main()
