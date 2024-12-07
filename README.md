# Sketch2Phylo
Sketch2Phylo uses minHash and Dashing to produce a phylogenetic tree from a set of genomes and their corresponding annotations.

## Dependencies:
The following Python packages are required:
* scipy
* matplotlib
* numpy 
* pandas
* biopython
* dendropy
* scikit-bio
* dashing (required for -d flag)
  ```bash
  git clone --recursive https://github.com/dnbaker/dashing
  cd dashing && make dashing
  alias dashing=$(pwd)/dashing
  ```

## Parameters
**`-pg, --genomepaths`:** (required) Path to input file holding paths separated by line breaks to genome files of interest. \
**`-pa, --annotationpaths`:** (required) Path to input file holding paths separated by line breaks to gene annotation files of interest. Must be in the same order as the genome paths. \
**`-l, --level`:** (optional) The level of specifity to preprocess the genome files at. Possible options are `genome`, `CDS`, `gene`. If this is not specified, it will produce results for all five levels. \
**`-t, --kmersize`:** (optional) The size of kmers to be used in the similarity comparison. Default is 31.\
**`-k, --bottomkparam`:** (optional) Used only for minHash. Specifies the sketch size when computing the bottom-k sketch. Default is 200. \
**`-c, --clean`:** (optional) Specifying the `-c` flag will clean intermediate files produced by the program. This includes the distance matrices produced by Dashing and minHash, along with the preprocessed .fna files.\
**`-fo, --fna_output`:** (optional) Specifies the path to output the processed .fna in the preprocessing stage. \
**`-o, --output`:** (required) Specifies the folder that will contain the output of the tree. \
**`-d, --dashing`:** (optional) Specifying the `-d` flag will produce the distance matrix with dashing instead of minHash.

## Running Sketch2Phylo
### Simple use case (minHash)
```bash
python3 wrapper.py -pg preprocessing_test/input_genome_paths.txt 
                   -pa preprocessing_test/input_gtf_paths.txt 
                   -o test -c
```
### Simple use case (Dashing)
```bash
python3 wrapper.py -pg preprocessing_test/input_genome_paths.txt
                   -pa preprocessing_test/input_gtf_paths.txt 
                   -o test_dashing -c -t 14 -d
```
In addition to running the program as an end-to-end pipeline, each part can also be ran separately for more control over the process.
### 1. preprocessing
#### Usage:
```bash
python3 extract_multiple.py <feature> <genome_list> <gtf_list> <output_folder>
```
#### Example use case:
```bash
python3 extract_multiple.py exon preprocessing_test/input_genome_paths.txt 
        preprocessing_test/input_gtf_paths.txt preprocessing_test/output_exons
```

The processed files will be stored in the path `output_folder`. \
The formatted file used by minHash is stored in `output_feature_paths` where `feature` is one of `genome`, `CDS`, `gene` and corresponds to the `-l` parameter in the wrapper.

### 2. minHash
#### Usage:
```bash
./minhash <list_name> <kmer size> <bottom-k parameter> <output_name (optional)> <-d (optional)>
```
#### Example use case:
This reproduces the Table II in the write-up.

```bash
./minhash list.txt 31 200
```

Adding the `-d` flag will compute a distance matrix variation by calculating -log(.) each of the entries from the estimated Jaccard similarity matrix. \
If an output is not provided, it will print the results directly to console. \
The `list_name` parameter is a formatted file where each contain contains the .fna file and its name separated by a white space. 
`kmer_size` and `bottom-k parameter` correspond to the `-t` and `-k` parameters in the wrapper above. 

### 3. consensus
#### Usage:
```bash
python3 consensus.py <tree_folder>
```
#### Example use case:
```bash
python3 consensus.py test
```

The program produces trees from different methods such as NN, UMPGA, and MinE, then forms a Consensus between them.\
Input is placed within `tree_folder/input/`, which is 3 similarity matrices at the different genome levels produced either by dashing or minHash. \
Results from each of these are output into the `tree_folder/result/` folder under their respective names.

## Reproducing results
To reproduce the consensus tree results, please run:
```bash
python3 wrapper.py -pg preprocessing_test/input_genome_paths.txt pa preprocessing_test/input_gtf_paths.txt -o test -c
```
This will produce the Consensus trees in `test/result/Consensus` in .nwk format, along with .png files of the trees, which is triplet of trees shown in the write-up. \
At the same time, the tree built by nearest neigbor, UPGMA(Unweighted Pair Group Method with Arithmetic Mean), average linkage heirachical cluster, and complete linkage heirachical cluster will be saved to result as well.
The reference tree is stored in `reference_tree.nwk`. \
In order to get the RF distance between the Consensus trees and the reference, we utilized Visual TreeCmp (https://eti.pg.edu.pl/TreeCmp/WEB). \
Selecting the "Ref-to-all comparison" we pasted the contents of the .nwk files for the consensus trees for cds, gene, and genome against the contents of `reference_tree.nwk` for comparison. This produces the RF distance table in the write-up.