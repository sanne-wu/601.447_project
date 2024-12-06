# Sketch2Phylo
# needed packages:
scipy, matplotlib, numpy, pandas, biopython, dendropy

# example usage:
`python3 wrapper.py -pg preprocessing_test/input_genome_paths.txt -pa preprocessing_test/input_gtf_paths.txt -o test -c`
### parameters
**`-pg, --genomepaths`:** (required) Path to input file holding paths separated by line breaks to genome files of interest. \
**`-pa, --annotationpaths`:** (required) Path to input file holding paths separated by line breaks to gene annotation files of interest. Must be in the same order as the genome paths. \
**`-l, --level`:** (optional) The level of specifity to preprocess the genome files at. Possible options are `genome`, `CDS`, `gene`, `transcript`, `exon`. If this is not specified, it will produce results for all five levels. \
**`-t, --kmersize`:** (optional) The size of kmers to be used in the similarity comparison. Default is 31.\
**`-k, --bottomkparam`:** (optional) Used only for minHash. Specifies the sketch size when computing the bottom-k sketch. Default is 200. \
**`-c, --clean`:** (optional) Specifying the `-c` tag will clean intermediate files produced by the program. This includes the distance matrices produced by Dashing and minHash, along with the preprocessed .fna files.\
**`-fo, --fna_output`:** (optional) Specifies the path to output the processed .fna in the preprocessing stage. \
**`-o, --output`:** (required) Specifies the folder that will contain the output of the tree.
**`-d, --dashing`:** (optional) Specifying the `-d` tag will produce the distance matrix with dashing instead of minHash.

Beyond running the program as an entire pipeline, each part can also be ran separately.
### 1. preprocessing:
#### usage:
`python3 extract_multiple.py <feature> <genome_list> <gtf_list> <output_folder>`
#### example usage:
`python3 extract_multiple.py exon preprocessing_test/input_genome_paths.txt preprocessing_test/input_gtf_paths.txt preprocessing_test/output_exons`

The processed files will be stored in the path `output_folder`. \
The formatted file used by minHash is stored in `output_feature_paths` where `feature` is one of `genome`, `CDS`, `gene`, `transcript`, `exon` and corresponds to the `-l` parameter in the wrapper.

### 2. minHash:
#### usage:
`./minhash <list_name> <kmer size> <bottom-k parameter> <output_name (optional)> <-d (optional)>` 
#### example usage:
This reproduces the Table I in the write-up.\
`./minhash list.txt 31 200`

Adding the `-d` flag will compute a distance matrix variation by calculating -log(.) each of the entries from the estimated Jaccard similarity matrix. \
If an output is not provided, it will print the results directly to console. \
The `list_name` parameter is a formatted file where each contain contains the .fna file and its name separated by a white space. 
`kmer_size` and `bottom-k parameter` correspond to the `-t` and `-k` parameters in the wrapper above. 
