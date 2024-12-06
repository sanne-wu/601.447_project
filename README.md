# Sketch2Phylo
### needed packages:
scipy, matplotlib, numpy, pandas, biopython, dendropy

### example usage:
`python3 wrapper.py -pg preprocessing_test/input_genome_paths.txt -pa preprocessing_test/input_gtf_paths.txt -o test -c`
### parameters
**-pg, --genomepaths:** (required) Path to input file holding paths separated by line breaks to genome files of interest. \
**-pa, --annotationpaths:** (required) Path to input file holding paths separated by line breaks to gene annotation files of interest. Must be in the same order as the genome paths. \
**-l, --level:** (optional) The level of specifity to preprocess the genome files at. Possible options are `genome`, `CDS`, `gene`, `transcript`, `exon`. If this is not specified, it will produce results for all five levels. \
**-t, --kmersize:** (optional) The size of kmers to be used in the similarity comparison. Default is 31.\
**-k, --bottomkparam:** (optional) Used only for minHash. Specifies the sketch size when computing the bottom-k sketch. Default is 200. \
**-c, --clean:** (optional) Specifying the `-c` tag will clean intermediate files produced by the program. This includes the distance matrices produced by Dashing and minHash, along with the preprocessed .fna files.\
**-fo, --fna_output:** (optional) Specifies the path to output the processed .fna in the preprocessing stage. \
**-o, --output:** (required) Specifies the folder that will contain the output of the tree.
**-d, --dashing:** (optional) Specifying the `-d` tag will produce the distance matrix with dashing instead of minHash.