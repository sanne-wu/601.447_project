import subprocess
import os
import argparse

def clean_after_error(out_path):
    subprocess.run(["rm", "-r", out_path])

def run_program(command_line, out_path, label, ex = False):
    result = subprocess.run(command_line, capture_output=True, text=True)
    if not ex:
        if result.stdout:
            print(result.stdout)
        if result.stderr or result.returncode != 0:
            print(f"error in {label}")
            print(result.stderr)
            clean_after_error(out_path)
            return 1
    return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pg", "--genomepaths", dest = "paths_to_genome", help="Path to input file holding paths to genome files", required=True)
    parser.add_argument("-pa", "--annotationpaths", dest = "paths_to_annotations", help="Path to input file holding paths to gtf files", required=True)
    parser.add_argument("-l", "--level", dest = "level_flag", help="Level of preprocessing Options: genome, CDS, gene")
    parser.add_argument("-t", "--kmersize", dest = "kmer_size", default = "31", help="Size of kmers to be used in similarity comparison")
    parser.add_argument("-k", "--bottomkparam", dest = "bottom_k_param", default = "200", help="Sketch size for minHash")
    parser.add_argument("-c", "--clean", dest = "clean", action="store_true", help="Clean intermediate files")
    parser.add_argument("-fo", "--fna_output", dest = "fna_out_path", help="Folder used to output preprocessing")
    parser.add_argument("-o", "--output", dest = "out_path", help="Folder used for tree", required=True)
    parser.add_argument("-d", "--dashing", dest = "dashing_flag", action="store_true", help = "Enable dashing for processing")

    try:
        args = parser.parse_args()
    except:
        print("error: missing required arguments")
        return 1

    # Create output directories
    os.makedirs(args.out_path, exist_ok=True)
    tree_input_path = os.path.join(args.out_path, "input")
    os.makedirs(tree_input_path, exist_ok=True)
    
    fna_out_path = args.fna_out_path or os.path.join(args.out_path, "fna_temp")
    os.makedirs(fna_out_path, exist_ok=True)

    levels = ["genome", "CDS", "gene"]
    if args.level_flag is not None:
        if args.level_flag not in levels:
            print("Error: Unknown preprocessing level")
            clean_after_error(args.out_path)
            return 1       
        levels = [args.level_flag]

    for l in levels:
        if run_program(["python3", "extract_multiple.py", l, args.paths_to_genome, args.paths_to_annotations, fna_out_path], args.out_path, "preprocessing"):
            return 1
        print(l, "preprocessing complete")
        preprocessing_output = os.path.join(fna_out_path, f"output_{l}_paths.txt")

        # Run similarity analysis using either MinHash or Dashing
        if not args.dashing_flag:
            minhash_output = f"minhash_out_{l}.csv"
            minhash_output_path = os.path.join(tree_input_path, minhash_output)
            if run_program(
                ["./minhash", preprocessing_output, args.kmer_size, 
                 args.bottom_k_param, minhash_output_path, "-d"],
                args.out_path,
                "minhash"
            ):
                return 1
            print(f"{minhash_output} complete.")
        else:
            # make sure the alias was added to ~/.bashrc according to the README installation guideline
            # read more in the "Dependencies" section
            dashing_output = f'dashing_out_{l}.tsv'
            dashing_output_path = os.path.join(tree_input_path, dashing_output)
            level_to_suffix = {
                "genome": "full_genome", 
                "CDS": "cds", 
                "gene": "gene"
            }
            dashing_cmd = [
                './dashing', 'dist',
                '-M', # compute Mash distance
                '-k', '15', # kmer size for E.coli, see Bonnie et al. 2024 for more details
                '-p', '40', # threads, assume running on salz servers
                '-O', dashing_output_path,
                f'*_{level_to_suffix[l]}.fna'
            ]
            if run_program(dashing_cmd, args.out_path, "dashing"):
                return 1
            print(f"{dashing_output} complete")

        os.remove(preprocessing_output)

    # Clean intermediate files if requested
    if args.clean:
        subprocess.run(["rm", "-r", fna_out_path])
    
    # Build phylogenetic trees
    consensus_location = os.path.join("tree_building", "consensus.py")
    if run_program(["python3", consensus_location, args.out_path], args.out_path, "tree building", True):
        return 1
        
    if args.clean:
        subprocess.run(["rm", "-r", tree_input_path])
    print("Tree building complete")
    return 0

if __name__ == "__main__":
    main()