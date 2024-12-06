import subprocess
import os
import argparse

def clean_after_error(out_path):
    subprocess.run(["rm", "-r", out_path])

def run_program(command_line, out_path, label):
    result = subprocess.run(command_line, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout)
    if result.stderr or result.returncode != 0:
        print(f"error in {label}")
        print(result.stderr)
        clean_after_error(out_path)
        return 1
    return 0;

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pg", "--genomepaths", dest = "paths_to_genome", help="Path to input file holding paths to genome files", required=True)
    parser.add_argument("-pa", "--annotationpaths", dest = "paths_to_annotations", help="Path to input file holding paths to gtf files", required=True)
    parser.add_argument("-l", "--level", dest = "level_flag", help="Level of preprocessing Options: genome, CDS, gene, transcript, exon")
    parser.add_argument("-t", "--kmersize", dest = "kmer_size", default = "31", help="Size of kmers to be used in similarity comparison")
    parser.add_argument("-d", "--dashing", dest = "dashing_flag", action="store_true", help = "Enable dashing for processing")
    parser.add_argument("-k", "--bottomkparam", dest = "bottom_k_param", default = "200", help="Sketch size for minHash")
    parser.add_argument("-c", "--clean", dest = "clean", action="store_true", help="Clean intermediate files")
    parser.add_argument("-fo", "--fna_output", dest = "fna_out_path", help="Folder used to output preprocessing")
    parser.add_argument("-o", "--output", dest = "out_path", help="Folder used for tree", required=True)

    try:
        args = parser.parse_args()
    except:
        print("error: missing required arguments")
        return 1

    # make the output path
    subprocess.run(["mkdir", args.out_path])
    # make the input folder for the tree
    tree_input_path = os.path.join(args.out_path, "input")
    subprocess.run(["mkdir", tree_input_path])
    # make the fna output path
    fna_out_path = args.fna_out_path
    if args.fna_out_path is None:
        fna_out_path = os.path.join(args.out_path, "fna_temp")
    subprocess.run(["mkdir", fna_out_path])

    levels = ["genome", "CDS", "gene", "transcript", "exon"]
    if args.level_flag is not None:
        if args.level_flag not in levels:
            print("error: unknown level")
            clean_after_error()
            return 1       
        levels = [args.level_flag]

    for l in levels:
        if run_program(["python3", "extract_multiple.py", l, args.paths_to_genome, args.paths_to_annotations, fna_out_path], args.out_path, "preprocessing"):
            return 1
        print(l, "preprocessing complete")
        preprocessing_output = os.path.join(fna_out_path, f"output_{l}_paths.txt")

        # minhash or dashing
        if not args.dashing_flag:
            minhash_output = f"minhash_out_{l}.csv"
            minhash_output_path = os.path.join(tree_input_path, minhash_output)
            if run_program(["./minhash", preprocessing_output, args.kmer_size, args.bottom_k_param, minhash_output_path, "-d"], args.out_path, "minhash"):
                return 1
            print(minhash_output, "complete")
        else:
            print("dashing")
            return 0
        subprocess.run(["rm", preprocessing_output])
    # clean the intermediate file
    
    if (args.clean):
        subprocess.run(["rm", "-r", fna_out_path])
    
    # excluded for testing    

    # trees!
    consensus_location = os.path.join("tree_building", "consensus.py")
    if run_program(["python3", consensus_location, args.out_path], args.out_path, "tree building"):
        return 1
    if (args.clean):
        subprocess.run(["rm", "-r", tree_input_path])
    print("tree building complete")
    return 0





if __name__ == "__main__":
    main()