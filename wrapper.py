import subprocess
import os
import argparse

def clean_after_error(out_path):
    subprocess.run(["rm", "-r", out_path])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pg", "--genomepaths", dest = "paths_to_genome", help="Path to input file holding paths to genome files", required=True)
    parser.add_argument("-pa", "--annotationpaths", dest = "paths_to_annotations", help="Path to input file holding paths to gtf files", required=True)
    parser.add_argument("-l", "--level", dest = "level_flag", help="Level of preprocessing Options: genome, CDS, gene, transcript, exon")
    parser.add_argument("-t", "--kmersize", dest = "kmer_size", default = "31", help="Size of kmers to be used in similarity comparison")
    parser.add_argument("-d", "--dashing", dest = "dashing_flag", action="store_true", help = "Enable dashing for processing")
    parser.add_argument("-k", "--bottomkparam", dest = "bottom_k_param", default = "200", help="Sketch size for minHash")
    parser.add_argument("-c", "--clean", dest = "clean", action="store_true", help="Clean intermediate files")
    parser.add_argument("-o", "--output", dest = "out_path", help="Folder used for tree", required=True)
    try:
        args = parser.parse_args()
    except:
        print("error: missing required arguments")
        return 1
    # # inputs to preprocessing
    print(args.paths_to_genome)
    print(args.paths_to_annotations)
    print(args.level_flag)
    print(args.kmer_size)
    print(args.dashing_flag)
    print(args.bottom_k_param)
    print(args.clean)
    print(args.out_path)


    # make the output path
    subprocess.run(["mkdir", args.out_path])
    # make the input folder for the tree
    tree_input_path = os.path.join(args.out_path, "input")
    subprocess.run(["mkdir", tree_input_path])

    levels = ["genome", "CDS", "gene", "transcript", "exon"]
    if args.level_flag is not None:
        if args.level_flag not in levels:
            print("error: unknown level")
            clean_after_error()
            return 1       
        levels = [args.level_flag]
    for l in levels:
        preprocessing_output = "preprocessed_out.txt"
        result = subprocess.run(["python3", "extract_multiple.py", l, args.paths_to_genome, args.paths_to_annotations, preprocessing_output], capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print("error in preprocessing")
            print(result.stderr)
            return 1
        print(l, "preprocessing complete")

        # minhash or dashing
        if not args.dashing_flag:
            minhash_output = f"minhash_out_{l}.csv"
            minhash_output_path = os.path.join(tree_input_path, minhash_output)
            result = subprocess.run(["./minhash", preprocessing_output, args.kmer_size, args.bottom_k_param, minhash_output_path, "-d"], capture_output=True, text=True)
            if result.stdout:
                print(result.stdout)
            if result.returncode != 0:
                print("error in minhash")
                print(result.stderr)
                clean_after_error()
                return 1
            print(minhash_output, "complete")
        else:
            print("dashing")
            return 0
    # clean the intermediate file
    subprocess.run(["rm", preprocessing_output])

    # trees!
    # consensus_location = os.path.join("tree_building", "consensus.py")
    # result = subprocess.run(["python3", consensus_location, args.out_path], capture_output=True, text=True)

    # # # inputs to minhash via command line
    # subprocess.run(["rm", preprocessing_output])

    
    # # subprocess.run(["rm", "minhash_out.csv"])

    return 0





if __name__ == "__main__":
    main()