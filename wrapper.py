import subprocess
import sys

def main():
    if len(sys.argv) < 5 or len(sys.argv) > 6:
        print("Usage: python3 wrapper.py <paths_to_genome> <paths_to_annotations> <level_flag = cds, gene> <kmer size> <bottom-k parameter (optional)>")
        sys.exit(1)
    # inputs to wrapper and each component
    paths_to_genome = sys.argv[1]
    paths_to_annotations = sys.argv[2]
    level_flag = sys.argv[3]
    kmer_size = sys.argv[4]
    bottom_k_param = str(200)
    if (len(sys.argv) == 6):
        bottom_k_param = sys.argv[5]
    
    # inputs to preprocessing
    preprocessing_output = "list.txt"



    # call preprocessing with subprocess command line
    ###


    # inputs to minhash via command line
    minhash_output = "minhash_out.csv"
    subprocess.run(["./minhash", preprocessing_output, kmer_size, bottom_k_param, minhash_output], capture_output=True, text=True)
    

    # remove intermediate file
    # subprocess.run(["rm", "minhash_out.csv"])

if __name__ == "__main__":
    main()