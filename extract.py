import sys

def parse_gtf(gtf_file, feature):
  features = []
  with open(gtf_file, 'r') as f:
    for line in f:
      if line.startswith("#") or feature not in line:
        continue
      entry = line.strip().split("\t")
      if entry[2] == feature:
        chr = entry[0]
        start = int(entry[3])
        end = int(entry[4])
        strand = entry[6]
        features.append((chr, start, end, strand))
  return features

def load_genome(fna_file):
    genome = {} #genome as a dict
    chr_num = ""
    seq = []

    with open(fna_file, 'r') as f:
      for line in f:
        if line.startswith(">"):
          if chr_num:
              genome[chr_num] = "".join(seq)
          chr_num = line.split()[0][1:].split("_")[0] 
          seq = []
        else:
          seq.append(line.strip())
      if chr_num: 
        genome[chr_num] = "".join(seq)
    return genome

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in reversed(seq))

def extract_and_write_features(genome, features, output_file):
    with open(output_file, 'w') as f:
      f.write(f">{output_file}_concatenated_features\n")
        
      for chr, start, end, strand in features:
        if chr in genome:
          seq = genome[chr][start-1:end] #adjust for 1-based indexing
          if strand == "-":
            seq = reverse_complement(seq)
          f.write(seq)

def main(genome_fna, gtf_file, output_file, feature):
  # if len(sys.argv) != 4:
  #   print("Usage: python3 extract.py <genome.fna> <annotation.gtf> <output_file>")
  #   sys.exit(1)
  
  # genome_fna = sys.argv[1]
  # gtf_file = sys.argv[2]
  # output_file = sys.argv[3]
    
  features = parse_gtf(gtf_file, feature)  
  genome = load_genome(genome_fna)
  extract_and_write_features(genome, features, output_file)
  
  
# if __name__ == "__main__":
#   main()


genome_fna = '/Users/haojun.xu/Desktop/Computational Genomics/Final Project/ASM584v2(K12)/ncbi_dataset/data/GCA_000005845.2/GCA_000005845.2_ASM584v2_genomic.fna'
gtf_file = '/Users/haojun.xu/Desktop/Computational Genomics/Final Project/ASM584v2(K12)/ncbi_dataset/data/GCA_000005845.2/genomic.gtf'
output_file = "/Users/haojun.xu/Desktop/Computational Genomics/Final Project/concatenated_exons_k12_alt.txt"
main(genome_fna, gtf_file, output_file, feature="CDS")  