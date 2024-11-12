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

def load_genome(fna_file, prokaryotes):
    genome = {} #genome as a dict
    chr_num = ""
    seq = []

    with open(fna_file, 'r') as f:
      for line in f:
        if line.startswith(">"):
          if chr_num:
            genome[chr_num] = "".join(seq)
          chr_num = line.split()[0][1:].split(" ")[0] if prokaryotes else line.split()[0][1:].split("_")[0]
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
      for chr, start, end, strand in features:
        if chr in genome:
          seq = genome[chr][start-1:end] #adjust for 1-based indexing
          if strand == "-":
            seq = reverse_complement(seq)
          f.write(seq)

def main():
  if len(sys.argv) != 5:
    print("Usage: python3 extract.py feature <genome.fna> <annotation.gtf> <output_file>")
    sys.exit(1)
    
  feature = sys.argv[1] #feature to extract, ex.CDS
  genome_fna = sys.argv[2]
  gtf_file = sys.argv[3]
  output_file = sys.argv[4]
  
  features = parse_gtf(gtf_file, feature)  
  genome = load_genome(genome_fna, True)
  print(genome.keys())
  extract_and_write_features(genome, features, output_file)
  
  
if __name__ == "__main__":
  main()
