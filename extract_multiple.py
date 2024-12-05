import sys
import os

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
  genome = {}
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
        seq = genome[chr][start-1:end]
        if strand == "-":
          seq = reverse_complement(seq)
        f.write(f">{chr}:{start}-{end} ({strand})\n{seq}\n")

def main():
  if len(sys.argv) != 5:
    print("Usage: python3 extract_line_by_line_with_names.py feature <genome_list.txt> <gtf_list.txt> <output_paths.txt>")
    sys.exit(1)

  feature = sys.argv[1]
  genome_list_file = sys.argv[2]
  gtf_list_file = sys.argv[3]
  output_paths_file = sys.argv[4]

  with open(genome_list_file, 'r') as genome_list, open(gtf_list_file, 'r') as gtf_list, open(output_paths_file, 'w') as out_paths:
    for genome_line, gtf_line in zip(genome_list, gtf_list):
      genome_line = genome_line.strip()
      gtf_line = gtf_line.strip()

      if not genome_line or not gtf_line:
        print(f"Skipping empty line: genome='{genome_line}', gtf='{gtf_line}'")
        continue

      # Parse genome_name and paths
      genome_name, genome_path = genome_line.split(" ", 1)
      gtf_name, gtf_path = gtf_line.split(" ", 1)

      if genome_name != gtf_name:
        print(f"Mismatch in genome and GTF names: {genome_name} vs {gtf_name}")
        continue

      output_file = f"{genome_name}_{feature}_sequences.fna"

      #process the genomes and GTF file
      features = parse_gtf(gtf_path, feature)
      genome = load_genome(genome_path, True)
      extract_and_write_features(genome, features, output_file)

      #write output path to <output_paths.txt>
      out_paths.write(f"{os.path.abspath(output_file)}\n")

if __name__ == "__main__":
    main()