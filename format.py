import sys
paths = sys.argv[1] 

with open(paths, 'r') as file_list:
  for line in file_list:
      path, identity = line.strip().split()
      with open(path, 'r') as original_file:
        content = original_file.read()
      new_content = f">{identity}\n" + content
      output_path = f"{identity}.fna"
      with open(output_path, 'w') as output_file:
          output_file.write(new_content)