import sys
paths = sys.argv[1] 

with open(paths, 'r') as file_list:
  for line in file_list:
      # Split each line into the file path and identifier
      path, identity = line.strip().split()
      
      # Read the contents of the file
      with open(path, 'r') as original_file:
          content = original_file.read()
      
      # Prepare the new content with the identifier header at the top
      new_content = f">{identity}\n" + content
      
      # Define the output file path with .fna extension
      output_path = f"{identity}.fna"
      
      # Write the new content to the .fna file
      with open(output_path, 'w') as output_file:
          output_file.write(new_content)