#!/usr/bin/python3

import subprocess
import sys

#.......................................
#             for edit                        
#.......................................
total_files = 1329 #number of files in directory (cd to directory and type "ls -1 | wc -l")
jobs = 83 #number of jobs to submit

bin=int(total_files/jobs)
res = total_files%bin
if(int(sys.argv[1]) > 83):
    print(f'Error: Maximum job number is: {jobs} but input readed is {sys.argv[1]}')
    sys.exit(-1)
def read_file_paths(file_path, start_line, end_line):
    """Read file paths from a text file, starting at `start_line` and ending at `end_line`."""
    file_paths = []
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            if i < start_line:
                continue
            if i >= end_line:
                break
            stripped_line = line.strip()
            if stripped_line:
                file_paths.append(stripped_line)
    return file_paths
    
def run_root_macro(i,macro, file_paths,minbias = sys.argv[2]):
    # Escape double quotes in file paths
    escaped_paths = [path.replace('"', r'\"') for path in file_paths]
    # Convert file paths to a comma-separated string
    file_paths_str = ",".join(escaped_paths)
    
    # Construct the ROOT command
    command = [
        "root", "-b", "-q",
        f"{macro}({i},\"{file_paths_str}\",{minbias})"
    ]
    
    # Print the command for debugging
    #print(f"Running command: {' '.join(command)}")
    
    # Execute the command
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")

if len(sys.argv) !=3:
  print("This python script needs 2 arguments but ", len(sys.argv) -1," was readed")
  sys.exit(-1)
    

macro = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/CorrFunc.cxx"
file_list_path = "/gpfs0/citron/users/bargl/ZDC/lhcf22/ppflow/files.txt"


# Specify the range of lines to read (0-indexed)
if (res!=0 and sys.argv[1]== jobs):
  start_line = int(sys.argv[1]) * bin
  end_line   = start_line + res 
else :  
  start_line = int(sys.argv[1]) * bin  #not included
  end_line = start_line + bin    #insdcluded
print(bin)
# Read file paths from the specified range
file_paths = read_file_paths(file_list_path, start_line, end_line)

# Run the ROOT macro with the file paths
run_root_macro(int(sys.argv[1]),macro, file_paths)
