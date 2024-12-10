 #!/bin/bash -l

search_dir=''
src_files=''
for entry in *.cpp
do
  echo -n "."
  src_files="$src_files"" ""$entry"
done

#  sbatch "$entry"

command="g++ -O3 -Ofast -std=c++20 -fopenmp $src_files -o jfg"

eval "$command"
echo "Compiled!"
