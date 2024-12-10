 #!/bin/bash -l

search_dir=''
src_files=''
for entry in *.cpp
do
  echo -n "."
  src_files="$src_files"" ""$entry"
done

#  sbatch "$entry"

command="g++ -Wall -g -std=c++11 -fopenmp -O3 $src_files -o jfg"
# command="g++ -g -std=c++11 -fopenmp -O3 $src_files -o jfg"

eval "$command"
echo "Compiled!"
