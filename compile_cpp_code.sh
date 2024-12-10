 #!/bin/bash -l

search_dir=''
src_files=''
echo "Takes a couple seconds..."
for entry in source_fibre_optimization/*.cpp
do
    echo "$entry"
  echo -n "."
  src_files="$src_files"" ""$entry"
done

#  sbatch "$entry"

command="g++ -O3 -Ofast -std=c++20 -fopenmp $src_files -o cactus_scripts/joint_fibre_optimizer.out"

eval "$command"
echo "Compiled!"
