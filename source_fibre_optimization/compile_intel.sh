search_dir=''
src_files=''
for entry in *.cpp
do
  echo -n "."
  src_files="$src_files"" ""$entry"
done

#  sbatch "$entry"
echo "-----------------"
echo $src_files
echo "-----------------"

# command="icpx -O3 -std=c++11 -fopenmp $src_files -o jfg"

command="icpx -O3 -std=c++11 -qopenmp $src_files -o jfg"

eval "$command"
echo "Compiled!"
