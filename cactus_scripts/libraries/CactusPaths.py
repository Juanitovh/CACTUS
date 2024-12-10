
from types import SimpleNamespace


cactus_paths_dict = {
    "cactus_code": str("/home/jppl/strand_optimization"),
    "global_optimizer": str("~/Documents/juan_parallel_grid/jfg"),
    "blender": str("/path/to/binary1"),
    "MCDC": str("/path/to/binary1"),
    "optimized_final":str("optimized_final.txt")
}

# Convert the dictionary to an object
cactus_paths = SimpleNamespace(**cactus_paths_dict)
