
from types import SimpleNamespace


cactus_paths_dict = {
    "cactus_code": str("/home/juanluis/CACTUS/cactus_scripts"),
    "global_optimizer": str("/home/juanluis/CACTUS/cactus_scripts/joint_fibre_optimizer.out"),
    "blender": str("/path/to/binary1"),
    "MCDC": str("/path/to/binary1"),
    "optimized_final":str("optimized_final.txt")
}

# Convert the dictionary to an object
cactus_paths = SimpleNamespace(**cactus_paths_dict)
