
import importlib

# List of unique libraries to check
libraries = [
    "collections", "icecream", "itertools", "libraries", "math", "matplotlib",
    "multiprocessing", "numba", "numpy", "rich", "scipy", "skimage", "sklearn",
    "argparse", "bz2", "copy", "cv2", "gc", "nibabel", "pyvista", "os",
    "pickle", "progressbar", "psutil", "seaborn", "sparse", "subprocess",
    "sys", "time", "tqdm", "warnings"
]

# Check versions of libraries using `__version__`
for lib in libraries:
    try:
        # Dynamically import the library
        module = importlib.import_module(lib)
        
        # Check for __version__
        if hasattr(module, "__version__"):
            print(f"{lib}=='{module.__version__}'")
    except:
        a=1

