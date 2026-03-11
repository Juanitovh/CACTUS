"""Dynamic path resolution for CACTUS binaries and directories."""

import os
import shutil
from pathlib import Path
from types import SimpleNamespace


def _find_optimizer() -> str:
    """Locate the C++ optimizer binary using a fallback chain.

    1. CACTUS_OPTIMIZER_PATH environment variable
    2. shutil.which("joint_fibre_optimizer") (binary on PATH)
    3. Known build location: <repo_root>/source_fibre_optimization/joint_fibre_optimizer
    4. Legacy location: <repo_root>/cactus_scripts/joint_fibre_optimizer.out
    """
    # 1. Environment variable
    env = os.environ.get("CACTUS_OPTIMIZER_PATH")
    if env and os.path.isfile(env):
        return env

    # 2. On PATH
    on_path = shutil.which("joint_fibre_optimizer")
    if on_path:
        return on_path

    # 3. Known build location (Makefile output)
    repo_root = Path(__file__).resolve().parent.parent.parent
    build_loc = repo_root / "source_fibre_optimization" / "joint_fibre_optimizer"
    if build_loc.is_file():
        return str(build_loc)

    # 4. Legacy location
    legacy = repo_root / "cactus_scripts" / "joint_fibre_optimizer.out"
    if legacy.is_file():
        return str(legacy)

    return str(build_loc)  # Return expected path even if not built yet


cactus_paths = SimpleNamespace(
    cactus_code=str(Path(__file__).resolve().parent),
    global_optimizer=_find_optimizer(),
    optimized_final="optimized_final.txt",
)
