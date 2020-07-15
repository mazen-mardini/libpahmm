__all__ = ["initialize", "PAHMM_TREE_EXEC_NAME", "run2"]

from subprocess import run, PIPE, CompletedProcess
from pathlib import Path
from typing import Union, List
import sys
from os import PathLike, symlink
import platform


# Path to paHMM-tree directory
# This directory contains: dlib, src, Makefile, etc.
PAHMM_TREE_DIR = Path("paHMM-dist")

# Executable names
if platform.system() != "Windows":
    PAHMM_TREE_EXEC_NAME = "HMM_" + platform.system()
else:
    PAHMM_TREE_EXEC_NAME = "HMM.exe"


def test_dir(path: Path, name: str, verbose: bool):
    if not path.is_dir():
        if verbose:
            print("Initialization failed. This", name, "directory doesn't exist:",
                  path, file=sys.stderr)
        return False

    return True


def run2(command: List[str], cwd: Union[bytes, str, PathLike, None],
         verbose: bool, **kwargs) -> CompletedProcess:
    return run(command,
               cwd=cwd,
               stderr=sys.stderr if verbose else PIPE,
               stdout=sys.stdout if verbose else PIPE, **kwargs)


def symlink2(src: Path, dst: Path, verbose: bool):
    try:
        if dst.is_symlink():
            dst.unlink(missing_ok=True)

        symlink(str(src), str(dst))
        return True
    except FileExistsError:
        if verbose:
            print("Initialization failed. Cannot create symbolic link because", dst,
                  "is a directory.", file=sys.stderr)
        return False


def initialize(verbose: bool = False) -> bool:
    """ Initialize testing suite.
    :return: True iff initialization ran successfully.
    """

    # Compile pure paHMM-dist/paHMM-tree tool
    if verbose:
        print("Compiling paHMM-tree tool...")

    if not test_dir(PAHMM_TREE_DIR, "paHMM-tree", verbose):
        return False
    
    if not Path(PAHMM_TREE_EXEC_NAME).exists():
        last_run_path = PAHMM_TREE_DIR / "last_run_system.txt"
        last_run_system = None

        if last_run_path.is_file():
            with open(last_run_path, "r") as f:
                last_run_system = f.read()

        if last_run_system != platform.system():
            if run2(["make", "clean"], PAHMM_TREE_DIR, verbose).returncode:
                return False

        if run2(["make"], PAHMM_TREE_DIR, verbose).returncode:
            return False

        if platform.system() != "Windows":
            executable_path = PAHMM_TREE_DIR / "paHMM-tree"
            executable_path_new = PAHMM_TREE_DIR / ("paHMM-tree_" + platform.system())
            executable_path.rename(executable_path_new)
            executable_path = executable_path_new
        else:
            executable_path = PAHMM_TREE_DIR / "paHMM-tree.exe"

        symlink2(executable_path, Path(PAHMM_TREE_EXEC_NAME), verbose)

        with open(last_run_path, "w") as f:
            f.write(platform.system())

    if verbose:
        print("Initialization done!")

    return True


if __name__ == '__main__':
    exit(int(not initialize(verbose=True)))
