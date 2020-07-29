# ==============================================================================
#  libpahmm - library for paHMM-Tree, a phylogenetic tree estimator
#
#  Copyright (c) 2020 Mazen Mardini.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#      This program is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#      GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses>.
# ==============================================================================

import setuptools
import distutils.command.build as b
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.sdist import sdist
from setuptools.command.bdist_egg import bdist_egg
import distutils.command.clean as cl
from subprocess import run
from shutil import copyfile, rmtree
from multiprocessing import cpu_count
import platform
from pathlib import Path

source_path = Path(__file__).parent.absolute()
pahmm_build_dir = source_path / "build"
pahmm_install_prefix_dir = pahmm_build_dir / "prefix"

if platform.system() == "Linux":
    lib_extension = "so"
elif platform.system() == "Darwin":
    lib_extension = "dylib"
elif platform.system() == "Windows":
    lib_extension = "dll"
else:
    lib_extension = "so"


def build(force: bool = False):
    """
    Build libpahmm
    """

    if not Path("python/pahmm/libpahmm." + lib_extension).exists() or force:
        if not pahmm_build_dir.exists():
            pahmm_build_dir.mkdir()

        if not pahmm_install_prefix_dir.exists():
            pahmm_install_prefix_dir.mkdir()

        run(["cmake", "-DBUILD_STATIC=OFF", ".."],
            cwd=str(pahmm_build_dir), check=True)
        run(["make", "-j", str(cpu_count())], cwd=str(pahmm_build_dir), check=True)
        copyfile(str(pahmm_build_dir / ("libpahmm." + lib_extension)), "python/pahmm/libpahmm." + lib_extension)


class PreBuildCommand(b.build):
    """Pre-installation for development mode."""

    def run(self):
        build(force=True)
        b.build.run(self)


class PreDevelopCommand(develop):
    """Pre-installation for development mode."""

    def run(self):
        build(force=True)
        develop.run(self)


class PreInstallCommand(install):
    """Pre-installation for installation mode."""

    def run(self):
        build()
        install.run(self)


class PreSdistCommand(sdist):
    """Pre-installation for sdist mode."""

    def run(self):
        build()
        sdist.run(self)


class PreBdistEggCommand(bdist_egg):
    """Pre-installation for binary distribution (egg) mode."""

    def run(self):
        build()
        bdist_egg.run(self)


class PostCleanCommand(cl.clean):
    """Post-routines after cleaning."""

    def run(self):
        cl.clean.run(self)
        self.safe_remove(Path(".eggs"))
        self.safe_remove(Path("build"))
        self.safe_remove(Path("dist"))
        self.safe_remove(Path("python/_pahmm_cffi.abi3.so"))
        self.safe_remove(Path("python/python_pahmm.egg-info"))
        self.safe_remove(Path("python/pahmm/libpahmm." + lib_extension))

    @staticmethod
    def safe_remove(path: Path, preserve_root=False):
        """Safely removes a tree path.

        It does so by ensuring that the path leads to an inode within
        the source directory.
        """
        path = path.absolute()

        if not path.exists():
            return

        relative_path = path.relative_to(source_path)

        if source_path not in path.parents:
            raise RuntimeError("Unsafe operation. Cannot remove inode outside of source directory: " + str(path))

        if path.is_dir():
            print(f"removing '{str(relative_path)}' (and everything under it)")
            rmtree(path)

            if preserve_root:
                path.mkdir()
        else:
            print(f"removing '{str(relative_path)}'")
            path.unlink(missing_ok=True)


if __name__ == "__main__":
    with open("README.md", "r") as fh:
        long_description = fh.read()

    with open('requirements.txt') as f:
        requirements = f.read().splitlines()

    setuptools.setup(
        name="python-pahmm",
        version="0.1.0",
        author="Mazen Mardini",
        author_email="mazen@mengate.se",
        description="Pairwise statistical phylogenetic distance estimation library using "
                    "pair hidden Markov models for Python and C.",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="",
        packages=setuptools.find_packages('python'),
        package_dir={'': 'python'},
        classifiers=[
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Programming Language :: Python :: 3",
            "Programming Language :: C",
            "Natural Language :: English",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Operating System :: POSIX :: Linux",
            "Operating System :: MacOS",
        ],
        setup_requires=requirements,
        cmdclass={
            'build': PreBuildCommand,
            'develop': PreDevelopCommand,
            'install': PreInstallCommand,
            'sdist': PreSdistCommand,
            'bdist_egg': PreBdistEggCommand,
            'clean': PostCleanCommand
        },
        cffi_modules=["auxiliary/python-pahmm-build.py:ffibuilder"],
        install_requires=requirements,
        python_requires='>=3.6',
        package_data={'pahmm': ["libpahmm." + lib_extension]}
    )
