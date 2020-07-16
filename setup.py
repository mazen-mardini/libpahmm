##==============================================================================
## libpahmm - library for paHMM-Tree, a phylogenetic tree estimator
##
## Copyright (c) 2020 Mazen Mardini.
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses>.
##==============================================================================

import setuptools
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.sdist import sdist
from setuptools.command.bdist_egg import bdist_egg
from setuptools.command.bdist_wheel import bdist_wheel
from subprocess import run
from python.paths import *
from shutil import copyfile
from multiprocessing import cpu_count
import sys


if sys.platform.startswith('linux'):
	lib_extension = "so"
elif sys.platform.startswith('darwin'):
	lib_extension = "dylib"

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


class PreDevelopCommand(develop):
	"""Pre-installation for development mode."""
	def run(self):
		build(force = True)
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


class PreBdistWheelCommand(bdist_wheel):
	"""Pre-installation for binary distribution (wheel) mode."""
	def run(self):
		build()
		bdist_wheel.run(self)


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
		description="",
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
			"License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
			"Operating System :: POSIX :: Linux",
			"Operating System :: MacOS",
		],
		setup_requires=requirements,
		cmdclass={
			'develop': PreDevelopCommand,
			'install': PreInstallCommand,
			'sdist': PreSdistCommand,
			'bdist_egg': PreBdistEggCommand,
			'bdist_wheel': PreBdistWheelCommand,
		},
		cffi_modules=["python/python-pahmm-build.py:ffibuilder"],
		install_requires=requirements,
		python_requires='>=3.6',
		package_data={'pahmm': ["libpahmm." + lib_extension]}
	)
