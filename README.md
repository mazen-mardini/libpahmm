# paHMM-Tree Library for Python and C
paHMM-Tree is short for "pairwise statistical phylogenetic distance estimation 
using pair hidden Markov models". It's a tool developed by Marcin Bogusz and Simon Whelan.
This library, "**libpahmm**" exposes the tool's capabilities to **Python** as well as **C**.

## What is libpahmm?
It's a library for performing **M**ultiple **S**equence **A**lignment (**MSA**). It's used to 
find evolutionary distances between a set of sequences of nucleotides or amino-acids 
(such as DNA, RNA, proteins, etc.). 
The library lets you configure the way distance calculations are made. These distances can 
then be used to create [phylogenetic trees](https://en.wikipedia.org/wiki/Phylogenetic_tree) 
and visualize the genetic relationship between various organisms.

## Why use libpahmm?
- It allows for more accurate distance calculations that take into account multiple substitutions 
per site on the sequences.
- It has good algorithmic time complexity: O(n^2 * L^2) where n is the number of sequences 
and L is their length (assuming all have the same length).
- It's implemented in C++ internally, most work is done by compiled code and not Python.
- It's free and open-source.

## Basic usage
Let's see how libpahmm can be used in Python and C. This following examples will calculate the 
distance between two sequences in fasta file containing RNA-data.

### Using Python
```python
from pahmm import *

# 1. Create a banding estimator object
be = BandingEstimator()

# 2a. (Alternative a) Load a file
be.set_file_input("covid-19_0.fasta")

# 2b. (Alternative b) Load a string
be.set_str_input(b"[sequences in FASTA-format here]")

# 3. Load the desired distance-calculation method
sequences = be.execute_hky85_model()

# 4. Calculate distance between two sequences
sequences.get_distance_from_names(b"seq1", b"seq2")
```

### Using C
1\. Create a banding estimator object:
```c++
EBCBandingEstimator *be = ebc_be_create();
```
2a\. (Alternative a) Load a file:
```c++
ebc_be_set_input_from_file(be, "covid-19_0.fasta");
```
2b\. (Alternative b) Load a string:
```c++
ebc_be_set_input(be, "[sequences in FASTA-format here]");
```
3\. Load the desired distance-calculation method:
```c++
EBCSequences *sequences = ebc_be_execute_hky85_modelv2(be);
```
4\. Calculate distance between two sequences:
```c++
double d = ebc_seq_get_distance(sequences, "seq1", "seq2");
```
5\. Clean up resources:
```c++
ebc_seq_free(sequences);
ebc_be_free(be);
```

## Installation
Python-library installation can easily be done as such:
```shell script
python3 -m pip install pahmm
```

**Windows notice:** The library build-scripts (CMakeList.txt and setup.py) have not been tested 
on Windows, only on Linux and MacOS. Therefore the library *does not have official Windows support* at the time being.

**C-library notice:** Currently there are no binary packages distributed for the C-library. If you wish to use libpahmm
as a C-library, you could compile it yourself by following the instructions under the Compilation-section of this README.

## Compilation
To compile libpahmm you'll need to install CMake. This can be done with the following command on Linux (Ubuntu):
```shell script
sudo apt install cmake
```
On MacOS I'd recommend using Brew for easy installation:
```shell script
brew install cmake
```

Then navigate to "libpahmm/", and enter the following command to build for Python 3:
```shell script
python3 setup.py install
```
Note: Instead of "install" you could write "develop" (installs, but doesn't copy any binaries to 
Python's site-packages folder) or "bdist_egg" to create a binary package, ready to be distributed.

To build and install the C-library, just run:
```shell script
mkdir build
cd build
cmake ..
make
sudo make install
```
To uninstall, look inside install_manifest.txt and delete everything that's listed there on your system. 
This can be done the follwing way. Make sure to check that install_manifest.txt contains the files you
actually want to delete:
```shell script
cat install_manifest.txt
```
If everything you want to remove everything listed in the file, then run the following command:
```shell script
sudo xargs rm < install_manifest.txt
```

## License
Libpahmm and paHMM-Tree are licensed under the GNU General Public License version 3.0 (GPLv3).

Dlib is licensed under Boost Software License - Version 1.0.

## Credits
- Marcin Bogusz and Simon Whelan for the development of paHMM-Tree: https://github.com/marbogusz/paHMM-Tree

- Lars Arvestad for refactoring paHMM-Tree: https://github.com/arvestad/paHMM-dist

- Mazen Mardini for the paHMM-Tree library code and more refactoring of paHMM-Tree: https://github.com/mazen-mardini/libpahmm
