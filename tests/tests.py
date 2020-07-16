# Tests to confirm that pahmm yields the same output as the
# pahmm-tree tool.

from initialize import *
from pahmm import *
import os
from typing import List, Union

NUCLEOTIDE_MODELS = ["GTR", "HKY85"]
AMINO_ACID_MODELS = ["JTT", "LG", "WAG"]

# There are separate tests for nucleotides and amino-acids.
# The tests are stored in lists. Each test will run once for
# each sample.
NUCLEOTIDE_TESTS = [
    {
        "model": "GTR", 
        "parameters": [],
        "alpha": None,
        "gamma_rate_categories": None
    },
    {
        "model": "GTR",
        "parameters": [0.1, 0.2, 0.3, 0.4, 0.5],
        "alpha": None,
        "gamma_rate_categories": None
    },
    {
        "model": "GTR",
        "parameters": [0.1, 0.2, 0.3, 0.4, 0.5],
        "alpha": 0.6,
        "gamma_rate_categories": 2
    },
    {
        "model": "HKY85",
        "parameters": [],
        "alpha": None,
        "gamma_rate_categories": None
    },
    {
        "model": "HKY85",
        "parameters": [0.1],
        "alpha": None,
        "gamma_rate_categories": None
    },
    {
        "model": "HKY85",
        "parameters": [0.1],
        "alpha": 0.6,
        "gamma_rate_categories": 2
    },
]

AMINO_ACID_TESTS = [
    {
        "model": "JTT",
        "parameters": [],
        "alpha": None,
        "gamma_rate_categories": None
    },
    {
        "model": "JTT",
        "parameters": [],
        "alpha": 0.6,
        "gamma_rate_categories": 2
    },
    {
        "model": "LG",
        "parameters": [],
        "alpha": None,
        "gamma_rate_categories": None
    },
    {
        "model": "LG",
        "parameters": [],
        "alpha": 0.6,
        "gamma_rate_categories": 2
    },
    {
        "model": "WAG",
        "parameters": [],
        "alpha": None,
        "gamma_rate_categories": None
    },
    {
        "model": "WAG",
        "parameters": [],
        "alpha": 0.6,
        "gamma_rate_categories": 2
    },
]


def read_distmat(distmat_path: str):
    adjacency_list = []
    
    try:
        with open(distmat_path, "rb") as distmat_file:
            seqs_count = int(distmat_file.readline())

            for i in range(seqs_count):
                parts = distmat_file.readline().split()
                seq_name = parts[0]
                distances = list(map(float, parts[1:]))
                adjacency_list.append((seq_name, distances))
    except FileNotFoundError:
        return []

    return adjacency_list


TOOL_MODEL_FLAGS = {
    "GTR":   "--GTR",
    "HKY85": "--HKY",
    "JTT":   "--JTT",
    "LG":    "--LG",
    "WAG":   "--WAG"
}

TOOL_MODEL_PARAM_FLAGS = {
    "GTR":   "--gtr_params",
    "HKY85": "--hky_params"
}


def test_fasta(fasta_path: str, model: str, parameters: List[float],
               alpha: Union[None, float], gamma_rate_categories: Union[None, int]):
    """Tests fasta samples.

    :param fasta_path: The samples path, must be a .fasta-file.
    :param model: The model.
    :param parameters: Parameters for the model.
    :param alpha: The global alpha parameter.
    :param gamma_rate_categories: Gamma rate ccategories.
    :return: A tuple: (Test status, A message)
    """

    # Get paHMM-tree tool distances
    command_line = ["./" + PAHMM_TREE_EXEC_NAME, "--in", fasta_path, TOOL_MODEL_FLAGS[model]]

    if parameters:
        command_line.append(TOOL_MODEL_PARAM_FLAGS[model])
        command_line.extend(list(map(str, parameters)))

    if alpha is not None:
        command_line.append("--estimateAlpha")
        command_line.append("0")
        command_line.append("--initAlpha")
        command_line.append(str(alpha))

    if gamma_rate_categories is not None:
        command_line.append("--rateCat")
        command_line.append(str(gamma_rate_categories))

    execution = run2(command_line, os.getcwd(), verbose=False)

    if execution.returncode:
        return False, "The pahmm-tree tool did not run successfully.\n" + execution.stdout

    distmat_path = fasta_path + ".paHMM-Tree.distmat"
    tool_adjacency_list = read_distmat(distmat_path)

    # Get pahmm library distances and compare
    be = BandingEstimator()
    be.set_file_input(fasta_path)
    be.alpha = alpha
    be.gamma_rate_categories = gamma_rate_categories

    try:
        if model == "GTR":
            seqs = be.execute_gtr_model(*parameters)
        elif model == "HKY85":
            seqs = be.execute_hky85_model(*parameters)
        elif model == "JTT":
            seqs = be.execute_jtt_model()
        elif model == "LG":
            seqs = be.execute_lg_model()
        elif model == "WAG":
            seqs = be.execute_wag_model()
        else:
            seqs = None  # Should never happen.
    except PAHMMError as error:
        return False, str(error)

    # Compare:
    if not seqs:
        return False, "Could not execute pahmm library model."

    if len(seqs) != len(tool_adjacency_list):
        return False, "The number of sequences read do not match."

    for i in range(len(seqs)):
        for j in range(i):
            lib_distance = seqs.get_distance_from_names(tool_adjacency_list[i][0], tool_adjacency_list[j][0])
            tool_distance = tool_adjacency_list[i][1][j]

            if abs(lib_distance - tool_distance) >= 0.00005:
                return False, f"Distance between '{tool_adjacency_list[i][0].decode('ascii')}' and " \
                              f"'{tool_adjacency_list[j][0].decode('ascii')}' did not match.\n" \
                              f"Library yields: {lib_distance}\n" \
                              f"Tool yields: {tool_distance}"

    # Test ran successfully
    return True, ""


def main():
    total_result = True

    for subdir, nucleotide in [("amino_acid", False), ("nucleotide", True)]:
        for dirpath, dirnames, filenames in os.walk("samples" + "/" + subdir):
            for filename in filenames:
                if not filename.endswith(".fasta"):
                    continue

                # Select the right tests
                if nucleotide:
                    tests = NUCLEOTIDE_TESTS
                else:
                    tests = AMINO_ACID_TESTS

                for test in tests:
                    model = test["model"]
                    parameters = test["parameters"]
                    alpha = test["alpha"]
                    gamma_rate_categories = test["gamma_rate_categories"]

                    print(f"Testing {filename} (model={model}, params={parameters}, "
                          f"a={alpha}, cat={gamma_rate_categories}): ", end="")

                    fasta_path = dirpath + "/" + filename
                    result, message = test_fasta(fasta_path, model, parameters, alpha, gamma_rate_categories)

                    if result:
                        print("\033[32m" + "Success" + "\033[39m")
                    else:
                        print("\033[31m" + "Failure" + "\033[39m")

                    if message:
                        print(message)

                    total_result = total_result and result

    if total_result:
        print("All tests ran successfully.")
    else:
        print("Some tests failed.")

    return True


if __name__ == '__main__':
    try:
        initialize(True)
        print("Testing PaHMM...")
        exit(int(not main()))
    except KeyboardInterrupt:
        print("Testing interrupted by user.")
