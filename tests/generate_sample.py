import sys
import random
import argparse

AMINOACID_ALPHABED = "ARNDCQEGHILKMFPSTWYV"
NUCLEIC_ACID_ALPHABED = "ACGT"

def generate_sample(count: int, length: int, alphabet: str) -> str:
    fasta = ""
    for n in range(1, count+1):
        sequence = "".join(random.choices(alphabet, k=length))
        fasta += f">S{n}\n{sequence}\n"
    return fasta

def main():
    parser = argparse.ArgumentParser(description='Generate random FASTA sample.')
    parser.add_argument('count', type=int,
                        help='number of sequences')
    parser.add_argument('length', type=int,
                        help='length of each sequence')
    parser.add_argument('alphabet', choices=["aa", "na"],
                        help='the alphabet (aa="Aminoacid", na="Nucleic Acid")')

    args = parser.parse_args()
    print(generate_sample(args.count, args.length,
                          AMINOACID_ALPHABED if args.alphabet == "aa" else NUCLEIC_ACID_ALPHABED))
    return 0


if __name__ == "__main__":
    sys.exit(main())