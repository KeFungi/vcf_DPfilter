import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-i", metavar='input', help="input .fasta file", type=str)
parser.add_argument("-o", metavar='output', help="output .nex file", type=str)

args = parser.parse_args()


count = SeqIO.convert(args.i, "fasta", args.o, "nexus", "DNA")
print("Converted %i records" % count)