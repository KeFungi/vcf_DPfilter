import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", metavar='input', help="input .phy file", type=str)
parser.add_argument("-o", metavar='output', help="output .fasta file", type=str)

args = parser.parse_args()

f = open(args.i)
f.readline()

fasta = open(args.o, "w")

phy_main = f.readlines()

for line in phy_main:
    chopped = line.rstrip().split(" ")
    seqname = ">" + chopped[0]
    seqcontent = chopped[1]
    fasta.write(seqname + '\n')
    fasta.write(seqcontent + '\n')

fasta.close()



