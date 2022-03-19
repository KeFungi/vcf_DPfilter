from Bio import SeqIO
import argparse

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", metavar='input.fasta', help="input .fasta file", type=str)
parser.add_argument("-s", metavar='site.txt', help="tab-separated CHROM POS, each site per line", type=str)
parser.add_argument("-p", metavar='PREFIX', help="output prefix", type=str)
parser.add_argument("--include-masked", help="include lowercase masked regions", action="store_true")

# args = parser.parse_args(['-i', 'reference.fasta', '-s', 'input.vcf.sites', '-p', 'ascbias'])
args = parser.parse_args()

# set IO
ref_file = SeqIO.parse(open(args.i, "r"), "fasta")
sites_file = open(args.s, "r")

basic_bases = ["A", "T", "C", "G"]
bases = basic_bases + ["N", "W", "M", "R", "Y", "K", "S"]
base_counter = {base: 0 for base in bases}
base_counter_var = {base: 0 for base in bases}


def sum_amb_count(base_count):
    A = base_count['A'] + 0.5*base_count['W'] + 0.5*base_count['M'] + 0.5*base_count['R']
    T = base_count['T'] + 0.5*base_count['W'] + 0.5*base_count['Y'] + 0.5*base_count['K']
    C = base_count['C'] + 0.5*base_count['M'] + 0.5*base_count['Y'] + 0.5*base_count['S']
    G = base_count['G'] + 0.5*base_count['R'] + 0.5*base_count['K'] + 0.5*base_count['S']
    return {'A': A, 'T': T, 'C': C, 'G': G}


for ref in ref_file:
    sites_file.seek(0)
    sites_var = \
        [site.strip().split('\t')[1] for site in sites_file if site.strip().split('\t')[0] == ref.name]

    if args.include_masked:
        full_seq = ref.seq.upper()
        var_seq = \
            ''.join([full_seq[int(pos) - 1] for pos in sites_var]).upper()
    else:
        full_seq = ref.seq
        var_seq = \
            ''.join([full_seq[int(pos) - 1] for pos in sites_var])

    for base in bases:
        base_counter[base] = base_counter[base] + full_seq.count(base)
        base_counter_var[base] = base_counter_var[base] + var_seq.count(base)

sites_file.close()

all_sum = sum_amb_count(base_counter)
var_sum = sum_amb_count(base_counter_var)
invar_sum = {base: round(all_sum[base]-var_sum[base]) for base in basic_bases}

with open(args.p + '.stamatakis', 'w') as f:
    f.writelines(f"{invar_sum['A']} {invar_sum['C']} {invar_sum['G']} {invar_sum['T']}")

with open(args.p + '.felsenstein', 'w') as f:
    f.writelines(f"{round((all_sum['A']+all_sum['T']+all_sum['C']+all_sum['G'])-(var_sum['A']+var_sum['T']+var_sum['C']+var_sum['G']))}\n")

