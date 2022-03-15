import argparse
import vcf
import re

def amb_convert(base):
    ambiguities = {"AA": "A", "TT": "T", "CC": "C", "GG": "G",
                   "AT": "W", "TA": "W", "AC": "M", "CA": "M",
                   "AG": "R", "GA": "R", "TC": "Y", "CT": "Y",
                   "TG": "K", "GT": "K", "CG": "S", "GC": "S",
                   "..": "N"}
    return ambiguities.get(base)


def consensus_from_gt(gt_bases):
    if gt_bases is None:
        return "N"
    else:
        return amb_convert(re.sub('[|/]', '', gt_bases))


# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", metavar='input.vcf.gz', help="input vcf.gz file", type=str)
parser.add_argument("-o", metavar='out.phy', help="output phylip file", type=str)
parser.add_argument("--snps", help="ignore non-snp sites", action='store_true')

#args = parser.parse_args(['-i', 'input.vcf.gz', '-o', 'out.phy', '--snps'])
args = parser.parse_args()

# setup IO
vcf_reader = vcf.Reader(filename=args.i)

# read and write sequences vertically (one column per sample) to .tmp file
## setup IO
temp_output = open(args.o + '.tmp', "w")

sample_IDs = vcf_reader.samples
n_sample=len(sample_IDs)
sample_ind = range(n_sample)
n_site=0

for record in vcf_reader:
    # skip MNP
    if not record.is_snp:
        continue

    # skip DELETION
    if args.snps:
        if '*' in record.alleles:
            continue

    # skip all missing site
    if record.num_called == 0:
        continue

    # skip invariable site
    if record.nucl_diversity == 0:
        continue

    site_seq = ''
    n_site = n_site + 1
    for sample in sample_ind:
        sample_seq = consensus_from_gt(record.samples[sample].gt_bases)
        site_seq = site_seq + sample_seq

    temp_output.writelines(site_seq + '\n')

temp_output.close()

## read .tmp file and write .phy
output = open(args.o, "w")
output.writelines(' ' + str(n_sample) + ' ' + str(n_site) + '\n')

for sample in sample_ind:
    sample_seq = ''
    temp_input = open(args.o + '.tmp', "r")
    for line in temp_input:
        sample_seq = sample_seq + line[sample]
    out_line = sample_IDs[sample] + ' ' + sample_seq
    output.writelines(out_line + '\n')

temp_input.close()
output.close()
