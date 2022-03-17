import argparse
import vcf
import re
import os
import Bio.bgzf

possible_base = {"A": set(["A", "W", "M", "R", "N"]),
                 "T": set(["T", "W", "Y", "K", "N"]),
                 "C": set(["C", "M", "Y", "S", "N"]),
                 "G": set(["G", "R", "K", "S", "N"])}


def amb_convert(base):
    ambiguities = {"AA": "A", "TT": "T", "CC": "C", "GG": "G",
                   "AT": "W", "TA": "W", "AC": "M", "CA": "M",
                   "AG": "R", "GA": "R", "TC": "Y", "CT": "Y",
                   "TG": "K", "GT": "K", "CG": "S", "GC": "S",
                   "..": "N"}
    return ambiguities[base]


def consensus_from_gt(gt_bases):
    if gt_bases is None:
        return "N"
    else:
        return amb_convert(re.sub('[|/]', '', gt_bases))


# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", metavar='input.vcf.gz', help="input vcf.gz file", type=str, required=True)
parser.add_argument("-o", metavar='out.phy', help="output phylip file", type=str, required=True)
parser.add_argument("--remove-raxml-invar", help="remove sites consider invariable when using --asc-corr in RAxML", action='store_true')
parser.add_argument("--skip-check", help="skip SNP, missing, invariable checks; override --snps", action='store_true')
parser.add_argument("--recode-vcf", metavar="out.vcf.gz", help="make new .vcf according to filtered sites")

# args = parser.parse_args(['-i', 'input.vcf.gz', '-o', 'out.phy', '--remove-raxml-invar'])
args = parser.parse_args()

# setup IO
vcf_reader = vcf.Reader(filename=args.i)
temp_output = open(args.o + '.tmp', "w")
if args.recode_vcf is not None:
    bgzip_output = Bio.bgzf.BgzfWriter(filename=args.recode_vcf)
    vcf_writer = vcf.Writer(bgzip_output, vcf_reader)

# read and write sequences vertically (one column per sample) to .tmp file
sample_IDs = vcf_reader.samples
n_sample=len(sample_IDs)
sample_ind = range(n_sample)
n_site=0

for record in vcf_reader:
    # whether conduct site checks
    if not args.skip_check:
        # skip all missing site
        if record.num_called == 0:
            continue

        # skip MNP
        if not record.is_snp:
            continue

        # skip DELETION
        if '*' in record.alleles:
            continue

        # skip invariable site
        if record.heterozygosity == 0:
            continue

    # make site sequence
    site_seq = ''
    for sample in sample_ind:
        sample_seq = consensus_from_gt(record.samples[sample].gt_bases)
        site_seq = site_seq + sample_seq

    # check if variable in RAxML
    if args.remove_raxml_invar:
        if True in [set(site_seq).issubset(base_set) for base_set in possible_base.values()]:
            continue

    # site included
    n_site = n_site + 1

    # write .vcf if requested
    if args.recode_vcf is not None:
        vcf_writer.write_record(record)

    # write to temp file
    temp_output.writelines(site_seq + '\n')

temp_output.close()  # end writing temp file

if args.recode_vcf is not None:
    bgzip_output.close()  # end writing .vcf file

# read .tmp file and write .phy
output = open(args.o, "w")
temp_input = open(args.o + '.tmp', "r")

# header line
output.writelines(f' {n_sample} {n_site}\n')

# sequence lines
for sample in sample_ind:
    sample_seq = ''
    temp_input.seek(0)
    for line in temp_input:
        sample_seq = sample_seq + line[sample]
    out_line = sample_IDs[sample] + ' ' + sample_seq
    output.writelines(out_line + '\n')

temp_input.close()  # end reading temp file
os.remove(args.o + '.tmp')  # remove temp file
output.close()  # end writing .phy
