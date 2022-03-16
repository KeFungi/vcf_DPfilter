import argparse
import pandas as pd
import vcf
import Bio.bgzf
from collections import namedtuple

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", metavar='input.vcf.gz', help="input vcf.gz file", type=str)
parser.add_argument("-o", metavar='output.vcf.gz', help="output vcf.gz file", type=str)
parser.add_argument("-c", metavar='cutout.tsv', help="cutoff tsv file; columns=sample, low cutoff, high cutoff", type=str)
parser.add_argument("--snps", help="ignore non-snp sites", action='store_true')
parser.add_argument("--nonvariant", help="include non-variable sites", action='store_true')

# args = parser.parse_args(['-i', 'input.vcf.gz', '-o', 'output.vcf.gz', '-c', 'cutout.tsv', '--snps'])
args = parser.parse_args()

# setup IO
DP_cut = pd.read_csv(args.c, sep='\t', header=None)
vcf_reader = vcf.Reader(filename=args.i)
bgzip_output = Bio.bgzf.BgzfWriter(filename=args.o)
vcf_writer = vcf.Writer(bgzip_output, vcf_reader)

# align sample id in .tsv and .vcf
DP_cut_ordered = DP_cut.set_index(0).loc[vcf_reader.samples, :]
sample_IDs = DP_cut_ordered.index.to_list()
n_sample=len(DP_cut_ordered)
sample_ind = range(n_sample)
low_bonds = DP_cut_ordered.iloc[:, 0].tolist()
high_bonds = DP_cut_ordered.iloc[:, 1].tolist()

# iter over sites
missing_gt = './.'
high_filter = [0] * n_sample  # records for high cuts
low_filter = [0] * n_sample  # records for low cuts
n_site = 0
for record in vcf_reader:
    # skip is all not called
    if record.num_called == 0:
        continue

    # skip non-variant sites unless --nonvariant
    if not args.nonvariant:
        if record.num_hom_ref == record.num_called:
            continue

    # skip non-snp site if --snps
    if args.snps:
        if not record.is_snp:
            continue
        elif '*' in record.alleles:
            continue

    new_CallData = namedtuple('CallData', record.FORMAT.split(':'))  # default FORMAT field
    n_site = n_site + 1
    for sample in sample_ind:
        # if not missing
        if record.samples[sample].called:
            # if meets high cutoff
            if record.samples[sample].data.DP > high_bonds[sample]:
                # make new Call object
                calldata = [missing_gt] + list(record.samples[sample].data[1:])
                record.samples[sample].data = new_CallData(*calldata)
                record.samples[sample].called = False
                # record high cut
                high_filter[sample] = high_filter[sample] + 1
            # if meets low cutoff
            if record.samples[sample].data.DP < low_bonds[sample]:
                # make new Call object
                calldata = [missing_gt] + list(record.samples[sample].data[1:])
                record.samples[sample].data = new_CallData(*calldata)
                record.samples[sample].called = False
                # record high cut
                low_filter[sample] = low_filter[sample] + 1

    # skip non-variant sites unless --nonvariant
    if not args.nonvariant:
        if record.num_hom_ref == record.num_called:
            continue

    # write site if not all missing
    if record.num_called > 0:
        vcf_writer.write_record(record)

bgzip_output.close()  # end file

# summarise cut rate
sum_high_filter = [cut/n_site for cut in high_filter]
sum_low_filter = [cut/n_site for cut in low_filter]

# print stats to std
sample_IDs.insert(0, "INDV")
low_bonds.insert(0, "low_Cutoff")
high_bonds.insert(0, "high_Cutoff")
sum_high_filter.insert(0, "high_CutRate")
sum_low_filter.insert(0, "low_CutRate")

print('\t'.join(map(str, sample_IDs)))
print('\t'.join(map(str, low_bonds)))
print('\t'.join(map(str, high_bonds)))
print('\t'.join(map(str, sum_high_filter)))
print('\t'.join(map(str, sum_low_filter)))
