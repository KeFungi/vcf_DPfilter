import argparse
import pandas as pd
import numpy as np
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
high_filters = []  # collector of records of high cuts
low_filters = []  # collector of records of low cuts
for record in vcf_reader:
    # skip non-snp site if --snps
    if args.snps:
        if not record.is_snp:
            continue
        elif '*' in record.alleles:
            continue

    # skip non-variant sites unless --nonvariant
    if not args.nonvariant:
        if record.num_hom_ref == record.num_called:
            continue

    new_CallData = namedtuple('CallData', record.FORMAT.split(':'))  # default FORMAT field
    row_high_filter = [False] * n_sample  # records for high cuts this site
    row_low_filter = [False] * n_sample  # records for low cuts this site
    for sample in sample_ind:
        # if not missing
        if record.samples[sample].data.DP is not None and record.samples[sample].data.GT != './.' and record.samples[sample].data.GT != '.|.':
            # if meets high cutoff
            if record.samples[sample].data.DP > high_bonds[sample]:
                # make new Call object
                calldata = [missing_gt] + list(record.samples[sample].data[1:])
                record.samples[sample].data = new_CallData(*calldata)
                record.samples[sample].called = False
                # record high cut
                row_high_filter[sample] = True
            # if meets low cutoff
            if record.samples[sample].data.DP < low_bonds[sample]:
                # make new Call object
                calldata = [missing_gt] + list(record.samples[sample].data[1:])
                record.samples[sample].data = new_CallData(*calldata)
                record.samples[sample].called = False
                # record high cut
                row_low_filter[sample] = True
    # append site cut records
    high_filters.append(row_high_filter)
    low_filters.append(row_low_filter)

    # skip non-variant sites unless --nonvariant
    if not args.nonvariant:
        if record.num_hom_ref == record.num_called:
            continue

    # write site if not all missing
    if record.num_called > 0:
        vcf_writer.write_record(record)

bgzip_output.close()  # end file

# summarise cut conuts
ind_high_filter = [filtered/len(high_filters) for filtered in np.array(high_filters).sum(axis=0)]
ind_low_filter = [filtered/len(low_filters) for filtered in np.array(low_filters).sum(axis=0)]
either_filter= [filtered/len(high_filters) for filtered in np.logical_or(np.array(high_filters), np.array(low_filters)).sum(axis=0)]

# print stats to std
sample_IDs.insert(0, "INDV")
low_bonds.insert(0, "low_Cutoff")
high_bonds.insert(0, "high_Cutoff")
ind_high_filter.insert(0, "high_CutRate")
ind_low_filter.insert(0, "low_CutRate")
either_filter.insert(0, "total_CutRate")

print('\t'.join(map(str, sample_IDs)))
print('\t'.join(map(str, low_bonds)))
print('\t'.join(map(str, high_bonds)))
print('\t'.join(map(str, ind_high_filter)))
print('\t'.join(map(str, ind_low_filter)))
print('\t'.join(map(str, either_filter)))
