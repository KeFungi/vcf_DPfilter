usage: vcf_DPfilter.py [-h] [-i input.vcf.gz] [-o output.vcf.gz]
                       [-c cutout.tsv] [--snps] [--nonvariant]

optional arguments:
  -h, --help        show this help message and exit
  -i input.vcf.gz   input vcf.gz file
  -o output.vcf.gz  output vcf.gz file
  -c cutout.tsv     cutoff tsv file; columns=sample_ID, low cutoff, high cutoff
  --snps            ignore non-snp sites
  --nonvariant      include non-variable sites
  
