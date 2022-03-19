vcf_DPfilter.py
filter genotype DP in .vcf with customized cutoff of each sample/individual

usage: vcf_DPfilter.py [-h] -i input.vcf.gz -o output.vcf.gz -c cutout.tsv
                       [--snps]

optional arguments:
  -h, --help        show this help message and exit
  -i input.vcf.gz   input vcf.gz file
  -o output.vcf.gz  output vcf.gz file
  -c cutout.tsv     cutoff tsv file; columns=sample, low cutoff, high cutoff
  --snps            ignore non-snp sites


vcf2phylip.py
make relaxed phylip from vcf file

usage: vcf2phylip.py [-h] -i input.vcf.gz -o out.phy [--remove-raxml-invar]
                     [--skip-check] [--recode-vcf out.vcf.gz]

optional arguments:
  -h, --help            show this help message and exit
  -i input.vcf.gz       input vcf.gz file
  -o out.phy            output phylip file
  --remove-raxml-invar  remove sites consider invariable when using --asc-corr
                        in RAxML
  --skip-check          skip SNP, missing, invariable checks; override --snps
  --recode-vcf out.vcf.gz
                        make new .vcf according to filtered sites

ref_ascbias.py 
make raxml ascertainment bias from fasta and variable (.vcf) sites

usage: ref_ascbias.py [-h] [-i input.fasta] [-s site.txt] [-p PREFIX]
                      [--include-masked]

optional arguments:
  -h, --help        show this help message and exit
  -i input.fasta    input .fasta file
  -s site.txt       tab-separated CHROM POS, each site per line
  -p PREFIX         output prefix
  --include-masked  include lowercase masked regions

 