import argparse
from pyfaidx import Fasta

parser = argparse.ArgumentParser(description='Convert CNVs into a VCF')
parser.add_argument('-v', help='CNVs from Control-FREEC', required=True)
parser.add_argument('-r', help='reference fasta', required=True)
parser.add_argument('-o', help='output VCF file', required=True)
args = parser.parse_args()

# open connection to reference fasta
fa = Fasta(args.r)

# VCF header
header = '##fileformat=VCFv4.2\n'
header += '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
header += '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Size of the SV in bp">\n'
header += '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
header += '##INFO=<ID=CN,Number=1,Type=Integer,Description="Predicted copy number">\n'
header += '##INFO=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
header += '##INFO=<ID=UNCERTAINTY,Number=1,Type=Float,Description="Percentage of uncertainty of the predicted genotype. It is NOT the uncertainty of a CNA">\n'
header += '##ALT=<ID=DEL,Description="Deletion">\n'
header += '##ALT=<ID=DUP,Description="Duplication">\n'
header += '##ALT=<ID=LOH,Description="Loss of heterozygosity">\n'
header += '##source=ControlFREEC\n'
header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'

# Start writing the output file
outf = open(args.o, 'w')
outf.write(header)

# write a record for each VCF
inf = open(args.v, 'r')
for line in inf:
    line = line.rstrip().split('\t')
    r_chr = 'chr' + line[0]  # Control-FREEC removes 'chr'...
    r_pos = int(line[1])
    r_end = int(line[2])
    sv_type = 'NA'
    if line[4] == 'loss':
        sv_type = 'DEL'
    elif line[4] == 'gain':
        sv_type = 'DUP'
    elif line[4] == 'neutral':
        sv_type = 'LOH'
    r_info = 'CN={};GT={};UNCERTAINTY={}'.format(line[3], line[5],
                                                 line[6])
    r_info += ';SVTYPE={};END={};SVLEN={}'.format(sv_type, r_end,
                                                  r_end - r_pos)
    ref_base = fa[r_chr][r_pos:(r_pos+1)].seq
    outf.write('{}\t{}\t.\t{}\t<{}>\t.\t.\t{}\n'.format(r_chr, line[1],
                                                        ref_base, sv_type,
                                                        r_info))

outf.close()
inf.close()
