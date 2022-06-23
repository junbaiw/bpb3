import argparse
import sys
import common

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description="""This program computes an intersection of one or more VCF files. """
                                             """The variants are matched by position only. """
                                             """VCF header and data fields will be taken from the first file""",
                                 add_help=False)
try:
    required_named = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required_named.add_argument("--input_files", help="list of VCF files to intersect",
                                type=argparse.FileType(), nargs="+", required=True, metavar="VCF_FILE")

    optional.add_argument("--output_file", help="output file name", type=argparse.FileType('w'),
                          metavar="VCF_FILE", default="-")

    optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

    args = parser.parse_args()
except IOError as err:
    print("Command line arguments error:", err, file=sys.stderr)
    exit(1)

first_vcf, header = common.read_vcf(args.input_files[0], return_header=True)
other_vcfs = map(common.read_vcf, args.input_files[1:])


def mut_pos(mut_info):
    return mut_info[0], mut_info[1]


good_muts = dict(zip(map(mut_pos, first_vcf), first_vcf))
for other_vcf in other_vcfs:
    all_pos = frozenset(map(mut_pos, other_vcf))
    for missing_mut in frozenset(good_muts.keys()) - all_pos:
        del good_muts[missing_mut]

common.write_vcf(args.output_file, good_muts.values(), header)
