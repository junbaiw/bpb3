import argparse
import sys
from .script_high.other import common
import argparse

def my_parser(parser): 
#parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#                                 description="""This program extracts regions near transcription start sites of """
#                                             """selected genes.""",
#                                 add_help=False)
#try:
    required_named = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required_named.add_argument("--gene_annotation", help="gene annotation in GTF format from genCode"
                                                          "(e.g. https://www.gencodegenes.org/releases/19.html)",
                                type=argparse.FileType(), required=True, metavar="GTF_FILE")

    optional.add_argument("--selected_genes", help="file with gene names to use. "
                                                   "The first column will be taken as the list of gene names. "
                                                   "If not given, all genes from the annotation file will be usedi, default is None",
                          type=argparse.FileType(),
                          metavar="FILE")

    optional.add_argument("--upstream_size", help="number of base pairs to take upstream of the TSS, default=1000", type=int,
                          metavar="NUMBER", default=1000)

    optional.add_argument("--downstream_size", help="number of base pairs to take downstream of the TSS, default=1000", type=int,
                          metavar="NUMBER", default=1000)

    optional.add_argument("--transform_chrom", help="transform chromosome names between "
                                                    "long (chr15) and short(15) formats, default is False",
                          action="store_true")

    optional.add_argument("--output_file", help="output file name in BED format, default is - ", type=argparse.FileType('w'),
                          metavar="FILE", default="-")

    #optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

#    args = parser.parse_args()
#except IOError as err:
#    print("Command line arguments error:", err, file=sys.stderr)
#    exit(1)
    return parser


def run(args):
 processed = set()
 selected_genes = None
 if args.selected_genes is not None:
    selected_genes = common.read_tsv(args.selected_genes)[0]
    selected_genes, dups = common.uniq_and_duplicates(selected_genes)
    if len(dups) > 0:
        print("Warning: duplicate genes in the selected genes file:", ", ".join(dups), file=sys.stderr)

 dups = set()
 for line in args.gene_annotation:
    line = line.replace("\n", "")
    if line[0] == "#":
        continue

    fields = line.split("\t")
    if fields[2] != "gene":
        continue

    info = fields[8]
    info_fields = info.split(";")
    gene_name = None
    for f in info_fields:
        if f.strip().startswith("gene_name "):
            gene_name = f.split("\"")[1]

    if gene_name is not None and (selected_genes is None or gene_name in selected_genes):
        if gene_name in processed:
            dups.add(gene_name)
        else:
            processed.add(gene_name)

        strand = fields[6]
        positive_strand = strand == "+"
        tss = int(fields[3 if positive_strand else 4])
        chrom = fields[0]
        if args.transform_chrom:
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            else:
                chrom = "chr" + chrom

        begin = tss - (args.upstream_size if positive_strand else args.downstream_size)
        end = tss + (args.downstream_size if positive_strand else args.upstream_size)

        if begin < 0:
            begin = 0

        args.output_file.write(chrom + "\t" + str(begin) + "\t" + str(end) + "\t" +
                               gene_name + ";" + strand + ";up" + str(args.upstream_size) +
                               ";down" + str(args.downstream_size) + "\n")
 print("Export at ", args.output_file.name )
 if len(dups) > 0:
    print("Warning: duplicate genes in the input GTF file:", ", ".join(dups), file=sys.stderr)

if __name__=='__main__':
  args= my_parser(argparse.ArgumentParser('python gene_regions.py')).parse_args()
  run(args)

