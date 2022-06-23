import argparse
import sys
from .script_high.other import common
import numpy as np
import scipy.stats

def my_parser(parser):
  #parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  #                                 description="This program determines which genes are differentially expressed "
  #                                             "based on RNA-seq data for two groups of samples. RPKM values are "
  #                                             "computed for each sample, optionally normalized, and Kolmogorov-Smirnov "
  #                                             "test is then applied "
  #                                             "to them to determine significant difference between "
  #                                             "distributions of values of the two groups. Order of optional "
  #                                             "normalizations: "
  #                                             "1) quantile normalization, 2) log transform, 3) z-score transform.",
  #                                 add_help=False)
  #try:
    required_named = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required_named.add_argument("--group_1_count_files", help="list of files with read counts for group 1 (e.g., tumor group). "
                                                              "Each file must have at least two columns: "
                                                              "gene name and count",
                                type=argparse.FileType(), nargs="+", required=True, metavar="COUNT_FILE")

    required_named.add_argument("--group_2_count_files", help="list of files with read counts for group 2 (e.g., normal group). "
                                                              "Each file must have at least two columns: "
                                                              "gene name and count",
                                type=argparse.FileType(), nargs="+", required=True, metavar="COUNT_FILE")

    required_named.add_argument("--gene_lengths", help="file with two columns: gene name and its length. "
                                                       "Only genes listed in this file will be considered in "
                                                       "the computation. Note that RPKM computation is affected by the "
                                                       "set of considered genes", type=argparse.FileType(),
                                metavar="LENGTH_FILE", required=True)

    optional.add_argument("--output_file", help="output file name. Two columns will be written: "
                                                "gene name and the corresponding KS test P-value. Only genes with "
                                                "P-values smaller than the threshold given by --p_value parameter will "
                                                "be written, default is -", type=argparse.FileType('w'),
                          metavar="FILE", default="-")

    optional.add_argument("--output_group_1_rpkm", help="output file name for gene expression in group 1. Two "
                                                        "columns "
                                                        "will be written: gene name and the corresponding median "
                                                        "RPKM "
                                                        "(quantile-normalized if --quantile_normalization is "
                                                        "specified) "
                                                        "across group 1 counts. default is None",
                          type=argparse.FileType('w'),
                          metavar="FILE")

    #added wang march 2021
    optional.add_argument("--output_group_2_rpkm", help="output file name for gene expression in group 2. Tow "
                                                        "columns "
                                                        "will be written: gene name and the corresponding median "
                                                        "RPKM "
                                                        "(quantile-normalized if --quantile_normalization is "
                                                        "specified ) "
                                                        "across group 2 counts. default is None", 
                           type=argparse.FileType('w'),
                            metavar="FILE")


    optional.add_argument("--p_value", help="P-value cutoff for the KS test or T-test, default= 0.05", type=float,
                          metavar="NUMBER", default=0.05)

    optional.add_argument("--quantile_normalization", help="apply quantile normalization to RPKM values of all "
                                                           "experiments, default is False",
                          action="store_true")

    optional.add_argument("--log_transform", help="apply log-transformation (ln(1 + x)) to all values, default is False",
                          action="store_true")

    optional.add_argument("--z_score", help="apply z-score transformation values of each experiment separately, default is False",
                          action="store_true")

    optional.add_argument("--min_fold_change", help="in addition to the KS-test or T-test, check that the minimum fold change in "
                                                    "median RPKM between the two groups is above the specified number. "
                                                    "If quantile normalization is "
                                                    "activated, the quantile normalized values are compared, default is None",
                          metavar="NUMBER", type=float)
    #addded by wang
    optional.add_argument("--min_medianRPKM", help ="check the median of RPKM in each group,"
                                                " and only keep genes with RPKM greater than the minimum value in bith groups."
                                               "If quantile normalization is activated, then the quantile normalized values are checked, default is None",
                          metavar="NUMBER", type=float)

    optional.add_argument("--output_all_values", help="put values (RPKM or their z-scores) for all input datasets as "
                                                      "additional columns in the output file, default is False",
                          action="store_true")
    #added by wang march 2021
    optional.add_argument("--test_method", help="Differential expression test methods: 0 for KS-test, 1 for T-test, default is 0 for KS-test",
                          metavar="NUMBER", type=int,default=0)

    #optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
    return parser

    #args = parser.parse_args()
  #except IOError as err:
  #  print("Command line arguments error:", err, file=sys.stderr)
  #  exit(1)
#  return parser


def read_gene_counts(open_file):
    cols = common.read_tsv(open_file, skip_comments=True)
    gene_counts = {}
    for i in range(len(cols[0])):
        if common.is_number(cols[1][i]):
            gene_name = cols[0][i]
            count = float(cols[1][i])
            if gene_name in gene_counts:
                raise RuntimeError("Duplicate gene name in file " + open_file.name + ": " + gene_name)

            gene_counts[gene_name] = count

    return gene_counts


def run(args):
 length_data = common.read_tsv(args.gene_lengths, skip_comments=True)
 gene_lengths = dict(zip(length_data[0], map(int, length_data[1])))

 #added wang
 group_1_counts = list(map(read_gene_counts, args.group_1_count_files))
 group_2_counts = list(map(read_gene_counts, args.group_2_count_files))

 if args.test_method==0:
   print('KS-test')
 elif args.test_method==1:
   print('T-test')
 for gene_values in group_1_counts + group_2_counts:
    genes_in_file = frozenset(gene_values.keys())
    for missing_gene in frozenset(gene_lengths.keys()) - genes_in_file:
        del gene_lengths[missing_gene]

 if len(gene_lengths) == 0:
    raise RuntimeError("No genes found that are present in all input files")

 # make sorted list of gene names
 gene_names = np.array(sorted(gene_lengths.keys()))
 # array of gene lengths, ordered by name
 gene_lengths = np.array([gene_lengths[n] for n in gene_names], float)

 # remove unused genes, transform values to arrays ordered by gene_names
 group_1_values = []
 for gene_counts in group_1_counts:
    group_1_values.append(np.array([gene_counts[n] for n in gene_names], float))

 group_2_values = []
 for gene_counts in group_2_counts:
    group_2_values.append(np.array([gene_counts[n] for n in gene_names], float))

 # transform counts to RPKM values
 for gene_values in group_1_values + group_2_values:
    # this is an approximation, because it only considers read counts inside genes.
    # it should be proportional to the total number of reads, which is supposed to
    # be used in the RPKM formula. Since we will only use it in the KS test, which
    # is insensitive to linear transformations, it works ok
    #print('RPKM transform')
    total_reads = float(sum(gene_values))
    gene_values += 1  # we add a pseudocount of 1 to resolve ties with 0 reads
    gene_values *= 1e9 / (total_reads * gene_lengths)

 value_kind = "RPKM"

 if args.quantile_normalization:
    # apply quantile normalization to RPKM values (all datasets simultaneously)
    value_kind = "quantile_normalized_" + value_kind
    #print('Quantile normalization')
    value_matrix = np.stack(group_1_values + group_2_values, axis=1)
    common.quantile_normalization(value_matrix)
    for i, gene_values in enumerate(group_1_values + group_2_values):
        np.copyto(gene_values, value_matrix[:, i])

 if args.min_fold_change is not None:
    # test the fold change between median values
    med_group_1 = np.median(np.stack(group_1_values, axis=1), axis=1)
    med_group_2 = np.median(np.stack(group_2_values, axis=1), axis=1)
    med_large = np.maximum(med_group_1, med_group_2)
    med_small = np.minimum(med_group_1, med_group_2)
    fold_change_large_enough = np.zeros_like(med_group_1, bool)
    fold_change_large_enough[(med_small == 0.0) & (med_large > 0)] = True
    non_zero = med_small > 0
    nz_good = (med_large[non_zero] / med_small[non_zero]) > args.min_fold_change
    fold_change_large_enough[non_zero] = nz_good
    #print('minímum Fold change')

 if args.min_medianRPKM is not None:
    #test the minmum median RPKM in each group
    med_group_1 = np.median(np.stack(group_1_values, axis=1), axis=1)
    med_group_2 = np.median(np.stack(group_2_values, axis=1), axis=1)
    med_group_1_large_enough= np.zeros_like(med_group_1,bool)
    med_group_1_large_enough= med_group_1>args.min_medianRPKM 
    med_group_2_large_enough= np.zeros_like(med_group_2,bool)
    med_group_2_large_enough= med_group_2>args.min_medianRPKM
    #print('minimum median RPKM')

 if args.output_group_1_rpkm is not None:
    args.output_group_1_rpkm.write("#gene\t")
    args.output_group_1_rpkm.write("median_of_" + value_kind + "_for_group_1\n")
    for i in range(len(gene_names)):
        med_group_1_val = np.median([values[i] for values in group_1_values])
        args.output_group_1_rpkm.write(gene_names[i] + "\t" + str(med_group_1_val) + "\n")
    #print('Export group 1 median ')

 #added by wang March 2021
 if args.output_group_2_rpkm is not None:
    args.output_group_2_rpkm.write("#gene\t")
    args.output_group_2_rpkm.write("median_of_" + value_kind + "_for_group_2\n")
    for i in range(len(gene_names)):
        med_group_2_val= np.median([values[i] for values in group_2_values])
        args.output_group_2_rpkm.write(gene_names[i] + "\t" + str(med_group_2_val) + "\n")
    #print('Export group 2 median') 

 if args.log_transform:
    # apply log transform to values
    value_kind = "log_of_" + value_kind
    #print('Log transform')
    for gene_values in group_1_values + group_2_values:
        logs = np.log(gene_values + 1)
        np.copyto(gene_values, logs)

 if args.z_score:
    # apply z-score transform to values
    value_kind = "z-score_of_" + value_kind
    #print('Zscore transform')
    for gene_values in group_1_values + group_2_values:
        z_scores = scipy.stats.zscore(gene_values)
        np.copyto(gene_values, z_scores)
 #added wang
 if args.test_method==0:
  args.output_file.write("#gene\tdifferential_expression_KS_test_p_value")
 elif args.test_method==1:
  args.output_file.write("#gene\tdifferential_expression_T_test_p_value")

 if args.output_all_values:
    for file in args.group_1_count_files + args.group_2_count_files:
        args.output_file.write("\t" + file.name + "_" + value_kind)

 args.output_file.write("\n")

 for i in range(len(gene_names)):
    if args.min_fold_change is not None and not fold_change_large_enough[i]:
        continue
    
    #added by wang 
    if args.min_medianRPKM is not None and ( not med_group_1_large_enough[i] and not med_group_2_large_enough[i]):
        continue
 
    group_1_gene_values = [values[i] for values in group_1_values]
    group_2_gene_values = [values[i] for values in group_2_values]
    if args.test_method ==0:
       D, p_value = scipy.stats.ks_2samp(group_1_gene_values, group_2_gene_values)
    elif args.test_method==1:
       D, p_value= scipy.stats.ttest_ind(group_1_gene_values, group_2_gene_values, equal_var=False)
    if p_value < args.p_value:
        args.output_file.write(gene_names[i] + "\t" + str(p_value))
        if args.output_all_values:
            for gene_values in group_1_values + group_2_values:
                args.output_file.write("\t" + str(gene_values[i]))

        args.output_file.write("\n")

if __name__== '__main__':
  args= my_parser(argparse.ArgumentParser('python differential_expression.py')).parse_args()
  run(args)
