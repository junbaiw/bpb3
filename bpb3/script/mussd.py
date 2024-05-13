import os
import sys
import argparse
import shutil
from .script_high.other import common
import glob
import numpy as np

def my_parser(parser):
#parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
#try:
    required_named = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required_named.add_argument("--patient_mutations", help="list of VCF files with somatic mutations, one per patient",
                                type=str, nargs="+", required=True, metavar="VCF_FILE")

    required_named.add_argument("--genome", help="reference genome in FASTA format", type=argparse.FileType(),
                                required=True,
                                metavar="FASTA_FILE")

    optional.add_argument("--regions", help="BED file with regions of interest. All mutations outside of these "
                                            "regions will be discarded, default is None", type=argparse.FileType(),
                          metavar="BED_FILE")

    optional.add_argument("--result_folder", help="folder to put the results (will be erased), default is result", type=str,
                          metavar="FOLDER", default="result")

    optional.add_argument("--cluster_distance", help="maximum distance between mutations in a hot region, default=30", type=int,
                          default=30, metavar="DISTANCE")

    optional.add_argument("--block_distance", help="maximum distance between hot regions in a block, default=1000", type=int,
                          default=1000, metavar="DISTANCE")

    optional.add_argument("--min_block_size", help="minimum number of mutations in a hot region, default=2", type=int,
                          default=2, metavar="NUMBER")

    optional.add_argument("--min_patients_in_block", help="minimum number of different patients' mutations "
                                                          "in a hot region, can be an integer or a float between 0 and "
                                                          "1, meaning the fraction of total number of patients, default=2",
                          type=float, default=2.0, metavar="NUMBER")

    optional.add_argument("--block_flank",
                          help="number of base pairs to add to each side of the block when extracting sequence, default=25",
                          type=int, default=25, metavar="NUMBER")

    optional.add_argument("--take_blocks_from",
                          help="do not compute mutation blocks, take them from the specified file instead. The file "
                               "format is the same as blocks_summary.tsv in the output",
                          type=argparse.FileType(), metavar="FILE")

    optional.add_argument("--patient_germline_mutations", help="list of VCF files with germline mutations, one per "
                                                               "patient. The order of files must correspond to the "
                                                               "order of --patient_mutations list. The reference "
                                                               "sequences for patients will "
                                                               "contain these mutations. default is None",
                          type=str, nargs="+", metavar="VCF_FILE")

#    optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

#    args = parser.parse_args()
#except IOError as err:
#    print("Command line arguments error:", err, file=sys.stderr)
#    exit(1)
    return parser


def ranges_intersect(a_left, a_right, b_left, b_right):
    return a_left <= b_right and b_left <= a_right

def find_selected_regions_in_range(chrom, pos_from, pos_to, regions_by_chrom):
    return [r for r in regions_by_chrom[chrom] if ranges_intersect(pos_from, pos_to, r[0], r[1])]

def mut_str(mut_info):
    return mut_info[0] + ":" + str(mut_info[1]) + " " + mut_info[2] + " -> " + mut_info[3]

def apply_mutation(seq, rel_pos, mut_info, id, p):
    ref_nucs = mut_info[2]
    alt_nucs = mut_info[3]
    if rel_pos + len(ref_nucs) > len(seq):
        print("Mutation", mut_str(mut_info), "in patient", p,
              "extends beyond its block ", id, ". Increase --block_flank if you want to include it", file=sys.stderr)

        exit(1)

    if seq[rel_pos:(rel_pos + len(ref_nucs))].upper() != ref_nucs.upper():
        print(seq, rel_pos)
        print("Mutation", mut_str(mut_info), "in patient", p,
              "has reference nucleotides",
              ref_nucs, "while the reference genome has", seq[rel_pos:(rel_pos + len(ref_nucs))],
              "at the same position. Make sure to provide the same reference genome as used "
              "in the VCF files")

        exit(1)

    seq = seq[:rel_pos] + alt_nucs.lower() + seq[(rel_pos + len(ref_nucs)):]
    return seq

def run(args):
 args.genome.close()
 regions_by_chrom = None
 if args.regions is not None:
    regions_data = common.read_tsv(args.regions)
    if len(regions_data) == 0:
        print("Warning: regions file,", args.regions.name, "is empty")
    else:
        regions_by_chrom = {}
        for i in range(len(regions_data[0])):
            chrom = regions_data[0][i]
            if chrom not in regions_by_chrom:
                regions_by_chrom[chrom] = []

            if len(regions_data) == 3:
                regions_by_chrom[chrom].append((int(regions_data[1][i]), int(regions_data[2][i])))
            else:
                regions_by_chrom[chrom].append((int(regions_data[1][i]), int(regions_data[2][i]), regions_data[3][i]))

        for chrom, regions in regions_by_chrom.items():
            regions_by_chrom[chrom] = sorted(regions)

        print(sum([len(rs) for chr, rs in regions_by_chrom.items()]), "regions of interest are selected in file",
              args.regions.name)


#def ranges_intersect(a_left, a_right, b_left, b_right):
#    return a_left <= b_right and b_left <= a_right

#def find_selected_regions_in_range(chrom, pos_from, pos_to):
#    return [r for r in regions_by_chrom[chrom] if ranges_intersect(pos_from, pos_to, r[0], r[1])]


 common.prepare_result_folder(args.result_folder)

 num_patients = len(args.patient_mutations)

 all_patients_vcf = set()
 for p in range(num_patients):
    if args.patient_mutations[p] in all_patients_vcf:
        print("VCF", args.patient_mutations[p], "is specified more than once", file=sys.stderr)
        exit(1)

    all_patients_vcf.add(args.patient_mutations[p])

 info = []
 info.append(("patients_count", num_patients))
 print(num_patients, "patients in the dataset")

 patient_chrom_muts_sorted = [{} for i in range(num_patients)]
 all_chrom = set()
 patient_mut_selected_in_regions = [0] * num_patients
 patient_mut_source = [0] * num_patients

 print("Sorting mutations...")
 for p in range(num_patients):
    for mut_info in common.read_vcf(args.patient_mutations[p]):
        chrom = mut_info[0]
        other_info = mut_info[1], mut_info[3], mut_info[4], mut_info[2]  # pos ref alt id
        if chrom not in patient_chrom_muts_sorted[p]:
            patient_chrom_muts_sorted[p][chrom] = []

        patient_chrom_muts_sorted[p][chrom].append(other_info)

    all_chrom = all_chrom.union(patient_chrom_muts_sorted[p].keys())
    for chrom in patient_chrom_muts_sorted[p].keys():
        patient_chrom_muts_sorted[p][chrom].sort(key=lambda x: x[0])
        patient_mut_source[p] += len(patient_chrom_muts_sorted[p][chrom])
        if regions_by_chrom is not None:
            if chrom not in regions_by_chrom:
                patient_chrom_muts_sorted[p][chrom] = []
            else:
                cur_region = 0
                cur_mut = 0
                selected_muts = []
                while cur_region < len(regions_by_chrom[chrom]) and cur_mut < len(patient_chrom_muts_sorted[p][chrom]):
                    mut = patient_chrom_muts_sorted[p][chrom][cur_mut]
                    reg = regions_by_chrom[chrom][cur_region]
                    mut_pos = mut[0]
                    if mut_pos > reg[1]:
                        cur_region += 1
                    elif mut_pos >= reg[0]:
                        selected_muts.append(mut)
                        cur_mut += 1
                    else:
                        cur_mut += 1

                patient_mut_selected_in_regions[p] += len(selected_muts)
                patient_chrom_muts_sorted[p][chrom] = selected_muts

 if regions_by_chrom is not None:
    selected_mutations_file = os.path.join(args.result_folder, "mutations_in_regions.tsv")
    with open(selected_mutations_file, "w") as f:
        for mut_patient, chrom_muts in enumerate(patient_chrom_muts_sorted):
            for chrom, muts in chrom_muts.items():
                for pos, ref, alt, id in muts:
                    f.write("\t".join([chrom, str(pos), ref, alt, id, "patient_" + str(mut_patient)]) + "\n")

 if args.patient_germline_mutations is not None:
    if len(args.patient_germline_mutations) != num_patients:
        print("The number of files with somatic mutations (%d) must be the same as the number of " \
              "files with germline mutations (%d)" % (num_patients, len(args.patient_germline_mutations)),
              file=sys.stderr)

        exit(1)

    print("Loading germline mutations to generate reference sequences for patients...")
    patient_chrom_germline_muts = [{} for i in range(num_patients)]
    for p in range(num_patients):
        for mut_info in common.read_vcf(args.patient_germline_mutations[p]):
            chrom = mut_info[0]
            other_info = mut_info[1], mut_info[3], mut_info[4], mut_info[2]  # pos ref alt id
            if chrom not in patient_chrom_germline_muts[p]:
                patient_chrom_germline_muts[p][chrom] = [[], [], [], []]

            for i in range(4):
                patient_chrom_germline_muts[p][chrom][i].append(other_info[i])

        for chrom in patient_chrom_germline_muts[p].keys():
            patient_chrom_germline_muts[p][chrom][0] = np.array(patient_chrom_germline_muts[p][chrom][0], int)
            for i in range(1, 4):
                patient_chrom_germline_muts[p][chrom][i] = np.array(patient_chrom_germline_muts[p][chrom][i], str)

 print("Computing blocks...")
 all_blocks = []
 all_chrom = frozenset(sorted(all_chrom))
 if 0 < args.min_patients_in_block < 1:
    args.min_patients_in_block = int(args.min_patients_in_block * num_patients)

 for chrom in all_chrom:
    all_muts_on_chrom = []
    for p in range(num_patients):
        if chrom in patient_chrom_muts_sorted[p]:
            for mut_info in patient_chrom_muts_sorted[p][chrom]:
                all_muts_on_chrom.append((mut_info[0], p, (chrom,) + mut_info))

    # first pass: select mutations that are dense enough
    all_muts_on_chrom.sort(key=lambda x: x[0])
    blocks = []
    for pat_mut_info in all_muts_on_chrom:
        if len(blocks) == 0:
            blocks.append([pat_mut_info])
        else:
            last_pos = blocks[-1][-1][0]
            pos = pat_mut_info[0]
            if pos - last_pos <= args.cluster_distance:
                blocks[-1].append(pat_mut_info)
            else:
                blocks.append([pat_mut_info])

    blocks = [b for b in blocks if len(b) >= args.min_block_size]
    blocks = [b for b in blocks if len(frozenset([pat_mut[1] for pat_mut in b])) >= args.min_patients_in_block]

    if len(blocks) == 0:
        continue

    # print ", ".join(map(str, blocks))

    if args.take_blocks_from is None:
        joined_blocks = [blocks[0]]
        for b in blocks[1:]:
            dist_between_blocks = b[0][0] - joined_blocks[-1][-1][0]
            if dist_between_blocks <= args.block_distance:
                joined_blocks[-1].extend(b)
            else:
                joined_blocks.append(b)
    else:
        joined_blocks = blocks

    all_blocks.extend(joined_blocks)

 if args.take_blocks_from is not None:
    block_data = common.read_tsv(args.take_blocks_from, skip=1)
    print(len(block_data[0]), "mutation blocks read from", args.take_blocks_from.name)
    blocks_info = []
    all_blocks_read = []
    for i in range(len(block_data[0])):
        bid = block_data[0][i]
        block_chrom = block_data[1][i]
        seq_start = int(block_data[2][i])
        seq_end = int(block_data[3][i])
        mut_start = seq_start + args.block_flank
        mut_end = seq_end - args.block_flank + 1

        read_block = []
        for block_muts in all_blocks:
            min_mut_pos, patient, mut_info = block_muts[0]
            proto_block_chrom = mut_info[0]
            max_mut_pos, _, _ = block_muts[-1]
            if proto_block_chrom == block_chrom and ranges_intersect(mut_start, mut_end, min_mut_pos, max_mut_pos):
                read_block.extend(pat_mut_info for pat_mut_info in block_muts
                                  if mut_start <= pat_mut_info[0] <= mut_end)

        if len(read_block) > 0:
            all_blocks_read.append(read_block)
            blocks_info.append((bid, block_chrom, seq_start, seq_end, mut_start, mut_end))
        else:
            print("The mutation block", bid, "has no mutations in the supplied dataset. Skipping it.")

    all_blocks = all_blocks_read
 else:
    blocks_info = []
    for block_idx in range(len(all_blocks)):
        block = all_blocks[block_idx]
        block_chrom = block[0][2][0]
        mut_start = block[0][0]
        mut_end = block[-1][0]
        seq_start = mut_start - args.block_flank
        seq_end = mut_end + args.block_flank - 1
        bid = "block_" + str(block_idx) + "_" + block_chrom + "_" + str(mut_start) + "_" + str(mut_end)
        blocks_info.append((bid, block_chrom, seq_start, seq_end, mut_start, mut_end))

 print("Patient list:")
 with open(os.path.join(args.result_folder, "patients_summary.tsv"), "w") as pat_file:
    pat_file.write("patient_id\tsource_vcf_file\ttotal_mutations")
    if regions_by_chrom is not None:
        pat_file.write("\tmutations_in_regions")

    pat_file.write("\tmutations_in_blocks\tnumber_of_blocks_with_mutations\n")

    for p in range(num_patients):
        muts_inside_blocks = len([m for b in all_blocks for m in b if m[1] == p])
        blocks_with_muts = len([b for b in all_blocks if len([m for m in b if m[1] == p]) > 0])
        print("[", p, "]", args.patient_mutations[p], ":", patient_mut_source[p], "mutations total ->")
        pat_file.write(str(p) + "\t" + args.patient_mutations[p] + "\t" + str(patient_mut_source[p]))
        if regions_by_chrom is not None:
            print(patient_mut_selected_in_regions[p], "mutations selected in regions ->")
            pat_file.write("\t" + str(patient_mut_selected_in_regions[p]))

        print(muts_inside_blocks, "mutations in", blocks_with_muts, "blocks")
        pat_file.write("\t" + str(muts_inside_blocks) + "\t" + str(blocks_with_muts) + "\n")

 print("Total loaded mutations:", sum(patient_mut_source))
 info.append(("total_mutations", sum(patient_mut_source)))
 if regions_by_chrom is not None:
    print(", total mutations selected in regions:", sum(patient_mut_selected_in_regions))
    info.append(("regions_count", sum(map(len, regions_by_chrom.values()))))
    info.append(("mutations_in_regions", sum(patient_mut_selected_in_regions)))

 print(", total mutations in blocks:", sum(map(len, all_blocks)))
 info.append(("blocks_count", len(all_blocks)))
 info.append(("mutations_in_blocks", sum(map(len, all_blocks))))


# def mut_str(mut_info):
#    return mut_info[0] + ":" + str(mut_info[1]) + " " + mut_info[2] + " -> " + mut_info[3]


#def apply_mutation(seq, rel_pos, mut_info, id, p):
#    ref_nucs = mut_info[2]
#    alt_nucs = mut_info[3]
#    if rel_pos + len(ref_nucs) > len(seq):
#        print("Mutation", mut_str(mut_info), "in patient", p,
#              "extends beyond its block ", id, ". Increase --block_flank if you want to include it", file=sys.stderr)
#
#        exit(1)

#    if seq[rel_pos:(rel_pos + len(ref_nucs))].upper() != ref_nucs.upper():
#        print(seq, rel_pos)
#        print("Mutation", mut_str(mut_info), "in patient", p,
#              "has reference nucleotides",
#              ref_nucs, "while the reference genome has", seq[rel_pos:(rel_pos + len(ref_nucs))],
#              "at the same position. Make sure to provide the same reference genome as used "
#              "in the VCF files")
#
#        exit(1)
#
#    seq = seq[:rel_pos] + alt_nucs.lower() + seq[(rel_pos + len(ref_nucs)):]
#    return seq


 print(len(all_blocks), "mutation blocks (ordered by number of patients):")
 #added wang
 num_patients_in_blocks =list( map(lambda block: len(frozenset([pat_mut[1] for pat_mut in block])), all_blocks))
 block_order = reversed(np.argsort(num_patients_in_blocks))
 #print(num_patients_in_blocks)

 blocks_summary_file = os.path.join(args.result_folder, "blocks_summary.tsv")
 with open(blocks_summary_file, "w") as bsf:
    bsf.write("block_id\tchrom\tstart_pos\tend_pos\tnumber_of_mutations\tnumber_of_patients\tmutation_distribution")
    if regions_by_chrom is not None:
        bsf.write("\tregion_names")

    bsf.write("\n")
    for block_idx in block_order:
        block = all_blocks[block_idx]
        block_tag, chrom, seq_start, seq_end, min_pos, max_pos = blocks_info[block_idx]

        with open(os.path.join(args.result_folder, block_tag + "_mutations.tsv"), "w") as mutf:
            mutf.write("chrom\tpos\tref\talt\tid\tpatient\n")
            for mut_pos, mut_patient, mut_info in block:
                ref_nucs = mut_info[2]
                alt_nucs = mut_info[3]
                id = mut_info[4]
                mutf.write("\t".join([chrom, str(mut_pos), ref_nucs, alt_nucs, id, "patient_" + str(mut_patient)]) +
                           "\n")

        # print seq_start, seq_end
        print("block", block_idx, "@", chrom, ":", min_pos, "-", max_pos)
        if regions_by_chrom is not None:
            rs = find_selected_regions_in_range(chrom, min_pos, max_pos,regions_by_chrom)
            if len(rs[0]) == 3:
                print("[", ", ".join([r[2] for r in rs]), "]")

        print(",", max_pos - min_pos + 1, "bp,", len(block), "mutations", ",")
        patients_in_block = sorted(frozenset([pat_mut[1] for pat_mut in block]))
        print(len(patients_in_block), "patients")

        block_definitions = []
        for p in patients_in_block:
            # seq_start - 1 because for some reason bedtools result is 1-off
            block_definitions.append((chrom, seq_start - 1, seq_end, block_tag + "_patient_" + str(p)))

        block_def_file = os.path.join(args.result_folder, block_tag + ".bed")
        common.write_tsv(block_definitions, block_def_file)
        reference_fasta = os.path.join(args.result_folder, block_tag + "_ref.fasta")
        common.get_fasta_by_bed(block_def_file, reference_fasta, args.genome.name)

        if not os.path.exists(reference_fasta):
            print("Reference FASTA could not be extracted. Is bedtools installed?", file=sys.stderr)
            exit(1)

        ref_seqs = common.read_fasta(reference_fasta)
        if len(ref_seqs) != len(block_definitions):
            print("Not all block sequences could be extracted", file=sys.stderr)
            exit(1)

        alt_seqs = []
        germline_seqs = []
        for s in range(len(ref_seqs)):
            id = block_definitions[s][3]
            p = int(id.split("_")[-1])
            chrom, min_pos, max_pos = block_definitions[s][:3]
            min_pos += 1
            name, seq = ref_seqs[s]

            if args.patient_germline_mutations is not None:
                seq = ref_seqs[s][1]
                if chrom in patient_chrom_germline_muts[p]:
                    germline_muts_in_block = (patient_chrom_germline_muts[p][chrom][0] >= min_pos) & \
                                             (patient_chrom_germline_muts[p][chrom][0] <= max_pos)

                    for i in np.flatnonzero(germline_muts_in_block):
                        mut_pos = patient_chrom_germline_muts[p][chrom][0][i]
                        ref_nucs = patient_chrom_germline_muts[p][chrom][1][i]
                        alt_nucs = patient_chrom_germline_muts[p][chrom][2][i]
                        mut_id = patient_chrom_germline_muts[p][chrom][3][i]
                        rel_pos = mut_pos - min_pos
                        seq = apply_mutation(seq, rel_pos, [chrom, mut_pos, ref_nucs, alt_nucs, mut_id], id, p)

                germline_seqs.append((name, seq))

            for mut_pos, mut_patient, mut_info in block:
                if mut_patient == p:
                    rel_pos = mut_pos - min_pos
                    seq = apply_mutation(seq, rel_pos, mut_info, id, p)

            alt_seqs.append((name, seq))

        alternate_fasta = os.path.join(args.result_folder, block_tag + "_alt.fasta")
        common.write_fasta(alt_seqs, alternate_fasta)

        if args.patient_germline_mutations is not None:
            common.write_fasta(germline_seqs, reference_fasta)

        mutation_distribution = dict([(patient, 0) for patient in patients_in_block])
        for pos, p, m in block:
            mutation_distribution[p] += 1

        bsf.write(block_tag + "\t" + chrom + "\t" + str(seq_start) + "\t" + str(seq_end) + "\t" + str(len(block)) +
                  "\t" + str(len(patients_in_block)) + "\t" + ",".join(map(str, mutation_distribution.values())))

        if regions_by_chrom is not None:
            bsf.write("\t" + ",".join([r[2] for r in rs]))

        bsf.write("\n")

 common.write_tsv(info, os.path.join(args.result_folder, "mutations_summary.tsv"))

if __name__=='__main__':
  args= my_parser(argparse.ArgumentParser('python mussd.py')).parse_args()
  run(args)

