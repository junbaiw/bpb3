import os
import urllib.request
import glob
import numpy as np
from .script_high.other import common
import argparse
#from bayespi_bar2_pipeline import melanoma_folder, patient_data_folder, normal_rnaseq_folder, source_data_folder, \
#    genome_folder, genome_fasta_file, gene_annotation_file


#genome_link = \
#    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
#annotation_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
#skin_cancer_data_link = "https://junbaiw.github.io/BayesPI-BAR2/" \
#                        "skin_cancer_icgc_and_normal_melanocytes_geo_data.tgz"
#rnaseq_file = os.path.join(source_data_folder, "exp_seq.tsv")
#rnaseq_file="exp_seq.tsv"


def my_parser(parser):
  required_named= parser.add_argument_group('required arguments')
  optional = parser.add_argument_group('optional arguments')

  required_named.add_argument("--out_data_folder", help = "Export data file folder", 
                              type=str,required=True,metavar="FOLDER" )

  required_named.add_argument('--patient_data_folder', help="patient data folder",
                              type=str,required=True, metavar="FOLDER")

  required_named.add_argument('--normal_rnaseq_folder', help='normal sample RNA-seq data folder, where the file namee shall be *.gene_counts.tsv if they will be selected for analysis',
                              type=str,required=True, metavar="FOLDER")

  required_named.add_argument('--source_data_folder', help='source data folder for downloaed files from ICGC data portal.',
                              type=str, required=True, metavar="FOLDER")

  required_named.add_argument('--genome_folder', help='Genome folder that contains genome information such as fasta fikle and gene annotation files et al',
                              type=str, required=True, metavar="FOLDER")

  required_named.add_argument('--genome_fasta_file', help='Genome fasta file',
                              type=str, required=True, metavar="FILE")

  required_named.add_argument('--gene_annotation_file', help='Gene annotation file from genCode with gtf format ',
                              type=str, required=True, metavar="FILE")

  required_named.add_argument('--rnaseq_file', help='RNA-seq file name for tumor samples ',
                               type=str, required=True, metavar="FILE" )

  optional.add_argument('--genome_link', help='Genome sequence link address, default is ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz',
                              type=str, metavar="WEB", 
                       default=  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz" )

  optional.add_argument('--annotation_link', help='Gene annotation link address, default is ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz',
                              type=str, metavar="WEB", 
                        default= "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz" )

  optional.add_argument('--skin_cancer_data_link', help='Skin cancer data link address, default is https://junbaiw.github.io/BayesPI-BAR2/" \
                             "skin_cancer_icgc_and_normal_melanocytes_geo_data.tgz ',
                              type=str, metavar="WEB", 
                         default= "https://junbaiw.github.io/BayesPI-BAR2/" \
                             "skin_cancer_icgc_and_normal_melanocytes_geo_data.tgz" )
  return parser

def dl_progress_maker(rem_file):
    def dl_progress(count, block_size, total_size):
        percent = int(count * block_size * 100 / total_size)
        print( "\rdownloading", rem_file, \
            "...[%d / %d Mb] %d%%" % ((count * block_size) / 1000000, total_size / 1000000, percent),)

    return dl_progress


#if __name__ == "__main__":
def run(args):
    if not os.path.exists(args.genome_folder):
        os.mkdir(args.genome_folder)
        print("Genome data will be stored in", args.genome_folder)
    else:
        print("Genome data exists in ", args.genome_folder)

    if not os.path.exists(args.genome_fasta_file):
        print("Reference genome file,", args.genome_fasta_file, ", not found. Downloading it from", args.genome_link)
        downloaded_file = args.genome_fasta_file + ".gz"
        urllib.request.urlretrieve(args.genome_link, downloaded_file, reporthook=dl_progress_maker(args.genome_link))
        print("unpacking...")
        os.system("gunzip " + downloaded_file)
    else:
        print("Genome file exist in ", args.genome_fasta_file)

    if not os.path.exists(args.gene_annotation_file):
        print("GENCODE gene annotation file,", args.gene_annotation_file, ", not found. Downloading it from", args.annotation_link)
        downloaded_file = args.gene_annotation_file + ".gz"
        urllib.request.urlretrieve(args.annotation_link, downloaded_file, reporthook=dl_progress_maker(args.annotation_link))
        print("unpacking...")
        os.system("gunzip " + downloaded_file)
    else:
        print("Gene annotation file exist in ", args.gene_annotation_file)

    if not os.path.exists(args.source_data_folder):
        print("Skin cancer source data is not found in folder", args.source_data_folder, ". Downloading it from", \
            args.skin_cancer_data_link)

        skin_cancer_data_folder = os.path.abspath(os.path.join(args.source_data_folder, ".."))
        print(skin_cancer_data_folder)
        if not os.path.exists(skin_cancer_data_folder):
            os.mkdir(skin_cancer_data_folder)

        downloaded_file = os.path.join(skin_cancer_data_folder, os.path.basename(args.skin_cancer_data_link))
        urllib.request.urlretrieve(args.skin_cancer_data_link, downloaded_file, reporthook=dl_progress_maker(args.skin_cancer_data_link))
        print("unpacking...")
        os.system("tar -zxf " + downloaded_file + " -C " + skin_cancer_data_folder)
    else:
        print("Source data folder exist in ", args.source_data_folder)

    if not os.path.exists(args.out_data_folder):
        os.mkdir(args.out_data_folder)
    else:
        print("Output data exist in ", args.out_data_folder)

    specimen = common.read_tsv(os.path.join(args.source_data_folder, "specimen.tsv"), skip=1)
    specimen_id = np.array(specimen[0], str)
    donor_id = np.array(specimen[4], str)
    specimen_type = np.array(specimen[6], str)
    specimen_interval = np.array(specimen[8], str)  # interval in months between initial diagnosis (?) and taking of
                                                    # the sample, used to find the latest sample of a patient in case
                                                    # there are several
    tumor_specimens = np.array([i for i,item in enumerate(specimen_type) if "tumour" in item.lower()])
    normal_specimens = np.array([i for i,item in enumerate(specimen_type) if "normal" in item.lower()])

    print("Total:", len(specimen_id), ", tumor:", len(tumor_specimens), ", normal:", len(normal_specimens), ", sum:", \
        len(tumor_specimens) + len(normal_specimens))

    tumor_specimen_id = frozenset(specimen_id[tumor_specimens])
    normal_specimen_id = frozenset(specimen_id[normal_specimens])
    donor_tumor = {}
    donor_normal = {}
    for i in range(len(specimen_id)):
        sid = specimen_id[i]
        did = donor_id[i]

        dd = donor_tumor if sid in tumor_specimen_id else donor_normal
        if did not in dd:
            dd[did] = []

        dd[did].append(sid)

    tumor_specimens_to_consider = set()
    donor_to_chosen_specimen = {}
    for did in donor_id:
        if did not in donor_tumor:
            print("No tumor specimens for", did)
        else:
            if len(donor_tumor[did]) == 1:
                spec = donor_tumor[did][0]
            else:
                # there are several specimens for this donor, find the one with greatest specimen_interval and use it
                intervals = specimen_interval[np.in1d(specimen_id, donor_tumor[did])]
                intervals_id = sorted(zip(map(int, intervals), donor_tumor[did]), reverse=True)
                spec = intervals_id[0][1]

            tumor_specimens_to_consider.add(spec)
            donor_to_chosen_specimen[did] = spec

        if did not in donor_normal:
           print("No normal specimens for", did)

    print("Considering", len(tumor_specimens_to_consider), "specimens:", ", ".join(sorted(tumor_specimens_to_consider)))
    print("Donors:", len(np.unique(donor_id)))

    unique_muts_file = os.path.join(args.out_data_folder, "ssm_nodup.tsv")
    if not os.path.exists(unique_muts_file):
        print("Removing duplicate mutations...")
        raw_muts_file = os.path.join(args.source_data_folder, "simple_somatic_mutation.open.tsv")
        seen_muts = set()
        with open(unique_muts_file, "w") as unique_muts:
            with open(raw_muts_file) as raw_muts:
                raw_muts.readline()
                for line in raw_muts:
                    fields = line.replace("\n", "").split("\t")
                    mut_id = fields[0]
                    donor_id = fields[1]
                    specimen_id = fields[3]
                    chrom = fields[8]
                    start = int(fields[9])
                    end = int(fields[10])
                    if start != end:
                        #print "non-SNP", mut_id, ": start:", start, ", end:", end
                        continue

                    assembly = fields[12]
                    if assembly != "GRCh37":
                        print("assembly", assembly, "for", mut_id)
                        continue

                    ref = fields[14]
                    mut_from = fields[15]
                    if len(mut_from) != 1:
                        continue

                    if ref != mut_from:
                        print("ref", ref, "is not the same as from", mut_from, "for", mut_id)
                        continue

                    mut_to = fields[16]
                    if len(mut_to) != 1:
                        continue

                    protocol = fields[33]
                    if protocol != "WGS":
                        #print "protocol:", protocol
                        continue

                    mut_info = specimen_id, chrom, start, ref, mut_to
                    if mut_info in seen_muts:
                        continue

                    seen_muts.add(mut_info)
                    unique_muts.write("\t".join([mut_id, donor_id, specimen_id, chrom, str(start), ref, mut_to]) + "\n")

    if not os.path.exists(args.patient_data_folder):
        print("Splitting the mutation set by patient...")
        os.mkdir(args.patient_data_folder)
        donor_muts = {}
        with open(unique_muts_file) as mf:
            for line in mf:
                mut_id, donor_id, specimen_id, chrom, start, ref, mut_to = line.replace("\n", "").split("\t")
                if specimen_id in normal_specimen_id:
                    print(mut_id, "from normal specimen", specimen_id)
                    continue

                if ref not in "ACGT" or mut_to not in "ACGT" or specimen_id not in tumor_specimens_to_consider:
                    continue

                if donor_id not in donor_muts:
                    donor_muts[donor_id] = []

                donor_muts[donor_id].append((chrom, start, mut_id, ref, mut_to))
        #add wang
        for did, muts in donor_muts.items():
            donor_folder = os.path.join(args.patient_data_folder, did)
            os.mkdir(donor_folder)
            with open(os.path.join(donor_folder, "icgc_mutations.tsv"), "w") as f:
                for mut_info in muts:
                    f.write("\t".join(mut_info) + "\n")

    #added by junbai for output file
    normal_count_output_folder= os.path.join(args.out_data_folder,os.path.basename(args.normal_rnaseq_folder))
    #normal_count_output_folder = os.path.join(args.out_data_folder, "normal_melanocyte_counts")
    if not os.path.exists(normal_count_output_folder):
        os.mkdir(normal_count_output_folder)

    gene_length_file = os.path.join(normal_count_output_folder, "gene_lengths.tsv")
    if not os.path.exists(gene_length_file):
        print("Preprocessing gene expression data")
        normal_count_files = []
        #added by junbai for input file
        #normal_count_source_files = glob.glob(os.path.join(args.source_data_folder, "normal_melanocyte_rnaseq",
        #                                                   "*.gene_counts.tsv"))
        normal_count_source_files = glob.glob(os.path.join(args.normal_rnaseq_folder,"*.gene_counts.tsv"))
    
        print(len(normal_count_source_files), "normal melanocyte expression datasets")
        for cf in normal_count_source_files:
            data = common.read_tsv(cf, skip=2)
            short_counts = os.path.join(normal_count_output_folder,
                                        os.path.basename(cf).replace(".gene_counts.tsv", ".only_counts.tsv"))

            common.write_tsv(common.transpose_list([data[0], data[6]]), short_counts)
            normal_count_files.append(short_counts)
            if not os.path.exists(gene_length_file):
                common.write_tsv(common.transpose_list([data[0], data[5]]), gene_length_file)

        donor_expression_data = {}
        with open(args.rnaseq_file) as expr:
            expr.readline()
            for line in expr:
                fields = line.replace("\n", "").split("\t")
                did = fields[0]
                specimen = fields[2]
                if specimen not in tumor_specimens_to_consider:
                    continue

                gene_name = fields[7]

                if gene_name == "?":
                    continue

                count = int(fields[9])
                if did not in donor_expression_data:
                    donor_expression_data[did] = []

                donor_expression_data[did].append((gene_name, count))

        donor_count_files = []
        for did, expr_table in donor_expression_data.items():
            donor_folder = os.path.join(args.patient_data_folder, did)
            if not os.path.exists(donor_folder):
                print("There is expression data for donor", did, ", but no mutations. Putting it into \"other\" folder")
                donor_folder = os.path.join(args.patient_data_folder, "other")
                if not os.path.exists(donor_folder):
                    os.mkdir(donor_folder)

            unique_gene_expr_table = {}
            for gene, count in expr_table:
                if gene in unique_gene_expr_table:
                    unique_gene_expr_table[gene] += count
                else:
                    unique_gene_expr_table[gene] = count

            count_file = os.path.join(donor_folder, did + "_expression.tsv")
            common.write_tsv(unique_gene_expr_table.items(), count_file)
            donor_count_files.append(count_file)

        print(len(donor_count_files), "patient expression datasets")

if __name__== '__main__':
  args= my_parser(argparse.ArgumentParser('python preprocess_icgc_data.py')).parse_args()
  run(args)


