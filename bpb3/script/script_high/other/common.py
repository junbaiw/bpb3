import os
import shutil
import sys
import csv
import re
import tempfile
import numpy as np
import io
import logging

nucleotide_standard_sequence = "ACGT"
nucleotide_complements = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}

def check_folder(out_folder):
 if not os.path.exists(out_folder):
    print("Creater , ", out_folder)
    os.mkdir(out_folder)
 else:
    print("Exists , ", out_folder)

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def setup_logger(out_folder, out_name,log_number):
   ''' set up logger files for different processes
   '''
   log_file= os.path.join(out_folder , out_name+'_bpb3.log_'+ str(log_number))
   if not os.path.exists(out_folder):
     os.makedirs(out_folder )
  
   logger= logging.getLogger(str(log_number))
   file_handler=logging.FileHandler(log_file,mode='w')
   formatter= logging.Formatter('%(asctime)s %(levelname)-8s %(message)s', datefmt='%a, %d %b %Y %H:%M:%S')
   file_handler.setFormatter(formatter)  
   logger.setLevel(logging.INFO)
   logger.propagate=False
   logger.addHandler(file_handler)
   
   #logging.addHandler(logging.FileHandler(log_file))
   #logging.basicConfig(level=logging.INFO,
   #                 format='%(asctime)s %(levelname)-8s %(message)s',
   #                 datefmt='%a, %d %b %Y %H:%M:%S',
   #                 filename=log_file,
   #                 filemode='w')
   print("Log file:", log_file)
   #logging print to screen
   #logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
   #loggin print to a file
   #logger1= logging.getLogger(str(log_number))
   #logger1.addHandler(logging.FileHandler(log_file))
   return logger


def reverse_complement(seq):
    return "".join(reversed([nucleotide_complements[c] for c in seq.upper()]))


def read_tsv(tsv_file, skip=0, sep="\t", return_header=False, skip_comments=False):
    if isinstance(tsv_file, io.IOBase):
        data = list(csv.reader(tsv_file, delimiter=sep))
    else:
        with open(tsv_file) as f:
            data = list(csv.reader(f, delimiter=sep))

    if skip_comments:
        data = [row for row in data if row[0][0] != "#"]

    if return_header:
        header = data[0]
        if skip == 0:
            skip = 1

    data = data[skip:]

    if len(data) == 0:
        return []

    num_cols = len(data[0])
    cols = [None] * num_cols
    for c in range(num_cols):
        cols[c] = [row[c] for row in data]

    if return_header:
        return cols, header

    return cols


def write_tsv(rows, tsv_file, prefix=None):
    with open(tsv_file, "w") as f:
        if prefix is not None:
            f.write(prefix + "\n")

        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


def read_fasta(fasta_file_name, max_seqs=None):
    sequences = []
    data = ""
    seq_name = None
    with open(fasta_file_name) as ff:
        for line in ff:
            line = line.replace("\n", "")
            if len(line) == 0:
                continue

            if line[0] == '>':
                if data != "":
                    sequences.append((seq_name, data))
                    data = ""
                    if max_seqs is not None and len(sequences) >= max_seqs:
                        break

                seq_name = line[1:]
            else:
                data += line.strip()

    if data != "":
        sequences.append((seq_name, data))

    return sequences


def write_fasta(seqs, fasta_file_name):
    with open(fasta_file_name, "w") as f:
        for name, seq in seqs:
            f.write(">" + name + "\n" + seq + "\n")


def prepare_result_folder(res_folder):
    if os.path.exists(res_folder):
        print("The result folder,", os.path.abspath(res_folder), ", exists and will be erased")

    clear_folder(res_folder)


def clear_folder(folder):
    if os.path.exists(folder):
        shutil.rmtree(folder)

    os.mkdir(folder)


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def read_lines(file):
    #added wang
    if isinstance(file, (str,bytes)):
        with open(file) as f:
            return [l.replace("\n", "") for l in f.readlines()]
    else:
        return [l.replace("\n", "") for l in file.readlines()]


def write_lines(lines, file_name):
    with open(file_name, "w") as f:
        f.writelines([str(l) + "\n" for l in lines])


def check_file_exists(file):
    if not os.path.exists(file):
        print(sys.stderr, os.path.abspath(file), "does not exist. Make sure that you have all files in the package")
        exit(1)


def quote(string):
    return "\"" + string + "\""


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def sort_by(list_to_be_sorted, values_to_be_used_for_sorting):
    return [x for (y, x) in sorted(zip(values_to_be_used_for_sorting, list_to_be_sorted), key=lambda pair: pair[0])]


def transpose_list(l):
    return zip(*l)


def normalize_string_for_file_name(s):
    return re.sub('[^0-9a-zA-Z]+', '_', s)


def read_vcf(file, return_header=False):
    if isinstance(file, str):
        open_file = open(file)
    else:
        open_file = file

    muts = []
    header = []
    for line in open_file:
        line = line.replace("\n", "")
        if line[0] == "#":
            header.append(line)
            continue

        fields = line.split("\t")
        # fields[0] = "chr" + fields[0]
        fields[1] = int(fields[1])  # position
        muts.append(fields)

    open_file.close()

    if return_header:
        return muts, header
    else:
        return muts


def write_vcf(open_file, vcf_data, header):
    for line in header:
        open_file.write(line + "\n")

    for mut_info in vcf_data:
        open_file.write("\t".join(map(str, mut_info)) + "\n")


def uniq_and_duplicates(a):
    seen = set()
    dups = set()
    for x in a:
        if x not in seen:
            seen.add(x)
        else:
            dups.add(x)

    return seen, dups


def make_unique(string_list):
    seen = set()
    dup_idx = 0
    for i in range(len(string_list)):
        s = string_list[i]
        unique_s = s
        while unique_s in seen:
            unique_s = s + "_dup" + str(dup_idx)
            dup_idx += 1
            string_list[i] = unique_s

        seen.add(unique_s)


def get_fasta_by_bed(bed_file, fasta_name, genome_file, use_name=True):
    if which("bedtools") is None:
        print >>sys.stderr, "Could not find bedtools executable. Make sure bedtools are installed"
        exit(1)

    os.system("bedtools getfasta -fi " + genome_file + " -bed " + bed_file + " -fo " + fasta_name +
              (" -name" if use_name else ""))


def seq_replace(seq, p, n):
    return seq[:p] + n + seq[p + len(n):]


def exact_wilcox_test(x, y, side):
    tf = tempfile.NamedTemporaryFile(delete=False)
    tf.close()
    write_lines(np.concatenate((x, y)), tf.name)
    command = "Rscript --default-packages=stats,utils -e \"data=read.table(\\\"" + tf.name + \
        """\\\"); wilcox.test(data[1:""" + str(len(x)) + ",1], data[" + str(len(x) + 1) + ":" + str(len(x) + len(y)) + \
        """,1], alternative='""" + side + "')\$p.value \" "

    #print command

    wilcox_p = os.popen(command).read()

    os.remove(tf.name)
    return float(wilcox_p.split(" ")[-1])


def quantile_normalization(a):
    order = a.argsort(axis=0)
    a.sort(axis=0)
    mean_vals = a.mean(axis=1)
    for c in range(a.shape[1]):
        a[order[:,c],c] = mean_vals


def seq_has_known(seq):
    seq = seq.upper()
    for nuc in nucleotide_standard_sequence:
        if nuc in seq:
            return True

    return False
