#from .other import common
from bpb3.script.script_high.other import common
#from bpb3.script.script_high.other.common import which, clear_folder, write_lines, quote, read_lines 
import sys
import time
import subprocess
import os
import argparse

def my_parser(parser):
#    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#                                     description="""This program runs commands from """
#                                                 """the given file in parallel.""", add_help=False)
#    try:
        required_named = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')

        required_named.add_argument("--commands_file", help="file name with commands to run, one command per line",
                                    type=str, required=True, metavar="FILE")

        optional.add_argument("--tmp_folder", help="folder to put the temporary files (will be erased)", type=str,
                              metavar="FOLDER", default="result")

        optional.add_argument("--name", help="name tag for jobs", type=str,
                              metavar="FOLDER", default="job")

#        optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

        add_parallel_options(parser)

#        args = parser.parse_args()
#        args.commands_file = open(args.commands_file)
#    except IOError as err:
##        print("Command line arguments error:", err, file=sys.stderr)
#        exit(1)

#    try:
#        runner = Parallel(args, args.tmp_folder)
#        commands = common.read_lines(args.commands_file)
#        runner.run_commands(commands, args.name)
#    except IOError as err:
#        print("Error running commands:", err, file=sys.stderr)
#        exit(1)
        return parser


def add_parallel_options(parser):
    parallel_options = parser.add_argument_group('parallelization arguments')

    parallel_options.add_argument("--use_cores", help="number of cores to use on one machine",
                                  type=int, default=4, metavar="NUMBER")

    parallel_options.add_argument("--use_slurm", help="use SLURM workload manager to distribute computations across "
                                                      "nodes",
                                  action="store_true")

    parallel_options.add_argument("--slurm_account", help="SLURM account name", type=str, metavar="NAME")

    parallel_options.add_argument("--max_nodes", help="maximum number of nodes to allocate when using SLURM",
                                  type=int, default=8, metavar="NUMBER")


def slurm_present():
    return common.which("sbatch") is not None


def single_node_parallel(commands, tmp_folder, num_processes, name):
    if len(commands) == 0:
        return

    print("Running", len(commands), "jobs in", min(num_processes, len(commands)), "parallel processes")
    common.clear_folder(tmp_folder)

    processes = []
    commands_processed = 0
    while len(commands) > 0:
        while len(processes) >= num_processes:
            for i in range(len(processes)):
                if processes[i].poll() is not None:
                    del processes[i]
                    break

            if len(processes) >= num_processes:
                time.sleep(0.01)

        with open(os.path.join(tmp_folder, name + "_" + str(commands_processed) + ".out"), "wb") as out, \
                open(os.path.join(tmp_folder, name + "_" + str(commands_processed) + ".err"), "wb") as err:
            processes.append(subprocess.Popen(commands[0], shell=True, stdout=out, stderr=err))

        commands_processed += 1
        del commands[0]

    for p in processes:
        #p.wait()
        p.communicate()


def make_sbatch_header(job_name, account_name, max_time_hours, memory_mbytes, num_cpus, num_nodes=1):
    return \
        """#!/bin/bash
        #SBATCH --job-name=""" + job_name + \
        ("" if account_name is None else ("\n#SBATCH --account=" + account_name)) + """
#SBATCH --time=""" + str(max_time_hours) + """:00:00
#SBATCH --mem-per-cpu=""" + str(memory_mbytes) + """M --partition=bigmem 
#SBATCH --cpus-per-task=""" + str(num_cpus) + """

"""
##SBATCH --nodes=""" + str(num_nodes) + """
#"""


def slurm_parallel(commands, tmp_folder, name, num_processes, max_nodes, account):
    if len(commands) == 0:
        return

    use_nodes = max_nodes
    if use_nodes > (len(commands) - 1) / num_processes + 1:
        use_nodes = (len(commands) - 1) / num_processes + 1

    print("Splitting", len(commands), "jobs into", use_nodes, "nodes (up to", num_processes, "processes each)" +
          " using SLURM")

    common.clear_folder(tmp_folder)

    node_commands = [[] for n in range(use_nodes)]
    while len(commands) > 0:
        for n in range(use_nodes):
            node_commands[n].append(commands[0])
            del commands[0]
            if len(commands) == 0:
                break

    done_files = []
    for n in range(use_nodes):
        node_tag = name + "_node_" + str(n)
        commands_file = os.path.abspath(os.path.join(tmp_folder, "commands_for_" + node_tag))
        common.write_lines(node_commands[n], commands_file)

        job_file = os.path.abspath(os.path.join(tmp_folder, "job_" + node_tag))
        with open(job_file, "w") as f:
            f.write(make_sbatch_header(name.replace(".", "") + str(n), account,
                                       8, 2500, min(num_processes, len(node_commands[n]))))

            f.write("\n")
            #f.write("python " + os.path.abspath(__file__) + " --use_cores " + str(num_processes) +
            #        " --commands_file " + commands_file +
            #        " --tmp_folder " + os.path.abspath(os.path.join(tmp_folder, "out_" + node_tag)) +
            #        " --name " + name + "\n")

            #changed by wang
            f.write("bpb3 parallel " + " --use_cores " + str(num_processes) +
                    " --commands_file " + commands_file +
                    " --tmp_folder " + os.path.abspath(os.path.join(tmp_folder, "out_" + node_tag)) +
                    " --name " + name + "\n")

            done_file = os.path.abspath(os.path.join(tmp_folder, node_tag + "_done"))
            f.write("echo done > " + done_file + "\n")
            done_files.append(done_file)

        os.system("sbatch -o " + common.quote(os.path.abspath(os.path.join(tmp_folder, node_tag + "_slurm-%j.out"))) +
                  " " + common.quote(job_file))

    print("Waiting for jobs to finish...")

    while True:
        missing = False
        for df in done_files:
            if not os.path.exists(df):
                missing = True
                break

        if not missing:
            break

        time.sleep(5)


class Parallel:
    def __init__(self, args, tmp_folder):
        self.cores = args.use_cores
        if args.use_slurm:
            if not slurm_present():
                raise IOError("SLURM workload manager not detected")

            self.slurm = True
        else:
            self.slurm = False

        self.slurm_account = args.slurm_account
        self.slurm_max_nodes = args.max_nodes
        self.tmp_folder = tmp_folder

    def run_commands(self, commands, name):
        if self.slurm:
            slurm_parallel(commands, self.tmp_folder, name, self.cores, self.slurm_max_nodes, self.slurm_account)
        else:
            single_node_parallel(commands, self.tmp_folder, self.cores, name)


def run(args):
   args.commands_file= open(args.commands_file)
   try:
     runner= Parallel(args, args.tmp_folder)
     commands = common.read_lines(args.commands_file)
     runner.run_commands(commands, args.name)
   except IOError as err:
     print("Error running commands: ", err, file=sys.stderr)
     exit(1)

if __name__ == "__main__":
  args= my_parser(argparse.ArgumentParser('python parallel.py')).parse_args()
  run(args)


