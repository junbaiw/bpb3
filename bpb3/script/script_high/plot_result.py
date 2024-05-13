import os
import glob
import collections
import argparse
import sys
import matplotlib as mlt
mlt.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .other import common

def my_parser(parser):
#parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False,
#                                 description="This produces heatmaps like the one in the package paper. For general "
#                                             "use, you may need to edit the layout, font sizes etc. so that everything "
#                                             "is visible")
#try:
    required_named = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required_named.add_argument("--background_folder", help="folder containing the computed background model "
                                                            "produced by background_affinity_changes.py",
                                type=str, required=True, metavar="FOLDER")

    required_named.add_argument("--foreground_folder", help="folder containing the computed delta-dba values for "
                                                            "patient data produced by bayespi_bar.py",
                                type=str, required=True, metavar="FOLDER")

    required_named.add_argument("--significant_pwms", help="file containing the table of PWMs which have significant "
                                                           "affinity changes, "
                                                           "produced by affinity_change_significance_test.py",
                                type=argparse.FileType(), required=True, metavar="FILE")

    optional.add_argument("--output_file", help="output file name. If not specified, the plot will be displayed "
                                                "on the screen", type=argparse.FileType('w'),
                          metavar="FILE")

    optional.add_argument("--title", help="plot title", type=str, metavar="TEXT")

    optional.add_argument("--show",
                          help="show the plot on the screen (this may not look exactly the same as the saved picture)",
                          action="store_true")

    optional.add_argument("--use_tags",
                          help="use the last column of the --significant_pwms file instead of PWM file names "
                               "(the second column) for plot labels",
                          action="store_true")

#    optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
#
#    args = parser.parse_args()
#except IOError as err:
#    print("Command line arguments error:", err, file=sys.stderr)
#    exit(1)
    return parser


#if args.output_file is None:
#    args.show = True


def load_mean_ddba(ddba_file_name):
    mean_ddba_data, header = common.read_tsv(ddba_file_name, return_header=True)
    col_names = mean_ddba_data[0]
    mean_ddba_data = np.array(mean_ddba_data[1:], float)
    header = header[1:]

    mean_ddba_distribution = {}
    for file_i in range(mean_ddba_data.shape[0]):
        file_name = header[file_i]
        if "_disc" in file_name:
            continue

        mean_ddba_distribution[file_name] = mean_ddba_data[file_i, :]

    print(ddba_file_name, ":", mean_ddba_data.shape[1], "mutations, ", mean_ddba_data.shape[0], "PWMs")
    return mean_ddba_distribution, col_names


def plot_block(block_pvals, tf_tags, mut_names, title, right_mark, colorbar):
    #add jbw
    cmap = plt.cm.get_cmap("PuOr_r").copy()
    max_abs = max(np.nanmax(block_pvals), -np.nanmin(block_pvals))
    print(max_abs)
    max_abs = 4
    cmap.set_under("whitesmoke")
    plt.pcolormesh(np.ma.masked_invalid(block_pvals), cmap=cmap, edgecolors="White", vmin=-max_abs, vmax=max_abs)
    plt.ylim(0, block_pvals.shape[0])

    if colorbar:
        cbar_ticks = np.linspace(-np.floor(max_abs), np.floor(max_abs), 5)
        #add jbw
        #cbar = plt.colorbar(cmap=cmap, ticks=cbar_ticks)
        cbar= plt.colorbar(ticks=cbar_ticks)
        #cbar.set_ticklabels(map(str, 10 ** -np.abs(cbar_ticks)))
        cbar.set_ticklabels(list(np.array(10**(-np.abs(cbar_ticks)),str)  ))        

        cbar.ax.set_xlabel('Expected\nprobability\n', fontsize=15)
        cbar.ax.xaxis.set_label_position('top')

        plt.figtext(right_mark, 0.74, "Positive effect", rotation="vertical", fontsize=15)
        plt.figtext(right_mark, 0.34, "Negative effect", rotation="vertical", fontsize=15)

    ax = plt.gca()
    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(block_pvals.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(block_pvals.shape[1]) + 1, minor=False)

    # ax.set_xticklabels(mut_names, minor=False, rotation=40, ha="right", fontsize=8)
    ax.set_xticklabels([])
    ax.set_yticklabels(tf_tags, minor=False,fontsize=6)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        #t.tick1On = False
        #t.tick2On = False
        t.tick1line.set_visible= False
        t.tick2line.set_visible= False

    for t in ax.yaxis.get_major_ticks():
        #t.tick1On = False
        #t.tick2On = False
        t.tick1line.set_visible= False
        t.tick2line.set_visible= False



    plt.xlabel("Patients", labelpad=10, fontsize=15)

    if title is not None:
        plt.title(title, fontsize=15)


def run(args):
 if args.output_file is None:
    args.show = True

 mean_ddba_distribution_background, _ = load_mean_ddba(os.path.join(args.background_folder,
                                                                   "bayespi_bar_result",
                                                                   "ddba",
                                                                   "ddba_integrated_using_mean_ddba.tsv"))

 block = common.read_tsv(args.significant_pwms, skip=1)
 mean_ddba_distribution_real, seq_id = \
    load_mean_ddba(os.path.join(args.foreground_folder, "ddba", "ddba_integrated_using_mean_ddba.tsv"))

 pwm_p_vals = np.array(block[2], float)
 pwm_names = np.array(block[1], str)
 pwm_tags = np.array(block[-1 if args.use_tags else 1], str)

 pwm_positive = np.array(block[0], str) == "Positive"
 mut_names = [sid.split("_")[-1] for sid in seq_id]
 patient_names = mut_names

 pwm_mut_p_vals = np.zeros((len(pwm_p_vals), len(mut_names)), float)
 pwm_mut_ddbas = np.zeros((len(pwm_p_vals), len(mut_names)), float)
 tags = []

 for i in range(len(pwm_p_vals)):
    tags.append(pwm_tags[i])
    pwm = pwm_names[i]
    positive = pwm_positive[i]
    if pwm not in mean_ddba_distribution_background:
        print("PWM file", pwm, "not found in the background distribution: the background is probably wrong")
        continue

    pwm_background = mean_ddba_distribution_background[pwm]
    pwm_background.sort()
    for j in range(len(mut_names)):
        ddba = mean_ddba_distribution_real[pwm][j]
        idx = np.searchsorted(pwm_background, ddba)
        p_val = float(idx) / len(pwm_background)
        if positive:
            p_val = 1 - p_val

        pwm_mut_p_vals[i, j] = p_val
        pwm_mut_ddbas[i, j] = ddba

 tags = np.array(tags, str)

 pwm_mut_p_vals = np.maximum(pwm_mut_p_vals, 1e-10)  # remove zeros from p values before log transform

 log_pwm_mut_p_vals = np.log10(pwm_mut_p_vals)
 pos_p_vals = -log_pwm_mut_p_vals[pwm_positive,]
 pos_tags = tags[pwm_positive]
 neg_p_vals = log_pwm_mut_p_vals[~pwm_positive,]
 neg_tags = tags[~pwm_positive]
 empty = np.zeros((1, len(mut_names)))
 empty[:] = np.nan
 empty_tags = np.array([""], str)
 # print neg_p_vals


 plt.subplots(figsize=(8, 6))
 plt.subplots_adjust(bottom=0.1, top=0.85, wspace=0.1, left=0.15, right=0.95)

 plot_block(np.vstack((neg_p_vals, empty, pos_p_vals)), np.hstack((neg_tags, empty_tags, pos_tags)), patient_names,
           args.title, 0.95, True)

 if args.output_file is not None:
    args.output_file.close()
    plt.savefig(args.output_file.name, transparent=False, dpi=300)

 if args.show:
    plt.show()

if __name__=='__main__':
  args=my_parser(argparse.ArgumentParser('python plot_result.py')).parse_args()
  run(args)


