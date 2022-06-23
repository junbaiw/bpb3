import os
import sys
from .other import common

bin_folder = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".../../../../", "bin"))
#print(os.path.dirname(os.path.abspath(__file__)))
common.check_file_exists(bin_folder)

platform_subfolders = {"darwin": "Mac", "linux": "Linux"}

platform = sys.platform
bin_subfolder = None
for platform_name, subfolder in platform_subfolders.items():
    if platform.startswith(platform_name):
        bin_subfolder = os.path.join(bin_folder, subfolder)
        break

if bin_subfolder is None:
    print("Unknown platform:", platform, ". This package has binary components for the following ",
          "platforms:", ", ".join(platform_subfolders.keys()), file=sys.stderr)
    exit(1)

common.check_file_exists(bin_subfolder)
bayespi_affinity = os.path.join(bin_subfolder, "bayesPI_affinity")
common.check_file_exists(bayespi_affinity)


def compute_affinity_command(fasta_file, pwm_file, background_iterations, output_folder, potential, seed):
    return bayespi_affinity + " -potential=" + ("0" if potential == "none" else "1") + \
           (" -override_potential=" + str(potential) if common.is_number(potential) else "") + \
           " -strand=2 -psam=" + common.quote(pwm_file) + " -seq=" + common.quote(fasta_file) + \
           " -shuffle=" + str(background_iterations) + " -seed=" + str(seed) + \
           " -out=" + common.quote(output_folder)
