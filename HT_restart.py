"""
This script is used to restart calculations with a given ID.
This only works if the calculation was at least set up correctly
by the main script.
To run calculations with a given ID which have not been set up,
use HT_full_restart.py
"""

# external imports
import os
import numpy as np
import time

# local imports
import src.utils.basic_utils as basic_utils


# number of cores per job
ncores = 32

# requested job memory (does nothing when doing local calculation)
memory = 100000  # (mb)

# select here which workflows to rerun for each material
script_name = ["yambo_IPA_broadening"]

# list of ids to rerun
data = ["agm004466380"]

# name of the calculation folder
calc_folder = "calc"

# get the base directory so each script knows where the source folder is
base_dir = os.getcwd()


def restart_material(mat_id):
    # create the job file and start the job
    try:
        filename = basic_utils.start_calc(
            base_dir,
            "batchjob",
            mat_id,
            script_name,
            file_name="Watcher",
            ncores=ncores,
            memory=memory,
            calc_folder=calc_folder,
            jobname="Watcher",
        )
    except:
        print(mat_id)
    os.chdir(base_dir)
    return True


# loop over all materials
curr_idx = 0
while curr_idx < len(data):
    curr_mat = data[curr_idx]
    os.chdir(base_dir)

    # check for control instructions
    os.chdir("CONTROL")
    if os.path.isfile("STOP"):
        print("Stopping!")
        break
    if os.path.isfile("PAUSE"):
        print("Paused!")
        time.sleep(60)
        continue

    # save the current number, in case the script crashes
    np.savetxt("CURR_IDX", [curr_idx])

    # find out if any jobs can be started
    # this probably only works on our cluster and has to be adapted
    # if the script is used on another cluster
    ncores_total = np.loadtxt("NCORES")
    ncores_curr = os.popen("bjobs -sum -J Watcher").read()
    ncores_curr = ncores_curr.split("\n")
    ncores_curr = ncores_curr[1].split()
    ncores_curr = sum(int(x) for x in ncores_curr)
    max_jobs = int(ncores_total / ncores)
    curr_jobs = int(ncores_curr / ncores)
    print(curr_jobs, max_jobs)

    if curr_jobs < max_jobs:
        restart_material(curr_mat)
        curr_idx = curr_idx + 1
    else:
        print("All slots filled")
        time.sleep(60)
