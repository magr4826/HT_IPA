"""
This workflow runs an IPA calculation based on a previous IPA convergence calculation.
The intent is to rerun the final IPA calculations with a different value of the broadening.
"""

# external imports
import sys
import os
import time
import numpy as np
import shutil
import pandas as pd
import json
from pymatgen.entries.computed_entries import ComputedStructureEntry

# local imports
from src.utils.calc_data_class import calc_data
import src.utils.yambo_helper as yambo_helper
import src.utils.yambo_write as yambo_write
import src.utils.qe_helper as qe_helper
import src.utils.qe_write as qe_write
import src.utils.qe_runner as qe_runner
import src.utils.basic_utils as basic_utils

# start timing
start_time = time.time()

# base directory
base_dir = str(sys.argv[1])

# number of cores
ncores = int(sys.argv[2])

# get the id from the args that are called with this script
id = str(sys.argv[3])

# workflow directory
wf_dir = os.path.join(os.getcwd(), "yambo_IPA_broadening")

if not os.path.exists(wf_dir):
    os.mkdir(wf_dir)
else:
    shutil.rmtree(wf_dir)
    os.mkdir(wf_dir)
os.chdir(wf_dir)

# initialize structure for a QuantumEspresso calculation
structure, name, ibrav = qe_helper.qe_init_structure(id, base_dir)


# check if convergence has been carried out
convergence_flag, calc_data_conv = qe_runner.qe_convergence_checker(
    id, wf_dir, ncores, conv_dir="qe_convergence_SG15"
)

# set the broadening and range parameters
wrange_min = 0
wrange_max = 20
cutoff_screening = int(0)
broadening = 0.3

with open(os.path.join(base_dir, "database", id + ".json"), "r") as file:
    entry = json.load(file)
cse = ComputedStructureEntry.from_dict(entry)
# check if the compound is a metal
# This also fails if the convergence algorithm hasn't been run before
if np.isnan(cse.parameters["ipa_eps_kppa"]):
    print(f"{id} is a metal!")
    params_up = {"ipa_eps_kppa": np.nan, "ipa_nbands": np.nan, "ipa_broad": broadening}
    data_up = {
        f"ipa_epsI_0": np.nan,
        f"ipa_epsR_0": np.nan,
        f"ipa_epsI_1": np.nan,
        f"ipa_epsR_1": np.nan,
        f"ipa_epsI_2": np.nan,
        f"ipa_epsR_2": np.nan,
        "ipa_direct_gap": 0,
        "ipa_indirect_gap": 0,
        "ipa_similarity": np.nan,
    }
    basic_utils.update_json(id, base_dir, data_up=data_up, params_up=params_up)
    shutil.rmtree("out")
    shutil.rmtree("SAVE")
    sys.exit()

converged_kppa = cse.parameters["ipa_eps_kppa"]
kppa = converged_kppa
converged_bands = cse.parameters["ipa_nbands"]
calc_data_nscf = calc_data(
    structure=calc_data_conv.structure,
    name=calc_data_conv.name,
    id=calc_data_conv.id,
    ibrav=calc_data_conv.ibrav,
    calc_type="nscf",
    kppa=converged_kppa,
    pw_cutoff=calc_data_conv.pw_cutoff,
    pseudo=calc_data_conv.pseudo,
)

# get k-meshes for pw-nscf calculation
mp_size = " ".join(str(1 * k) for k in calc_data_nscf.k_points_grid)
kmesh_pw = os.popen(f"kmesh.pl {mp_size} wann").read()
kmesh_float = kmesh_pw.split("\n")
kmesh_float = [line.lstrip().split() for line in kmesh_float][:-1]
kmesh_float = np.array([[float(kp) for kp in row] for row in kmesh_float])
kmesh_float[:, 0] = kmesh_float[:, 0] + 1 / 64
kmesh_float[:, 1] = kmesh_float[:, 1] + 2 / 64
kmesh_float[:, 2] = kmesh_float[:, 2] + 3 / 64
kmesh_float = [np.append(row, 1) for row in kmesh_float]

kmesh_pw = "\n".join(" ".join(str(x) for x in y) for y in kmesh_float)
nk_grid = len(kmesh_float)
kmesh_pw = f"{nk_grid}\n" + kmesh_pw

calc_data_nscf.k_points_grid = kmesh_pw

filename_nscf = qe_runner.qe_pw_run(
    calc_data_nscf,
    qe_write.write_pw,
    ncores=ncores,
    kwargs={
        "occupations": "fixed",
        "symmorphic": True,
        "kpoints_mode": "crystal",
        "n_bands": converged_bands,
    },
)


# delete any previous YAMBO database, if it exists
if os.path.exists("SAVE"):
    shutil.rmtree("SAVE")

# create the yambo database
os.chdir(f"out/{id}.save")
os.system(f"p2y -O ../../")
os.chdir(wf_dir)

# yambo setup
os.system(f"yambo")
path_to_rsetup = os.path.join(os.getcwd(), "r_setup")

# check if the compound is a metal
setup_id = 0
with open("r_setup", "r") as setup_file:
    if "Metallic Bands" in setup_file.read():
        params_up = {
            "ipa_eps_kppa": np.nan,
            "ipa_nbands": np.nan,
            "ipa_broad": broadening,
        }
        data_up = {
            f"ipa_epsI_0": np.nan,
            f"ipa_epsR_0": np.nan,
            f"ipa_epsI_1": np.nan,
            f"ipa_epsR_1": np.nan,
            f"ipa_epsI_2": np.nan,
            f"ipa_epsR_2": np.nan,
            "ipa_direct_gap": 0,
            "ipa_indirect_gap": 0,
            "ipa_similarity": np.nan,
        }
        basic_utils.update_json(id, base_dir, data_up=data_up, params_up=params_up)
        shutil.rmtree("out")
        shutil.rmtree("SAVE")
        sys.exit()

# cartesian directions
cart_dirs = np.eye(3, dtype=int)
row_idx = yambo_helper.get_eps_axis(id, base_dir)
r = 0

# create a reference calculation
f_name = yambo_write.write_RPA_scissor(
    cutoff_screening,
    scissor=0,
    direction=cart_dirs[r, :],
    wrange=[wrange_min, wrange_max],
    wsteps=2001,
    broadening=broadening,
)
os.system(f"mpirun -np {ncores} yambo -F {f_name}.in -J {f_name}")


# read out the imaginary part of the dielectric function
eps_df = pd.read_csv(
    f"o-RPA_dir_100_sc_{cutoff_screening}.eps_q1_inv_rpa_dyson",
    comment="#",
    names=["E", "epsi", "epsr", "epsi_o", "epsr_o"],
    delim_whitespace=True,
)

# delete the dipoles
shutil.rmtree(f_name)

direct_gap, _ = yambo_helper.get_direct_gap_parameters(path_to_rsetup)
indirect_gap, _ = yambo_helper.get_indirect_gap_parameters(path_to_rsetup)

# update the json
epsI = eps_df["epsi"].values
epsR = eps_df["epsr"].values
params_up = {
    "ipa_eps_kppa": kppa,
    "ipa_nbands": converged_bands,
    "ipa_broad": broadening,
}
data_up = {
    f"ipa_epsI_{r}": epsI,
    f"ipa_epsR_{r}": epsR,
    "ipa_direct_gap": direct_gap,
    "ipa_indirect_gap": indirect_gap,
    "ipa_similarity": np.nan,
}
basic_utils.update_json(id, base_dir, data_up=data_up, params_up=params_up)

# save the results
np.savetxt("kppa_conv_ipa_shifted.csv", [kppa])
np.savetxt("gaps_rpa.csv", [direct_gap, indirect_gap])

# run the IPA for the other directions, if necessary
for r in row_idx:
    if r == 0:
        pass
    else:
        f_name = yambo_write.write_RPA_scissor(
            cutoff_screening,
            scissor=0,
            direction=cart_dirs[r, :],
            wrange=[wrange_min, wrange_max],
            wsteps=2001,
            broadening=broadening,
        )
        os.system(f"mpirun -np {ncores} yambo -F {f_name}.in -J {f_name}")
        eps_df = pd.read_csv(
            f"o-{f_name}.eps_q1_inv_rpa_dyson",
            comment="#",
            names=["E", "epsi", "epsr", "epsi_o", "epsr_o"],
            delim_whitespace=True,
        )
        epsI = eps_df["epsi"].values
        epsR = eps_df["epsr"].values
        data_up = {f"ipa_epsI_{r}": epsI, f"ipa_epsR_{r}": epsR}
        basic_utils.update_json(id, base_dir, data_up=data_up)
        # delete the dipoles
        shutil.rmtree(f_name)

shutil.rmtree("out")
shutil.rmtree("SAVE")

# save calculation time
with open("timing.txt", "a+") as f:
    f.write(
        f"{os.path.basename(__file__):<25}  {(time.time() - start_time):7.2f} s  {ncores} cores\n"
    )
