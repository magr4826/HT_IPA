"""
This workflow runs a simple convergence algorithm to converge
the number of k-points for an IPA calculation.
"""

# external imports
import sys
import os
import re
import time
import numpy as np
import shutil
import pandas as pd

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
wf_dir = os.path.join(os.getcwd(), "yambo_IPA_convergence")

if not os.path.exists(wf_dir):
    os.mkdir(wf_dir)
else:
    shutil.rmtree(wf_dir)
    os.mkdir(wf_dir)
os.chdir(wf_dir)

# initialize structure for a QuantumEspresso calculation
structure, name, ibrav = qe_helper.qe_init_structure(id, base_dir)

# check if the number of electrons in the structure is even
if structure.composition.total_electrons % 2 == 1:
    params_up = {"ipa_eps_kppa": np.nan, "ipa_nbands": np.nan}
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
    sys.exit()

# check if convergence has been carried out
convergence_flag, calc_data_conv = qe_runner.qe_convergence_checker(
    id, wf_dir, ncores, conv_dir="qe_convergence_SG15"
)

calc_data_bands = calc_data_conv
calc_data_bands.calc_type = "nscf"

# find out how many bands are necessary
num_elec = qe_helper.qe_get_electrons(calc_data_bands)
num_bands_pw = np.max([100, int(num_elec / 2 + 10)])
while True:
    filename_bands = qe_runner.qe_pw_run(
        calc_data_bands,
        qe_write.write_pw,
        ncores=ncores,
        kwargs={"occupations": "fixed", "n_bands": num_bands_pw},
    )

    # read out gap
    manifolds, indirect_gap, direct_gap = qe_helper.qe_identify_manifolds(
        f"out/{id}.xml", tolerance=0.1
    )
    _, max_band_idx = yambo_helper.get_eps_band_range(
        f"out/{id}.xml", num_elec, 20 - direct_gap, direct_gap, 0
    )
    if max_band_idx != 0:
        break
    num_bands_pw = num_bands_pw + 50


# get k-meshes for pw-nscf calculation
mp_size = " ".join(str(1 * k) for k in calc_data_conv.k_points_grid)
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

calc_data_nscf = calc_data_conv
calc_data_nscf.calc_type = "nscf"
calc_data_nscf.k_points_grid = kmesh_pw

filename_nscf = qe_runner.qe_pw_run(
    calc_data_nscf,
    qe_write.write_pw,
    ncores=ncores,
    kwargs={
        "occupations": "fixed",
        "symmorphic": True,
        "kpoints_mode": "crystal",
        "n_bands": max_band_idx,
    },
)


# set the broadening and range parameters
wrange_min = 0
wrange_max = 20
cutoff_screening = int(0)
broadening = 0.1
# delete any previous YAMBO database, if it exists
if os.path.exists("SAVE"):
    shutil.rmtree("SAVE")

# create the yambo database
os.chdir(f"out/{id}.save")
os.system(f"p2y -O ../../")
os.chdir(wf_dir)

# yambo setup
os.system(f"yambo")
path_to_rsetup = os.getcwd()

# check if the compound is a metal
setup_id = 0
with open("r_setup", "r") as setup_file:
    if "Metallic Bands" in setup_file.read():
        params_up = {"ipa_eps_kppa": np.nan, "ipa_nbands": np.nan}
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
eps1 = eps_df["epsi"]
sc = 0
delta_kppa = 1500
kppa = calc_data_conv.kppa
calc_data_curr = calc_data(
    structure,
    name,
    id=id,
    ibrav=ibrav,
    kppa=kppa,
    pw_cutoff=calc_data_conv.pw_cutoff,
    pseudo=calc_data_conv.pseudo,
)

# rename the file so the name of the next output file does not get changed
os.rename(
    f"o-RPA_dir_100_sc_{cutoff_screening}.eps_q1_inv_rpa_dyson",
    f"epsilon_{kppa}_{cutoff_screening}",
)

# converge the k-grid
while sc < 0.9:
    # delete the previous YAMBO database
    shutil.rmtree("SAVE")
    calc_data0 = calc_data(
        structure,
        name,
        id=id,
        ibrav=ibrav,
        kppa=kppa,
        pw_cutoff=calc_data_conv.pw_cutoff,
        pseudo=calc_data_conv.pseudo,
    )
    # check whether any k_points were actually added
    # if not, increase kppa again
    while calc_data_curr.k_points_grid == calc_data0.k_points_grid:
        calc_data0 = calc_data(
            structure,
            name,
            id=id,
            ibrav=ibrav,
            kppa=kppa,
            pw_cutoff=calc_data_conv.pw_cutoff,
            pseudo=calc_data_conv.pseudo,
        )
        kppa += delta_kppa
        calc_data_curr = calc_data(
            structure,
            name,
            id=id,
            ibrav=ibrav,
            kppa=kppa,
            pw_cutoff=calc_data_conv.pw_cutoff,
            pseudo=calc_data_conv.pseudo,
        )

    mp_size = " ".join(str(1 * k) for k in calc_data_curr.k_points_grid)
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

    calc_data_nscf = calc_data_conv
    calc_data_nscf.calc_type = "nscf"
    calc_data_nscf.k_points_grid = kmesh_pw

    # run the calculation with the improved k-point grid
    filename_shifted = qe_runner.qe_pw_run(
        calc_data_nscf,
        qe_write.write_pw,
        ncores=ncores,
        kwargs={
            "occupations": "fixed",
            "symmorphic": True,
            "kpoints_mode": "crystal",
            "n_bands": max_band_idx,
        },
    )

    # create the yambo database
    os.chdir(f"out/{id}.save")
    os.system(f"p2y -O ../../")
    os.chdir(wf_dir)

    # yambo setup
    os.system(f"yambo")

    # create a reference calculation
    f_name = yambo_write.write_RPA_scissor(
        cutoff_screening,
        scissor=0,
        direction=cart_dirs[r, :],
        wrange=[wrange_min, wrange_max],
        wsteps=2001,
        broadening=0.1,
    )
    os.system(f"mpirun -np {ncores} yambo -F {f_name}.in -J {f_name}")

    # read out the imaginary part of the dielectric function
    eps_df = pd.read_csv(
        f"o-RPA_dir_100_sc_{cutoff_screening}.eps_q1_inv_rpa_dyson",
        comment="#",
        names=["E", "epsi", "epsr", "epsi_o", "epsr_o"],
        delim_whitespace=True,
    )
    eps2 = eps_df["epsi"]

    os.rename(
        f"o-RPA_dir_100_sc_{cutoff_screening}.eps_q1_inv_rpa_dyson",
        f"epsilon_{kppa}_{cutoff_screening}",
    )

    # compare the two dielectric functions
    sc = 1 - np.trapz(np.abs(eps2.values - eps1.values)) / np.trapz(eps1.values)
    print(f"The similarity is: {sc}")
    eps1 = eps2

    # delete the dipoles
    shutil.rmtree(f_name)

    # check if the material is a metal after increasing the k-point grid:
    all_files = os.listdir()
    r_setups = [file for file in all_files if "r_setup_" in file]
    setup_id = 0
    for filename in r_setups:
        curr_id = int(re.findall(r"\d+", filename)[0])
        if curr_id > setup_id:
            setup_id = curr_id
            path_to_rsetup = filename
    with open(path_to_rsetup, "r") as setup_file:
        if "Metallic Bands" in setup_file.read():
            print("Material is a metal!")
            params_up = {"ipa_eps_kppa": np.nan, "ipa_nbands": np.nan}
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

direct_gap, _ = yambo_helper.get_direct_gap_parameters(path_to_rsetup)
indirect_gap, _ = yambo_helper.get_indirect_gap_parameters(path_to_rsetup)

# update the json
eps_df = pd.read_csv(
    f"epsilon_{kppa}_{cutoff_screening}",
    comment="#",
    names=["E", "epsi", "epsr", "epsi_o", "epsr_o"],
    delim_whitespace=True,
)
epsI = eps_df["epsi"].values
epsR = eps_df["epsr"].values
params_up = {
    "ipa_eps_kppa": kppa,
    "ipa_nbands": max_band_idx,
    "ipa_broad": broadening,
}
data_up = {
    f"ipa_epsI_{r}": epsI,
    f"ipa_epsR_{r}": epsR,
    "ipa_direct_gap": direct_gap,
    "ipa_indirect_gap": indirect_gap,
    "ipa_similarity": sc,
}
basic_utils.update_json(id, base_dir, data_up=data_up, params_up=params_up)

# save the results
np.savetxt("kppa_conv_ipa_shifted.csv", [kppa])
np.savetxt("gaps_rpa.csv", [direct_gap, indirect_gap])

# run the IPA calculation for the other directions, if necessary
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
            broadening=0.1,
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
