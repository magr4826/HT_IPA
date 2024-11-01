# external imports
import os
import re


def write_RPA_scissor(
    cutoff_screening,
    scissor=0.0,
    direction=[1, 0, 0],
    wrange=[0, 20],
    wsteps=2001,
    broadening=0.1,
):
    """
    Creates and adjusts the input file for a Yambo RPA calculation that uses
    a scissor operator to compute the dielectric function in the current directory.
    The function assumes that the save folder is one folder up.
    INPUT:
        cutoff_screening:   number of G-vectors in the screening, i.e. the energy cutoff in Ry
        scissor:            scissor operator in eV
        direction:          direction along which the dielectric function is calculated
        wrange:             energy range
        wsteps:             number of steps along the x-axis
        broadening:         value of the broadening in eV
    Output
        f_name:             path to the input file (without the ending).
    """

    # input file name
    f_name = (
        f"RPA_dir_{direction[0]:d}{direction[1]:d}{direction[2]:d}_"
        + f"sc_{int(cutoff_screening):d}"
    )

    # create the input file
    os.system(f"yambo -o c -k hartree -F {f_name}.in -Q -V all")

    # read the input file
    with open(f"{f_name}.in", "r") as f:
        RPA_str = f.read()

    # adjust the input file
    RPA_str = re.sub(
        r"NLogCPUs=0", "NLogCPUs=2", RPA_str
    )  # reduce the number of log file
    RPA_str = re.sub(
        r"NGsBlkXd= 1[ \t]+RL", f"NGsBlkXd = {cutoff_screening:d} Ry", RPA_str
    )
    RPA_str = re.sub(
        r"% QpntsRXd\n[ \t]+\d+[ \t]+\|[ \t]+\d+[ \t]+\|",
        "% QpntsRXd\n 1 | 1 |",
        RPA_str,
    )  # it is always only the first q-point, because we do optics ...
    RPA_str = re.sub(
        r"% EnRngeXd\n[ \t]+[0-9]+.[0-9]+[ \t]+\|[ \t]+[0-9]+.[0-9]+[ \t]+\|",
        f"% EnRngeXd\n {wrange[0]:.6f} | {wrange[1]:.6f} |",
        RPA_str,
    )
    RPA_str = re.sub(
        r"% DmRngeXd\n[ \t]+[0-9]+.[0-9]+[ \t]+\|[ \t]+[0-9]+.[0-9]+[ \t]+\|",
        f"% DmRngeXd\n {broadening} | {broadening} |",
        RPA_str,
    )
    RPA_str = re.sub(r"ETStpsXd= \d+", f"ETStpsXd= {wsteps:d}", RPA_str)
    RPA_str = re.sub(r"% XfnQP_E\n[ \t]+\d+.\d+", f"% XfnQP_E\n {scissor:.6f}", RPA_str)
    RPA_str = re.sub(
        r"% LongDrXd\n[ \t]+[0-9]+.[0-9]+[ \t]+\|[ \t]+[0-9]+.[0-9]+[ \t]+\|[ \t]+[0-9]+.[0-9]+[ \t]+\|",
        f"% LongDrXd\n {direction[0]:.6f} | {direction[1]:.6f} | {direction[2]:.6f} | ",
        RPA_str,
    )
    RPA_str = re.sub(r"FFTGvecs=\s+\d+\s+RL\s+", "", RPA_str)
    # write the adjusted input file
    with open(f"{f_name}.in", "w") as f:
        f.write(RPA_str)

    return f_name
