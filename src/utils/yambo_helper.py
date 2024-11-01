# external imports
import re
import numpy as np

# local imports
from src.utils.basic_utils import ev2ha
import src.utils.qe_helper as qe_helper


def get_direct_gap_parameters(path_to_rsetup):
    """
    Reads the direct gap and the associated k-point + band indices from the r_setup file
    INPUT:
        path_to_rsetup:     path to the r_setup file (full or relative path)
    OUTPUT:
        direct_gap:         direct gap in eV
        kpt_bnd_idx:        k-point and bands of of VB and CB at the direct gap
    """

    # read the input file
    with open(path_to_rsetup, "r") as f:
        setup_str = f.read()

    # array for gap location in the bandstructure
    kpt_bnd_idx = np.zeros(4, dtype=int)

    # get the vbm and cbm index of the direct bandgap
    kpt_bnd_idx[2] = int(
        re.findall("\d+", re.findall(r"Filled Bands[ \t]+:[ \t]+\d+", setup_str)[0])[0]
    )
    kpt_bnd_idx[3] = kpt_bnd_idx[2] + 1

    # get the k-point index of the direct bandgap
    kpt_bnd_idx[0] = int(
        re.findall(
            "\d+",
            re.findall(r"Direct Gap localized at k[ \t]+:[ \t]+\d+", setup_str)[0],
        )[0]
    )
    kpt_bnd_idx[1] = kpt_bnd_idx[0]

    # get the direct bandgap
    direct_gap = float(
        re.findall(
            "\d+.\d+", re.findall(r"Direct Gap[ \t]+:[ \t]+\d+.\d+", setup_str)[0]
        )[0]
    )

    return direct_gap, kpt_bnd_idx


def get_indirect_gap_parameters(path_to_rsetup):
    """
    Reads the indirect gap and the associated k-point + band indices from the r_setup file
    INPUT:
        path_to_rsetup:     path to the r_setup file (full or relative path)
    OUTPUT:
        indirect_gap:       indirect gap in eV
        kpt_bnd_idx:        k-point and bands of of VBM and CBM
    """

    # read the input file
    with open(path_to_rsetup, "r") as f:
        setup_str = f.read()

    # array for gap location in the bandstructure
    kpt_bnd_idx = np.zeros(4, dtype=int)

    # get the vbm and cbm index of the direct bandgap
    kpt_bnd_idx[2] = int(
        re.findall("\d+", re.findall(r"Filled Bands[ \t]+:[ \t]+\d+", setup_str)[0])[0]
    )
    kpt_bnd_idx[3] = kpt_bnd_idx[2] + 1

    # get the k-point index of the indirect bandgap
    matches = re.findall(
        "\d+",
        re.findall(r"Indirect Gap between kpts[ \t]+:[ \t]+\d+[ \t]+\d+", setup_str)[0],
    )
    kpt_idx = [int(m) for m in matches]
    kpt_bnd_idx[0] = kpt_idx[0]
    kpt_bnd_idx[1] = kpt_idx[1]

    # get the direct bandgap
    indirect_gap = float(
        re.findall(
            "\d+.\d+", re.findall(r"Indirect Gap[ \t]+:[ \t]+\d+.\d+", setup_str)[0]
        )[0]
    )

    return indirect_gap, kpt_bnd_idx


def get_eps_band_range(xml_path, num_elec, delta_energy, direct_gap, scissor):
    """
    Obtain the minimum bumber of bands required to calculate the
    dielectric function up to (vbm+direct_gap+delta_energy).
    INPUT:
        xml_path:      path to the QE xml-output
        num_elec:      number of electrons
        delta_energy:  energy difference from direct gap to maximum energy for which to calculate the DF
        direct_gap:    direct gap
        scissor:       optional scissor
    OUTPUT:
        wrange:         energy for which the DF can be calculated correctly
        max_band_idx:   number of necessary bands
    Function returns 0,0 if QE calculation doesn't have enough bands
    """

    # go to the qe output folder and parse the eigenvalues for later on
    k_points, eigenvalues, vbm, cbm = qe_helper.qe_get_eigenvalues(xml_path, num_elec)
    num_kpt = len(k_points)

    # obtain a good band range for a reasonable frequency range of the dielectric function
    max_band_idx = []
    for i in range(num_kpt):
        # check if there are enough empty bands, otherwise return 0
        if (
            np.where(eigenvalues[i, :] > ev2ha(vbm + direct_gap + delta_energy))[0].size
            == 0
        ):
            return 0, 0
        else:
            max_band_idx.append(
                np.where(eigenvalues[i, :] > ev2ha(vbm + direct_gap + delta_energy))[0][
                    0
                ]
            )
    max_band_idx = np.max(max_band_idx)

    # associated frequency range
    correct_gap = direct_gap + scissor
    wrange = [
        np.max([0, correct_gap - delta_energy]),
        correct_gap + delta_energy,
    ]

    return wrange, max_band_idx


def get_eps_axis(id, base_dir):
    """
    Evaluates for which cartesian directions the dielectric tensor has to be calculated.
    INPUT:
        id:         ID of the material
        base_dir:   path to the base directory
    OUTPUT:
        row_idx:    array of independet directions of the dielectric tensor
    """

    # load the class
    _, _, ibrav = qe_helper.qe_init_structure(id, base_dir)
    # the ibrav variable is weird. Low-sysmmetry systems have to be handled differently
    if ibrav == 0:
        pass
    else:
        ibrav = ibrav["ibrav"]

    # go over all symmetries
    # (we ignore offdiagonal elements of the dielectric tensor)
    if ibrav in [1, 2, 3, -3]:
        row_idx = [0]  # x
    elif ibrav in [4, 5, -5, 6, 7]:
        row_idx = [0, 2]  # x, z
    else:
        row_idx = [0, 1, 2]  # x, y, z

    return row_idx
