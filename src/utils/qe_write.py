"""
Functions that write input files for QE (but don't actually run them)
"""

# external imports
import numpy as np

# local imports
import src.utils.basic_utils as basic_utils
import src.utils.qe_helper as qe_helper


def write_pw(
    calc_data,
    occupations="smearing",
    disk_io="medium",
    diagonalization="david",
    symmorphic=False,
    kpoints_mode="automatic",
    nosym=False,
    shift=[0, 0, 0],
    n_bands=0,
):
    """
    Writes the input file for a pw.x calculation
    INPUT:
        calc_data:      calc_data of the calculation which should be run
        For the other inputs, see the pw.x documentation
    OUTPUT:
        output_filename:       Name of the written input file
    """

    output_filename = (
        "pw_"
        + calc_data.identifier
        + "_"
        + calc_data.calc_type
        + "_k"
        + str(calc_data.kppa)
        + "_E"
        + str(calc_data.pw_cutoff)
    )
    pseudo = {}
    for elem in calc_data.structure.types_of_species:
        pseudo[elem.name] = elem.name + ".upf"
    k_points_grid = calc_data.k_points_grid
    control = {
        "calculation": calc_data.calc_type,
        "prefix": calc_data.identifier,
        "outdir": f"out/",
        "pseudo_dir": calc_data.pseudo,
        "disk_io": disk_io,
    }
    if n_bands == 0:
        system = {
            "ecutwfc": calc_data.pw_cutoff,
            "occupations": occupations,
            "degauss": basic_utils.ev2ha(0.025) / 2,  # 25meV smearing to "emulate" 300K
            "smearing": "gaussian",
            "nbnd": int(
                np.ceil(
                    np.max(
                        [
                            qe_helper.qe_get_electrons(calc_data),
                            8,
                        ]
                    )
                )
            ),
            "force_symmorphic": symmorphic,
            "nosym": nosym,
        }
    else:
        system = {
            "ecutwfc": calc_data.pw_cutoff,
            "occupations": occupations,
            "degauss": basic_utils.ev2ha(0.025) / 2,  # 25meV smearing to "emulate" 300K
            "smearing": "gaussian",
            "nbnd": n_bands,
            "force_symmorphic": symmorphic,
            "nosym": nosym,
        }
    if calc_data.ibrav != 0:
        system.update(calc_data.ibrav)
    electrons = {
        "conv_thr": 1e-10,
        "diagonalization": diagonalization,
        "diago_full_acc": True,
        "diago_thr_init": 5e-6,
    }
    input = qe_helper.qe_PWInput(
        calc_data.structure,
        pseudo=pseudo,
        kpoints_grid=k_points_grid,
        control=control,
        system=system,
        electrons=electrons,
        kpoints_mode=kpoints_mode,
        kpoints_shift=shift,
    )
    input.write_file(output_filename + ".in")
    return output_filename


def write_nscf_yambo(
    calc_data,
    n_bands,
    occupations="smearing",
    kpoints_shift=[0, 0, 0],
    diagonalization="david",
):
    """
    Writes the input file for a QuantumEspresso nscf calculation.
    calc_data is the corresponding calc_data_class object.
    n_bands should be set high enough for optics and gw calculations
    occupations can be set to tetrahedron for dos style calculation.
    kpoints_shift is used for the fine grid in the double grid method for yambo.
    Returns the path to the input file (without the ending).
    """

    output_filename = (
        "pw_"
        + calc_data.identifier
        + "_"
        + calc_data.calc_type
        + "_k"
        + str(calc_data.kppa)
        + "_E"
        + str(calc_data.pw_cutoff)
    )
    pseudo = {}
    for elem in calc_data.structure.types_of_species:
        pseudo[elem.name] = elem.name + ".upf"
    k_points_grid = calc_data.k_points_grid
    control = {
        "calculation": calc_data.calc_type,
        "prefix": calc_data.identifier,
        "outdir": f"out/",
        "pseudo_dir": calc_data.pseudo,
        "disk_io": "low",
    }
    system = {
        "ecutwfc": calc_data.pw_cutoff,
        "occupations": occupations,
        "degauss": basic_utils.ev2ha(0.025) / 2,  # 25meV smearing to "emulate" 300K
        "smearing": "gaussian",
        "nbnd": n_bands,
        "force_symmorphic": True,
    }
    if calc_data.ibrav != 0:
        system.update(calc_data.ibrav)
    electrons = {
        "conv_thr": 1e-10,
        "diago_full_acc": True,
        "diago_thr_init": 5e-6,
        "diagonalization": diagonalization,
    }
    input = qe_helper.qe_PWInput(
        calc_data.structure,
        pseudo=pseudo,
        kpoints_grid=k_points_grid,
        kpoints_shift=kpoints_shift,
        control=control,
        system=system,
        electrons=electrons,
    )
    input.write_file(output_filename + ".in")
    return output_filename
