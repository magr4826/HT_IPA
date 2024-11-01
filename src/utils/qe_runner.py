"""
Functions that can start pw.x calculations
"""

# external imports
import os
import pickle
import shutil

# local imports
from src.utils.calc_data_class import calc_data
import src.utils.qe_write as qe_write
import src.utils.qe_helper as qe_helper


def qe_pw_run(
    calc_data,
    input_name,
    ncores,
    kwargs={},
):
    """
    This function writes and starts a pw.x calculation and does some rudimentary error handling.
    INPUT:
        calc_data:      The calculation that should be run
        input_name:     Which function should be used to write the input (usually a function in qe_write)
        ncores:         How many cores should be used by the calculation
        kwargs:         Keyword arguments to be passed to the function input_name
    OUTPUT:
        filename:       name of the input file which was executed
    """

    # write the input file for pw.x and start the calculation
    filename = input_name(calc_data, **kwargs)
    os.system(f"mpirun -np {ncores} pw.x -inp {filename}.in > {filename}.out")

    # if there is an error, try to use the more stable paro diagonalization
    # (errors can happen with nscf calculations, we never observed one in scf calculations)
    with open(filename + ".out") as f:
        if "Error" in f.read():
            kwargs["diagonalization"] = "paro"
            filename = input_name(calc_data, **kwargs)
            os.system(f"mpirun -np {ncores} pw.x -inp {filename}.in > {filename}.out")
            open("paro.txt", "a").close()  # create a file so we know that paro is used

    # check if the "eigenvalues not converged" appear
    # but we need to separate two cases
    # 1st case: last iteration of a scf calculation
    # 2nd case: nscf calculation
    with open(filename + ".out") as f:
        # parse the output file
        if "nscf" in filename:
            out_str = f.read()  # keep the whole output file
        else:
            out_str = f.read()  # only we the part after the last iteration
            out_str = out_str.split("iteration #")[
                -1
            ]  # we only keep the part after the last iteration

        # check if the eigenvalues did not converge
        if "eigenvalues not converged" in out_str:
            kwargs["diagonalization"] = "paro"
            filename = input_name(calc_data, **kwargs)
            os.system(f"mpirun -np {ncores} pw.x -inp {filename}.in > {filename}.out")
            open("paro.txt", "a").close()  # create a file so we know that paro is used

    return filename


def qe_convergence_checker(id, wf_dir, ncores, conv_dir="qe_convergence_SG15"):
    """
    Check if convergence for a material has been performed
    If yes:     Copy output of converged scf calculation to workflow dir
                Return converged convergence parameters
    If no:      Perform one scf with default parameters in workflow dir
                Return default (unconverged) convergence parameters
    This function should be called at the start of most QE workflows
    INPUT:
        id:             Materials Project id of the material
        wf_dir:         Directory of the workflow which called this function
        ncores:         Number of cores used (in case convergence hasnt been run before)
        conv_dir:       Name of the workflow directory in which the convergence was run (in case there are various versions)
    OUTPUT:
        conv_flag:      Boolean, whether convergence ran before or not
        calc_data_conv: calc_data of the converged calculation, if it ran before, or the default parameters otherwise
    """

    # check if convergence has been succesfully carried out
    conv_flag = False
    conv_dir = os.path.join(wf_dir, os.pardir, conv_dir)
    if os.path.exists(conv_dir):
        if os.path.exists(os.path.join(conv_dir, "conv_calc_data.pckl")):
            conv_flag = True
            print(
                f"\nConvergence has been carried out, using converged parameters and results from\n{conv_dir}",
                flush=True,
            )
        else:
            print(
                f"\nConvergence has started but not successfully finished, starting calculation with initial values for\n{conv_dir}",
                flush=True,
            )
    else:
        print(
            f"\nConvergence has not been started. Starting calculation with initial values for\n{conv_dir}",
            flush=True,
        )

    if conv_flag:
        # load the converged calc-data object
        f = open(os.path.join(conv_dir, "conv_calc_data.pckl"), "rb")
        calc_data_conv = pickle.load(f)
        f.close()

        # Create the necessary folders
        if not os.path.exists(os.path.join(wf_dir, "out")):
            os.mkdir(os.path.join(wf_dir, "out"))
        if not os.path.exists(os.path.join(wf_dir, "out", id + ".save")):
            os.mkdir(os.path.join(wf_dir, "out", id + ".save"))

        # copy the results of the scf-calculation
        shutil.copy(
            os.path.join(conv_dir, "out", id + ".xml"),
            os.path.join(wf_dir, "out", id + ".xml"),
        )  # not sure if this is even necessary
        shutil.copy(
            os.path.join(conv_dir, "out", id + ".save", "data-file-schema.xml"),
            os.path.join(wf_dir, "out", id + ".save", "data-file-schema.xml"),
        )
        shutil.copy(
            os.path.join(conv_dir, "out", id + ".save", "charge-density.dat"),
            os.path.join(wf_dir, "out", id + ".save", "charge-density.dat"),
        )
        if "PAW" in conv_dir:
            shutil.copy(
                os.path.join(conv_dir, "out", id + ".save", "paw.txt"),
                os.path.join(wf_dir, "out", id + ".save", "paw.txt"),
            )
        return conv_flag, calc_data_conv

    else:
        base_dir = os.path.join(wf_dir, os.pardir, os.pardir, os.pardir)
        structure, name, ibrav = qe_helper.qe_init_structure(id, base_dir)
        calc_data_conv = calc_data(
            structure,
            name,
            id=id,
            ibrav=ibrav,
            calc_type="scf",
            kppa=1500,
            pseudo=os.path.join(base_dir, "pseudo", "SG15"),
        )

        qe_pw_run(calc_data_conv, qe_write.write_pw, ncores, kwargs={})
        return conv_flag, calc_data_conv
