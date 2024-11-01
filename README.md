# High Throughput IPA (HT_IPA)
Code to run High Throughput IPA calculations on the Alexandria database [1]

**REQUIREMENTS:**
All version numbers are the versions that were tested - other version might work, but success is not guaranteed.

- Python > v3.9
    Installed with miniconda and the following packages (plus their dependencies):
        numpy = v1.26.4
        pandas = v2.2.1
        scipy = v1.13.0
        matplotlib = v3.8.4
        pymatgen = v2024.5.1
        ase = v3.22.1
- QuantumEspresso = v7.1 ("/bin" directory loaded into the path, compiled for MPI)
- Yambo = v5.2 ("/bin" directory loaded into the path, compiled for MPI)
- Wannier90 = v3.1 ("/bin" and "/utility" directory loaded into the path, compiled for MPI, used only for the kmesh.pl utility)


**USAGE:**

1. Convert the Alexandria database into .pckl-files using the conversion.ipynb notebook and place them in the input_data folder
2. Create an empty folder named "database" and a folder named "CONTROL", which contains the file "NCORES". In this file, specify the total number of cores you want to use (e.g., 640)
3. Specify the name of the .pckl-file and the workflows to run in the HT_watcher.py script
4. If necessary, adapt the HT_watcher.py script and the start_calc function in basic_utils.py to your cluster
5. Run the HT_watcher.py script
6. If necessary, restart failed calculations using HT_restart.py or HT_full_restart.py


**DATABASE:**
The main results of the calculations are saved in .json-files the database folder.
The database follows the structure of the original database.
We add the following data and parameters:

data:
    ipa_epsI_0:         Imaginary part of the xx-component of the dielectric tensor
    ipa_epsR_0:         Real part of the xx-component of the dielectric tensor
    ipa_indirect_gap:   Indirect gap determined on the converged k-point grid for the IPA calculation
    ipa_direct_gap:     Direct gap determined on the converged k-point grid for the IPA calculation
    ipa_eps_similarity: Similarity coefficient between ipa_epsI_0 using the final and the penultimate k-point density
    ipa_epsI_1:         Imaginary part of the yy-component of the dielectric tensor
    ipa_epsR_1:         Real part of the yy-component of the dielectric tensor
    ipa_epsI_2:         Imaginary part of the zz-component of the dielectric tensor
    ipa_epsR_2:         Real part of the zz-component of the dielectric tensor
The latter four properties are only present if they are different from eps_I_0 and eps_R_0 based on symmetry.
The gaps differ from the original gaps in the Alexandria database mostly due to the use of different codes, methods, and pseudopotentials.
params:
    pw_conv_k:          Converged K-point density in inverse Angstrom for the ground state calculation
    pw_conv_cutoff:     Converged cutoff in Rydberg for the ground state calculation
    ipa_eps_kppa:       Converged K-point density in inverse Angstrom for the IPA calculation
    ipa_nbands:         Number of bands used for the IPA calculation
    ipa_broad:          Value of the broadening in eV used for the IPA calculation


[1] https://alexandria.icams.rub.de/