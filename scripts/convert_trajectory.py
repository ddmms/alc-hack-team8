#!/Users/maciej.bartkowiak/python/builder/bin/python3.14

import os
os.environ.update(
    OMP_NUM_THREADS = '1',
    OPENBLAS_NUM_THREADS = '1',
    MKL_NUM_THREADS = '1',
    VECLIB_MAXIMUM_THREADS = '1',
    NUMEXPR_NUM_THREADS = '1'
)

########################################################
# This is an automatically generated MDANSE run script #
########################################################

from MDANSE.Framework.Converters.Converter import Converter

########################################################
# Job parameters                                       #
########################################################

parameters = {
    'atom_aliases': '{"": {"Ti": "Ti", "Sr": "Sr", "O": "O"}}',  # Atom mapping
    'fold': False,  # Fold coordinates into box
    'n_steps': '0',  # Number of time steps (0 for automatic detection)
    'output_files': ('../trajectories/MDANSE_trajectory', 32, (8, 128), 'gzip', 'no logs', 32768),  # MDANSE output trajectory. Values are: filename, data type, chunk size, compression, logfile output.
    'time_step': '50',  # Time step
    'time_unit': 'fs',  # Time step unit
    'trajectory_file': '../trajectories/small-sto-T50-nve.extxyz',  # An MD trajectory file supported by ASE
}

########################################################
# Setup and run the analysis                           #
########################################################

if __name__ == "__main__":
    ase = Converter.create("ASE")
    # Progress bars only available if tqdm available.
    # Install with `cli` optional dependency.
    ase.run(parameters, status=True, prog_bar=True)
