Molecular dynamics parameters
===============================

General information
-------------------

The entry script for ChemLab is in `src/start_simulation.py` which provides a number of command-line options.
The complete list can be found here: :doc:`start_simulation`.

All options described here can be either passed by `--<option_name>` in the command-line or placed in the `params` file.
The parameters from `params` file can be loaded by:

.. code-block:: sh

   $ python src/start_simulation.py @params

General options
^^^^^^^^^^^^^^^

.. glossary::

      conf
         Input coordinate file. Currently only `.gro` file format is supported.

      topol
         Input topology file. This is a GROMACS-like file format. See: :doc:`topology`.

      node_grid
         If provided then use a custom node grid, format: `node_x,node_y,node_z`.

      skin
         The parameter for Verlet list algorithm.

      output_prefix
         Add prefix to all output files.

      trj_collect
         Collect trajectory every `n`-steps.

      energy_collect
         Collect system information like potential energy, temperature and others every `n`-steps.

      topol_collect
         Collect topology information every `n`-steps. The information are stored in HDF5 file along with trajectory
         information.

      reactions
         The reaction settings file.

      exclusion_list
         The file with the list of particle pairs to be excluded from non-bonded interactions.

      max_force
         The maximum value of the force in the system. If the force computed on any particle is larger than this value
         then the `max_force` is used (preserving the direction of force vector).


Running
^^^^^^^

.. glossary::

      run
         The number of simulation steps to run

      start_ar
         The simulation step when the chemical reactions are enabled

      stop_ar
         The simulation step when the chemical reactions are disabled

      rng_seed
         The seed for random number generator.

      gen_velocity
         Generate velocities according to Maxwell-Boltzmann distribution