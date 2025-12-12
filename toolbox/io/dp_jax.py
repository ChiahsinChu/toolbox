# SPDX-License-Identifier: LGPL-3.0-or-later
"""DeepMD JAX molecular dynamics module.

This module provides classes for running molecular dynamics simulations
using DeepMD JAX implementation with trajectory dumping functionality.
"""

import sys
import time
from typing import List, Optional

import numpy as np
from ase import Atoms, io
from deepmd_jax.md import Simulation as _Simulation


class TrajDump:
    """Class for dumping trajectory data during simulation.
    
    This class handles writing trajectory data at specified intervals
    during molecular dynamics simulations.
    """
    
    def __init__(
        self,
        atoms: Atoms,
        fname: str,
        interval: int,
        vel: bool = False,
        **kwargs,
    ) -> None:
        """Initialize TrajDump.
        
        Parameters
        ----------
        atoms : ase.Atoms
            Atoms object for the system
        fname : str
            Output filename
        interval : int
            Dump interval in timesteps
        vel : bool, optional
            Whether to dump velocities, by default False
        **kwargs
            Additional keyword arguments for writing
        """
        self.fname = fname
        self.interval = interval
        self.vel = vel
        self.atoms = atoms

        self.write_settings = kwargs

    def write(self, positions, cell):
        """Write current frame to file.
        
        Parameters
        ----------
        positions : array_like
            Atomic positions
        cell : array_like
            Unit cell vectors
        """
        self.atoms.set_positions(positions)
        self.atoms.set_cell(cell)
        io.write(
            self.fname,
            self.atoms,
            **self.write_settings,
        )


class Simulation(_Simulation):
    """DeepMD JAX molecular dynamics simulation class.
    
    This class extends the base DeepMD JAX Simulation class to provide
    trajectory dumping and logging functionality.
    
    Example
    -------
    # setup simulation
    sim = Simulation(
        model_path="model.pkl",  # Has to be an 'energy' or 'dplr' model
        box=box,  # Angstroms
        type_idx=type_idx,  # here the index-element map (e.g. 0-Oxygen, 1-Hydrogen) must match the dataset used to train the model
        mass=[15.9994, 1.0078, 195.08],  # Oxygen, Hydrogen
        routine="NVT",  # 'NVE', 'NVT', 'NPT' (Nosé-Hoover)
        dt=0.5,  # femtoseconds
        initial_position=initial_position,  # Angstroms
        temperature=330,  # Kelvin
        report_interval=1,  # Report every 100 steps
        seed=np.random.randint(1, 1e5),  # Random seed
    )

    sim.run(
        n_steps,
        [
            TrajDump(atoms, "pos_traj.xyz", 10, append=True),
            TrajDump(atoms, "vel_traj.xyz", 10, vel=True, append=True),
        ],
    )

    """

    def __init__(
        self,
        model_path,
        box,
        type_idx,
        mass,
        routine,
        dt,
        initial_position,
        log_file: Optional[str] = "deepmd_jax.stdout",
        **kwargs,
    ):
        """Initialize Simulation.
        
        Parameters
        ----------
        model_path : str
            Path to trained model file
        box : array_like
            Simulation box dimensions in Angstroms
        type_idx : array_like
            Index-element mapping for atom types
        mass : array_like
            Atomic masses
        routine : str
            Ensemble type ('NVE', 'NVT', 'NPT')
        dt : float
            Timestep in femtoseconds
        initial_position : array_like
            Initial atomic positions in Angstroms
        log_file : Optional[str], optional
            Log file path, by default "deepmd_jax.stdout"
        **kwargs
            Additional keyword arguments
        """
        super().__init__(
            model_path,
            box,
            type_idx,
            mass,
            routine,
            dt,
            initial_position,
            **kwargs,
        )
        self.log_file = log_file
        if log_file is not None:
            # export all stdout to log_file
            self._stdout = sys.stdout
            with open(log_file, "w", encoding="utf-8") as f:
                sys.stdout = f

    def __del__(self):
        """Restore stdout when object is deleted."""
        if self.log_file is not None:
            sys.stdout.close()
            sys.stdout = self._stdout

    def _initialize_run(self, steps):
        """Initialize simulation run.
        
        This function resets trajectory for each new run and
        initializes run variables. If the simulation has not
        been run before, it includes the initial state.
        
        Parameters
        ----------
        steps : int
            Total number of steps to run
        """
        print(f"# Running {steps} steps...")
        self._offset = self.step - int(self._is_initial_state)
        self._tic_of_this_run = time.time()
        self._tic_between_report = time.time()
        self._error_code = 0
        self._print_report()
        self._is_initial_state = False

    def run(self, steps, dump_list: List[TrajDump]):
        """Run the simulation for a number of steps."""
        self._initialize_run(steps)
        remaining_steps = steps
        while remaining_steps > 0:
            # run the simulation for a jit-compiled chunk of steps
            next_chunk = min(
                self.report_interval - self.step % self.report_interval,
                self._step_chunk_size,
                remaining_steps,
            )
            states = (
                self._state,
                self._typed_nbrs if self._static_args["use_neighbor_list"] else None,
                self._error_code,
                self._neighbor_update_profile,
            )
            states_new, traj = self._multiple_inner_step_fn(states, next_chunk)
            state_new, typed_nbrs_new, error_code, profile = states_new
            self._error_code |= error_code

            if self._error_code & 16:
                print("# Warning: Nan or Inf encountered in simulation. Terminating.")
                remaining_steps = 0

            # If there is any hard overflow, we have to re-run the chunk
            if not (
                self._error_code == 0 or self._error_code == 4 or self._error_code == 16
            ):
                self._resolve_error_code()
                continue

            # If nothing overflows, update the tracked state and record the trajectory
            self._keep_nbr_or_lattice_up_to_date()
            self._state = state_new
            self._typed_nbrs = typed_nbrs_new
            self._neighbor_update_profile = profile
            if "NPT" in self._routine:
                self._current_box = self._state.box
            pos_traj, vel_traj, box_traj = traj
            self.step += next_chunk
            remaining_steps -= next_chunk

            # Report at preset regular intervals
            if self.step % self.report_interval == 0 or remaining_steps == 0:
                self._print_report()

            for dump in dump_list:
                if self.step % dump.interval == 0 or remaining_steps == 0:
                    cell = np.concatenate([np.array(box_traj[-1]), [90, 90, 90]])
                    dump.write(pos_traj[-1] if not dump.vel else vel_traj[-1], cell)

        self._print_run_profile(steps, time.time() - self._tic_of_this_run)
        self._keep_nbr_or_lattice_up_to_date()
