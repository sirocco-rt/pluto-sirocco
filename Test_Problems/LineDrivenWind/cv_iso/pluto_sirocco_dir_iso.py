#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Script: pluto_sirocco_dir_iso.py
# Author: Nick Higginbottom
# Revised By: Amin Mosallanezhad (a.mosallanezhad@soton.ac.uk)
# Last Modified: 2024-11-19
# -----------------------------------------------------------------------------

import subprocess
import sys
import glob
import os
import shutil
from astropy import constants as c
from astropy import units as u
from scipy.integrate import quad
import pyPLUTO.pload as pp
import numpy as np
import pluto_sirocco_sub as pss  # Importing custom module for PLUTO-Python subroutines
import re
from pluto_sirocco_config import data  # Import configuration data


# Attempt to import code units from definitions.h using the get_units function from pluto_python_sub
# This is necessary because PLUTO expects geometrical parameters in code units,
# and the translation between code and physical units is set by the UNIT_ commands.
try:
    UNIT_DENSITY, UNIT_LENGTH, UNIT_VELOCITY = pss.get_units()
except:
    print("Unable to open definitions.h file - big problems")
    sys.exit()

# Set up initial parameters for the simulation run
t0 = 1.0  # Initial run time for the PLUTO simulation (in seconds)
dt = 2.0  # Time increment between calls to PLUTO (in seconds)

# Ensure t0 is at least 1.0
if t0 == 0.0:
    print("We need to run for at least one second, dummy")
    t0 = 1.0

init_sir_cycles = 5   # Initial number of cycles that SIROCCO will perform
sir_cycles = 2        # Number of cycles in SIROCCO once started

# Open a log file to record the run progress
logfile = open("run_logfile", "w")

# Handle restarts by checking for command-line arguments
if len(sys.argv) < 2:
    # No restart argument provided; start from the beginning
    istart = 0
    sim_time = t0
else:
    # Attempt to parse the restart cycle number from the command-line argument
    try:
        istart = int(sys.argv[1])
    except ValueError:
        print(f"Invalid argument: {sys.argv[1]} is not an integer.")
        sys.exit(1)

    if istart > 0:
        # Restarting from a previous cycle
        root = f"{istart:08d}"
        directory = "cycle" + root
        print(f"We will be trying to restart using files in {directory}")

        # Copy necessary files from the previous cycle's directory
        try:
            for filename in os.listdir(directory):
                source_file = os.path.join(directory, filename)
                if os.path.isfile(source_file):
                    shutil.copy(source_file, '.')
        except Exception as e:
            print(f"Cannot restart: {e}")
            sys.exit(1)

        # Update simulation time based on the last PLUTO output
        print("Last run finished at ", pp.pload(istart).SimTime)
        sim_time = pp.pload(istart).SimTime + dt
    else:
        # Starting from the first cycle
        istart = 0
        sim_time = t0

# In the absence of force multiplier files, use a k-alpha formulation for the initial PLUTO run
data["k"] = 0.59
data["alpha"] = -0.6  # Alpha is normally negative

# Copy initial flux files needed for the simulation
try:
    shutil.copy("directional_flux_theta_init.dat", "directional_flux_theta.dat")
    shutil.copy("directional_flux_phi_init.dat", "directional_flux_phi.dat")
    shutil.copy("directional_flux_r_init.dat", "directional_flux_r.dat")
except Exception as e:
    print(f"Error copying directional flux files: {e}")
    sys.exit(1)

# Main simulation loop
for i in range(istart, 750):
    # For cycles beyond the initial run, set k and alpha to dummy values
    if i > 0:
        data["k"] = 999
        data["alpha"] = 999

    # Format cycle number for consistent file naming
    root = f"{i:08d}"
    logfile.write(f"Making a PLUTO input file for cycle {i}\n")
    print(f"Making a PLUTO input file for cycle {i}\n")

    # Create a PLUTO input file using the current simulation time and data
    pss.pluto_input_file(sim_time, data)

    # Determine the command to run PLUTO, considering the number of processors
    if i == 0:
        # First cycle; no restart
        if data["nproc_pluto"] == 1:
            cmd = ["./pluto"]
        else:
            cmd = ["mpirun", "-np", str(data["nproc_pluto"]), "./pluto"]
    else:
        # Subsequent cycles; restart from the previous cycle
        if data["nproc_pluto"] == 1:
            cmd = ["./pluto", "-restart", str(i)]
        else:
            cmd = ["mpirun", "-np", str(data["nproc_pluto"]), "./pluto", "-restart", str(i)]

    logfile.write("Running PLUTO run\n")
    print("Running PLUTO run\n")

    logfile.write("Command line: " + ' '.join(cmd) + "\n")
    print("Command line: " + ' '.join(cmd) + "\n")

    # Execute the PLUTO simulation
    try:
        with open('pluto_log', 'w') as log_file:
            result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            print("Error: PLUTO command failed. Check pluto_log for details.")
            sys.exit(1)
    except Exception as e:
        print(f"Failed to run PLUTO: {e}")
        sys.exit(1)

    logfile.write("Finished PLUTO run\n")
    print("Finished PLUTO run\n")

    # Update simulation time and cycle number for the next iteration
    ifile = i + 1
    sim_time = sim_time + dt
    root = f"{ifile:08d}"
    dbl = f"data.{ifile:04d}.dbl"
    directory = "cycle" + root

    logfile.write(f"Turning dbl file {dbl} into a SIROCCO model file\n")
    print(f"Turning dbl file {dbl} into a SIROCCO model file\n")

    # Convert PLUTO output to a SIROCCO-readable model file
    sir_model_file = pss.pluto2sir_rtheta(ifile)
    logfile.write(f"Made a SIROCCO model file called {sir_model_file}\n")
    print(f"Made a SIROCCO model file called {sir_model_file}\n")

    # Create a SIROCCO parameter file for the simulation
    logfile.write("Making a SIROCCO input file\n")
    print("Making a SIROCCO input file\n")
    pss.sirocco_input_file(root, data, cycles=init_sir_cycles + i * sir_cycles)
    logfile.write("Successfully made a SIROCCO input file\n")
    print("Successfully made a SIROCCO input file\n")

    # Copy the SIROCCO parameter file to a generic name (input.pf) to maintain consistency
    try:
        shutil.copy(f"{root}.pf", "input.pf")
    except Exception as e:
        print(f"Error copying {root}.pf to input.pf: {e}")
        sys.exit(1)

    logfile.write(f"Command line: cp {root}.pf input.pf\n")
    print(f"Command line: cp {root}.pf input.pf\n")

    # Running the SIROCCO simulation
    if i > 0:
        # For cycles beyond the initial run, restart SIROCCO using the previous windsave file
        logfile.write("Restarting SIROCCO\n")
        print("Restarting SIROCCO\n")

        # Update densities in the windsave file using the new model
        logfile.write("Running modify_wind on the old windsave\n")
        print("Running modify_wind on the old windsave\n")
        cmd = ["modify_wind" + data["RAD_CODE_VER"], "-model_file", sir_model_file, "input"]
        logfile.write("Command line: " + ' '.join(cmd) + "\n")
        print("Command line: " + ' '.join(cmd) + "\n")

        # Execute modify_wind to update the windsave file
        try:
            with open('mod_wind_log', 'w') as log_file:
                result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
            if result.returncode != 0:
                print("Error: modify_wind command failed. Check mod_wind_log for details.")
                sys.exit(1)
        except Exception as e:
            print(f"Failed to run modify_wind: {e}")
            sys.exit(1)

        # Copy the new windsave file to input.wind_save
        try:
            shutil.copy("new.wind_save", "input.wind_save")
        except Exception as e:
            print(f"Error copying new.wind_save to input.wind_save: {e}")
            sys.exit(1)

        logfile.write("Copied new.wind_save to input.wind_save\n")
        print("Copied new.wind_save to input.wind_save\n")

        # Prepare the command to run Python in restart mode
        cmd = ["mpirun", "-np", str(data["nproc_sir"]), data["RAD_CODE"] + data["RAD_CODE_VER"], "-f", "-r", "input.pf"]
    else:
        # First run of Python; no restart
        cmd = ["mpirun", "-np", str(data["nproc_sir"]), data["RAD_CODE"] + data["RAD_CODE_VER"], "-f", "input.pf"]

    logfile.write("Running SIROCCO\n")
    logfile.write("Command line: " + ' '.join(cmd) + "\n")
    print("Running SIROCCO\n")
    print("Command line: " + ' '.join(cmd) + "\n")

    # Execute the Python simulation
    try:
        with open('sirocco_log', 'w') as log_file:
            result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            print("Error: SIROCCO command failed. Check sirocco_log for details.")
            sys.exit(1)
    except Exception as e:
        print(f"Failed to run SIROCCO: {e}")
        sys.exit(1)

    logfile.write("Finished SIROCCO\n")
    print("Finished SIROCCO\n")

    # Process SIROCCO outputs to generate files needed by PLUTO
    cmd = ["mpirun", "-np", str(data["nproc_sir"]), "rad_hydro_files" + data["RAD_CODE_VER"], "input"]


    logfile.write("Command line: " + ' '.join(cmd) + "\n")
    print("Command line: " + ' '.join(cmd) + "\n")
    try:
        with open('rad_hydro_files_output', 'w') as log_file:
            result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            print("Error: rad_hydro_files command failed. Check rad_hydro_files_output for details.")
            sys.exit(1)
    except Exception as e:
        print(f"Failed to run rad_hydro_files: {e}")
        sys.exit(1)

    # Copy the resulting Python data files with the current root name
    try:
        shutil.copy("py_heatcool.dat", f"{root}_py_heatcool.dat")
        shutil.copy("py_driving.dat", f"{root}_py_driving.dat")
        shutil.copy("py_pcon_data.dat", f"{root}_py_pcon_data.dat")
    except Exception as e:
        print(f"Error copying py_* files: {e}")
        sys.exit(1)

    # # Calculate accelerations using the output from Python
    # pss.accel_calc(ifile)

    # Run the external CAK (Castor-Abbott-Klein) code for radiation force calculations
    cmd = ["mpirun", "-np", str(data["nproc_cak"]), "./cak_v3"]
    print("Running CAK")
    print("Command line: " + ' '.join(cmd))
    try:
        with open('cak_output', 'w') as log_file:
            result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            print("Error: cak_v3 command failed. Check cak_output for details.")
            sys.exit(1)
    except Exception as e:
        print(f"Failed to run cak_v3: {e}")
        sys.exit(1)
    print("Finished CAK\n")

    # Create prefactors file for the Blondin heating and cooling rates
    print("Creating prefactors file for the Blondin heating and cooling rates")
    pss.pre_calc(ifile)
    print("Finished creating prefactors file\n")

    # Clean up old directories to manage storage space
    old_directory = f"cycle{ifile - 5:08d}"
    if os.path.exists(old_directory):
        try:
            shutil.rmtree(old_directory)
        except Exception as e:
            print(f"Error removing old directory {old_directory}: {e}")

    # Create a new directory for the current cycle to store output files
    try:
        os.mkdir(directory)
    except FileExistsError:
        # If the directory exists, rename it and create a new one
        try:
            shutil.rmtree(directory + "_old")
        except FileNotFoundError:
            pass
        shutil.move(directory, directory + "_old")
        os.mkdir(directory)

    # Copy all relevant files into the cycle directory for storage
    try:
        shutil.copy("dbl.out", directory)
        shutil.copy("pluto.ini", directory)
        shutil.copy("restart.out", directory)
        shutil.copy("grid.out", directory)
        shutil.copy("definitions.h", directory)
        shutil.copy(dbl, directory)
        # Copy all Python data files matching the pattern py_*.dat
        for file in glob.glob("py_*.dat"):
            shutil.copy(file, directory)
        shutil.copy("M_UV_data.dat", directory)
        shutil.copy("prefactors.dat", directory)
        shutil.copy("directional_flux_r.dat", directory)
        shutil.copy("directional_flux_theta.dat", directory)
        shutil.copy("directional_flux_phi.dat", directory)
        shutil.copy("input.wind_save", directory)
        # Move the model and parameter files into the cycle directory
        shutil.move(sir_model_file, directory)
        shutil.move(f"{root}.pf", directory)
        shutil.move(f"{root}_py_heatcool.dat", directory)
        shutil.move(f"{root}_py_driving.dat", directory)
        shutil.move(f"{root}_py_pcon_data.dat", directory)
    except Exception as e:
        print(f"Error copying files to {directory}: {e}")
        sys.exit(1)

    logfile.write("Finished tidying up\n")
    print("Finished tidying up\n")

# Close the log file after the simulation completes
logfile.close()
print("Fin")
