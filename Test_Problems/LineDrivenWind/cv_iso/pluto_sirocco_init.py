#!/usr/bin/env python

import argparse
import subprocess
import sys
import shutil
import os
from pluto_sirocco_config import data  # Import data from python_config
import pluto_sirocco_sub as pss

# -----------------------------------------------------------------------------
# Script: pluto_sirocco_init.py
# Purpose: Initialize and run the SIROCCO simulation pipeline with PLUTO.
# Notes: Sets nproc_sir dynamically from the config for mpirun.
# Author: Amin Mosallanezhad, a.mosallanezhad@soton.ac.uk
# -----------------------------------------------------------------------------


# Set simulation parameters
data["k"] = 0.59
data["alpha"] = -0.6  # Alpha is normally negative


init_sir_cycles = 30  # Initial cycles for the SIROCCO simulation
sim_time = 0.0
ifile = 0      # File index for data.0000.dbl
root = f"{ifile:08d}"  # Root identifier for file naming

# Logging setup
logfile_path = "processing_log.txt"
with open(logfile_path, "w") as logfile:
    logfile.write("Starting SIROCCO processing pipeline\n")

def log_and_print(message):
    """Log a message to both the terminal and logfile."""
    print(message)
    with open(logfile_path, "a") as logfile:
        logfile.write(message + "\n")

def cleanup():
    """Clean up temporary files in case of errors."""
    files_to_remove = ["input.pf", "sirocco_log", "rad_hydro_files_output"]
    for file in files_to_remove:
        if os.path.exists(file):
            os.remove(file)
    log_and_print("Cleaned up temporary files.")

# Resource monitoring
log_and_print(f"Using {data['nproc_sir']} processors.")

# Step 1: Generate PLUTO input file
log_and_print("Generating PLUTO input file.")
pss.pluto_input_file(sim_time, data)

# Step 2: Run PLUTO executable
log_and_print("Running PLUTO executable.")
try:
    subprocess.run(["./pluto"], check=False, stderr=subprocess.DEVNULL)
except Exception as e:
    log_and_print(f"Warning: Issue encountered with PLUTO - {e}")

# Step 3: Generate Python model file
dbl = f"data.{ifile:04d}.dbl"  # Path to the initial dbl file
log_and_print(f"Converting {dbl} into a Python model file.")
sirocco_model_file = pss.pluto2sir_rtheta(ifile)
log_and_print(f"Generated Python model file: {sirocco_model_file}")

# Step 4: Create SIROCCO parameter file (.pf)
log_and_print("Creating SIROCCO input parameter file (.pf).")
pss.sirocco_input_file(root, data, cycles=init_sir_cycles)
log_and_print("SIROCCO parameter file created successfully.")


# Step 5: Copy parameter file to input.pf
log_and_print("Copying parameter file to input.pf.")
source_file = f"{root}.pf"  # This should be '00000000.pf'
destination_file = "input.pf"

# Debug: Check if source_file exists
if os.path.exists(source_file):
    log_and_print(f"Source file {source_file} found. Proceeding to copy...")
    try:
        shutil.copy(source_file, destination_file)
        log_and_print(f"Successfully copied {source_file} to {destination_file}")
    except Exception as e:
        log_and_print(f"Error copying {source_file} to {destination_file}: {e}")
        cleanup()
        sys.exit(1)
else:
    log_and_print(f"Error: Source file {source_file} does not exist. Cannot copy to {destination_file}.")
    cleanup()
    sys.exit(1)


# Step 6: Run Python with mpirun and redirect output to python_log
# log_and_print("Running SIROCCO simulation with mpirun.")

try:
    cmd = [
        "mpirun",
        "-np",
        str(data["nproc_sir"]),
        "sirocco" + data["RAD_CODE_VER"],
        "-f",
        "-p",
        "2",
        "input.pf",
    ]

    log_and_print(f"Executing command: {' '.join(cmd)}")

    with open("sirocco_log", "w") as sirocco_log_file:
        subprocess.run(cmd, check=True, stdout=sirocco_log_file, stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    log_and_print(f"Error: SIROCCO simulation failed with exit code {e.returncode}")
    log_and_print(f"Standard Error: {e.stderr.decode()}")
    cleanup()
    sys.exit(1)

# Step 7: Run rad_hydro_files and capture output in rad_hydro_files_output
log_and_print("Running rad_hydro_files and saving output.")
try:
    cmd = ["mpirun", "-np", str(data["nproc_sir"]), "rad_hydro_files" + data["RAD_CODE_VER"], "input"]
    with open('rad_hydro_files_output', 'w') as output_file:
        subprocess.run(cmd, stdout=output_file, stderr=subprocess.STDOUT, check=True)
except subprocess.CalledProcessError as e:
    log_and_print("Error running rad_hydro_files. Check rad_hydro_files_output for details.")
    cleanup()
    sys.exit(1)

# Step 8: Rename flux files
log_and_print("Renaming flux files.")
try:
    shutil.move("directional_flux_r.dat", "directional_flux_r_init.dat")
    shutil.move("directional_flux_theta.dat", "directional_flux_theta_init.dat")
    shutil.move("directional_flux_phi.dat", "directional_flux_phi_init.dat")
except Exception as e:
    log_and_print(f"Error renaming flux files: {e}")
    cleanup()
    sys.exit(1)

log_and_print("Completed SIROCCO run and file management successfully.")
