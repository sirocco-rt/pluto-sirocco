#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Script: pluto_sirocco_config.py
# Purpose: Define physical and numerical parameters for the simulation pipeline.
# Notes: Update the parameter values as per the specific simulation requirements.
# Author: [Amin Mosallanezhad, a.mosallanezhad@soton.ac.uk]
# -----------------------------------------------------------------------------

# Import necessary modules
from astropy import constants as c
import numpy as np

# Define the data dictionary with the initial setup

# Number of processors
data = {
    "nproc_sir": 50,           # Number of cores for Python code
    "nproc_cak": 50,           # Number of cores for CAK code
    "CAK_EXEC": "cak_v4",   # Executable version/name; copy it into the run folder
    "nproc_pluto": 50,         # Number of cores for Pluto code
    "nproc_radhydro": 50,      # Number of cores for rad_hydro_files code
    "nproc_spec": 50,          # Number of cores for the spectral-only runs

    # Radial grid parameters
    "R_MIN": 0.87,              # Radial grid minimum (units depend on context)
    "R_MAX": 8.7,              # Radial grid maximum
    "N_R": 128,                # Number of radial bins
    "s_r": 1.05,               # Radial scaling factor

    # Theta grid parameters
    "T_MIN": np.radians(0.5),   # Theta grid minimum in radians (polar cutout)
    "T_MAX": np.radians(90.0),  # Theta grid maximum in radians
    "N_T": 96,                  # Number of theta bins
    "s_theta": 0.95,            # Theta scaling factor

    # Temperature and radiation parameters
    "T_X": 160000.0,            # Temperature in Kelvin (unused in current runs)
    "F_X": 0.1,                 # X-ray flux proportion
    "L_STAR": 9.05e34,          # Luminosity of the central source in erg/s
    "rad_force": 1,             # Radiation force enabled (1=True, 0=False)

    # Density parameters
    "RHO_0": 1e-9,              # Midplane density in g/cm^3
    "RHO_ALPHA": 0.0,           # Midplane density drop-off (unused)
    "R_0": 8.7e8,               # Density scaling radius in cm

    # Physical constants
    "MU": 0.6,                  # Mean molecular mass (atomic mass units)
    "T_ISO": 40000.0,           # Isothermal temperature in Kelvin
    "DFLOOR": 1e-22,            # Density floor in g/cm^3 
    "GAMMA": 1.66667,           # Adiabatic index

    # System properties
    "system_type": 'star',      # Geometry type
    "CENT_RADIUS": 8.7e8,       # Radius of the central object in cm
    "CENT_MASS": 0.6 * c.M_sun.cgs.value,  # Mass of the central object in g

    # Central source and boundary layer
    "cent_spectype": "none",    # Central Sirocco source disabled
    "boundary_layer": "no",     # Boundary layer disabled
    "L_BL": 9.05e34,            # Boundary layer luminosity in erg/s
    "T_BL": 40000.0,            # Boundary layer temperature in Kelvin
    "NPHOT": 1e7,               # Number of photons for Sirocco

    # Disk properties
    "disk_radiation": "yes",    # Disk radiation enabled
    "SIR_DISK_MDOT": 3.14e-8,    # Disk accretion rate in solar masses per year
    "DISK_MDOT": 3.14e-8 * c.M_sun.cgs.value / (365.25 * 24 * 3600),  # Disk accretion rate in g/s
    "DISK_TRUNC_RAD": 8.7e9,    # Maximum radius of the disk in cm

    # Wind and line transfer properties
    # Volume filling factor used when building every SIROCCO wind.
    # Spectrum-only runs inherit it from wind_save: rebuild/reconverge the
    # wind after changing this value before calculating new spectra.
    "FILLING_FACTOR": 0.01,       # 1=smooth; 0 < f < 1=microclumped
    "wind_radiation": "yes",     # Wind radiation enabled
    "line_trans": "macro",       # Line transfer mode
    "RAD_CODE": "sirocco",
    "RAD_CODE_VER":"-1.3",

    # -------------------------------------------------------------------------
    # PER-CYCLE IONISATION-STATE ARCHIVE
    #
    # Independent of the SPEC_RUN block below: this does NOT depend on
    # whether spectral-only test runs are switched on, and does NOT only
    # apply to the handful of cycles selected by SPEC_START/SPEC_END/
    # SPEC_STEP. It runs (or not) for every single main-loop cycle.
    #
    # Normally, each main-loop cycle's full archive folder (cycleNNNNNNNN/,
    # including its wind_save) is deleted wholesale once it is more than
    # SPEC_KEEP-style N cycles old (currently hard-coded to 5 cycles in
    # pluto_sirocco_dir_ideal.py), since the wind_save binaries are large.
    # That means the ionisation state at any earlier cycle is normally lost
    # for good.
    #
    # "yes" -> after each main-loop cycle's normal archive step, run
    #          windsave2table on that cycle's wind_save, then copy
    #          everything just archived for that cycle (the windsave2table
    #          ion tables plus the same static/sirocco_* files already kept in
    #          cycleNNNNNNNN/) into ION_TABLE_DIR/cycleNNNNNNNN/ --
    #          EXCLUDING the wind_save binary, data.NNNN.dbl and grid.out.
    #          The .dbl is PLUTO's own per-cycle snapshot and is never
    #          pruned from the main run directory, and grid.out is
    #          identical, unchanging geometry metadata -- copying either
    #          into ions_data/ as well would just duplicate data already
    #          available elsewhere; a later plotting script can read them
    #          straight from the main run directory by cycle number.
    #          Folders under ION_TABLE_DIR are never pruned; they
    #          accumulate for the whole run, giving a lightweight,
    #          permanent, per-cycle record of the ionisation state for
    #          later analysis.
    # "no"  -> skip this step entirely; behaviour is unchanged from before.
    "ION_TABLE_ALL_CYCLES": "yes",
    "ION_TABLE_DIR":        "ions_data",

    # -------------------------------------------------------------------------
    # SPECTRAL-ONLY TEST RUNS  s
    # -------------------------------------------------------------------------

    "SPEC_RUN":        "no",       # "yes" enables; "no" leaves the pipeline untouched
    "SPEC_MODE":       "run",       # "run" = build folder AND execute inline (Iridis 6)
                                    # "prepare" = build folder only (inspect the .pf first)
    "SPEC_DIR":        "spec_data",

    # Cycle selection: SPEC_START, SPEC_START+SPEC_STEP, ... up to SPEC_END

    "SPEC_START":      300,
    "SPEC_END":        1250,
    "SPEC_STEP":       5,

    # Spectral run parameters
    "SPEC_CYCLES":     30,          # Spectrum_cycles (Ionization_cycles is always 0)
    "SPEC_NPHOT":      1e7,         # Photons_per_cycle for the spectral run
    "SPEC_EXTRA_PROC": "none",      # Wind_heating.extra_processes for the spectral run
    "SPEC_ATOMIC_DATA": "auto",
    # Spectrum block
    "SPEC_NWAVE":      10000,
    "SPEC_WAVEMIN":    1450,         # Angstroms (full UV: below Lyman limit
                                    # 912 through He II 1640 / C III] 1909 region)
    "SPEC_WAVEMAX":    1650,
    "SPEC_ANGLES":     [5, 10, 15, 20, 25, 30, 35, 40, 45,
                        50, 55, 60, 65, 70, 75, 80, 85],

    # wind_save post-processing / disk-usage management for SPEC_DIR.

    "SPEC_WINDSAVE2TABLE": "yes",   # "yes"/"no" -- run windsave2table per cycle

    # How many of the most-recently-selected spectral cycles keep their
    # *.wind_save binaries in SPEC_DIR:
    #   "all"   -> never delete any wind_save (disk-hungry, full fidelity)
    #   N (int) -> keep only the N most recent (default 5); older cycles keep
    #              their .spec, .pf, log and windsave2table tables, only the
    #              wind_save binaries are removed.
    "SPEC_KEEP_WINDSAVE": 5,
}

# End of script
