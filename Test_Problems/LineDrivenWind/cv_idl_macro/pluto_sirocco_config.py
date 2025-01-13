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
    "nproc_sir": 40,           # Number of cores for Python code
    "nproc_cak": 40,           # Number of cores for CAK code
    "nproc_pluto": 40,         # Number of cores for Pluto code

    # Radial grid parameters
    "R_MIN": 0.87,             # Radial grid minimum (units depend on context)
    "R_MAX": 8.7,              # Radial grid maximum
    "N_R": 128,                # Number of radial bins
    "s_r": 1.05,               # Radial scaling factor

    # Theta grid parameters
    "T_MIN": np.radians(0.0),   # Theta grid minimum in radians
    "T_MAX": np.radians(90.0),  # Theta grid maximum in radians
    "N_T": 96,                  # Number of theta bins
    "s_theta": 0.95,            # Theta scaling factor

    # Temperature and radiation parameters
    "T_x": 160000.0,            # Temperature in Kelvin (unused in current runs)
    "f_uv": 0.9,                # UV flux proportion
    "f_x": 0.1,                 # X-ray flux proportion
    "L_star": 9.05e34,          # Luminosity of the central source in erg/s
    "rad_force": 1,             # Radiation force enabled (1=True, 0=False)

    # Density parameters
    "RHO_0": 1e-9,              # Midplane density in g/cm^3
    "RHO_ALPHA": 0.0,           # Midplane density drop-off (unused)
    "R_0": 8.31e8,              # Density scaling radius in cm

    # Physical constants
    "MU": 0.6,                  # Mean molecular mass (atomic mass units)
    "T_ISO": 40000.0,           # Isothermal temperature in Kelvin
    "ALPHARAD": -0.6,           # Alpha for radiation force
    "KDAR": 0.59,               # Radiation force factor
    "DFLOOR": 1e-24,            # Density floor in g/cm^3
    "GAMMA": 5.0 / 3.0,         # Adiabatic index

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
    "wind_radiation": "yes",     # Wind radiation enabled
    "line_trans": "macro",      # Line transfer mode
    "SIROCCO_VER":"-0.1"
}

# End of script
