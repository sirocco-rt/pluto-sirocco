#!/usr/bin/env python

import subprocess
import glob
from astropy.io import ascii
from astropy.table import Table
from astropy import constants as c
from astropy import units as u
from scipy.integrate import quad
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy import interpolate

import pyPLUTO.pload as pp
import numpy as np


# -----------------------------------------------------------------------------
# Script: pluto_sirocco_sub.py
# Author: Nick Higginbottom
# Revised By: Amin Mosallanezhad (a.mosallanezhad@soton.ac.uk)
# Last Modified: 2024-12-20
# -----------------------------------------------------------------------------


def get_units(fname='definitions.h'):
    """
    Retrieves code units from the definitions.h file.
    """
    UNIT_DENSITY = UNIT_LENGTH = UNIT_VELOCITY = 1.0
    with open(fname, 'r') as inp:
        for line in inp:
            data = line.split()
            if len(data) > 1:
                if data[1] == 'UNIT_DENSITY':
                    UNIT_DENSITY = float(data[2])
                elif data[1] == 'UNIT_LENGTH':
                    UNIT_LENGTH = float(data[2])
                elif data[1] == 'UNIT_VELOCITY':
                    UNIT_VELOCITY = float(data[2])
    return UNIT_DENSITY, UNIT_LENGTH, UNIT_VELOCITY



# This subroutine makes a pluto input file
def pluto_input_file(tlim, data):
    """
    Creates a PLUTO input file (pluto.ini) based on provided parameters.

    Parameters:
    tlim (float): The stopping time for the simulation.
    data (dict): Dictionary containing simulation parameters.
    """
    output = open('pluto.ini', 'w')
    output.write("[Grid]\n\n")
    # output.write("X1-grid 2 {} 184 s  4.0  72 u {}\n".format(data["R_MIN"], data["R_MAX"]))
    # output.write("X2-grid 2 {} 72  u  0.785398163397448  184  s {}\n".format(data["T_MIN"], data["T_MAX"]))
    output.write("X1-grid 1 {} {} r {}  {}\n".format(data["R_MIN"], data["N_R"], data["R_MAX"], data["s_r"]))
    output.write("X2-grid 1 {} {} r {}  {}\n".format(data["T_MIN"], data["N_T"], data["T_MAX"], data["s_theta"]))
    output.write("X3-grid 1    0.0    1      u    1.0\n\n")
    output.write("[Chombo Refinement]\n\n")
    output.write("Levels           4\n")
    output.write("Ref_ratio        2 2 2 2 2\n")
    output.write("Regrid_interval  2 2 2 2 \n")
    output.write("Refine_thresh    0.3\n")
    output.write("Tag_buffer_size  3\n")
    output.write("Block_factor     8\n")
    output.write("Max_grid_size    64\n")
    output.write("Fill_ratio       0.75\n\n")
    output.write("[Time]\n\n")
    output.write("CFL              0.4\n")
    output.write("CFL_max_var      1.1\n")
    output.write("tstop            {}\n".format(tlim))
    output.write("first_dt         1e-4\n\n")
    output.write("[Solver]\n\n")
    output.write("Solver         hll\n\n")
    output.write("[Boundary]\n\n")
    output.write("X1-beg        outflow\n")
    output.write("X1-end        outflow\n")
    output.write("X2-beg        polaraxis\n")
    output.write("X2-end        reflective\n")
    output.write("X3-beg        outflow\n")
    output.write("X3-end        outflow\n\n")
    output.write("[Static Grid Output]\n\n")
    if data["rad_force"] == 1:
        output.write("uservar    23 XI T T_r ch cc lc bc xh ch_pre cc_pre lc_pre bc_pre xh_pre XI_pre ne nh gr gt gp dv_ds t M_max1 M_max2\n")
    else:
        output.write("uservar    14    XI T ch cc lc bc xh ch_pre cc_pre lc_pre bc_pre xh_pre ne nh\n")
    output.write("dbl        1000000000000   -1   single_file\n")
    output.write("flt       -1.0  -1   single_file\n")
    output.write("vtk        1000000000000  -1   single_file\n")
    output.write("dbl.h5    -1.0  -1\n")
    output.write("flt.h5    -1.0  -1\n")
    output.write("tab       -1.0  -1   \n")
    output.write("ppm       -1.0  -1   \n")
    output.write("png       -1.0  -1\n")
    output.write("log        1000\n")
    output.write("analysis  -1.0  -1\n\n")
    output.write("[Chombo HDF5 output]\n\n")
    output.write("Checkpoint_interval  -1.0  0\n")
    output.write("Plot_interval         1.0  0 \n\n")
    output.write("[Parameters]\n\n")
    output.write("MU                          {:.2f}\n".format(data["MU"]))
    output.write("RHO_0                       {:4.2e}\n".format(data["RHO_0"]))
    output.write("R_0                         {:4.2e}\n".format(data["R_0"]))
    output.write("RHO_ALPHA                   {:.2f}\n".format(data["RHO_ALPHA"]))
    output.write("CENT_MASS                   {:4.2e}\n".format(data["CENT_MASS"]))
    output.write("DISK_MDOT                   {:4.2e}\n".format(data["DISK_MDOT"]))
    output.write("T_ISO                       {:4.2e}\n".format(data["T_ISO"]))
    output.write("L_star                      {:4.2e}\n".format(data["L_star"]))
    output.write("f_x                         {:.2f}\n".format(data["f_x"]))
    output.write("f_uv                        {:.2f}\n".format(data["f_uv"]))
    output.write("T_x                         {:4.2e}\n".format(data["T_x"]))
    output.write("KRAD                        {:4.2e}\n".format(data["k"]))
    output.write("ALPHARAD                    {:4.2e}\n".format(data["alpha"]))
    output.write("DFLOOR                      {:4.2e}\n".format(data["DFLOOR"]))
    output.write("GAMMA                       {:4.2e}\n".format(data["GAMMA"]))
    output.close()
    return


# -----------------------------------------------------------------------------
# Subroutine: pluto2sir_rtheta
# Purpose  : Converts PLUTO simulation output into a SIROCCO-readable file
#            format, and recalculating cartesian velocity
#            components. The processed data are then flattened and
#            written to an output file with a prescribed format and precision.
#
# Notes    :
#   - Ghost cells are added in both the r and theta directions.
#   - Density and temperature in ghost cells are set to zero.
#   - Velocities are converted to SIROCCO-readable Cartesian components.
#   - The output file is formatted with specific column headers and precision.
#
# -----------------------------------------------------------------------------
def pluto2sir_rtheta(ifile):
    """
    Converts a PLUTO .dbl file into a file that can be read in by SIROCCO.

    Parameters:
    ifile (int): The file index number of the PLUTO output.

    Returns:
    fname (str): The filename of the generated SIROCCO-readable file.
    """
    D = pp.pload(ifile)
    UNIT_DENSITY, UNIT_LENGTH, UNIT_VELOCITY = get_units('definitions.h')

    # Extract necessary data from the PLUTO file and scale accordingly
    pluto_r_inner = D.x1r * UNIT_LENGTH  # Inner radius of each shell (cell edge positions)
    pluto_theta_inner = D.x2r            # Theta coordinates of inner edges of cell

    pluto_r_c = D.x1 * UNIT_LENGTH       # Central point of each PLUTO cell - where the velocity is defined
    pluto_theta_c = D.x2                 # Central theta of each PLUTO cell

    pluto_vr_c = D.vx1 * UNIT_VELOCITY   # Radial velocity in each PLUTO cell
    pluto_vt_c = D.vx2 * UNIT_VELOCITY   # Theta velocity in each PLUTO cell
    pluto_vy_c = D.vx3 * UNIT_VELOCITY   # Phi velocity in each PLUTO cell

    pluto_density = D.rho * UNIT_DENSITY
    pluto_temperature = D.T

    # Set up arrays to hold velocities in x and z components, as expected by Python code
    pluto_vx_c = np.zeros(np.shape(pluto_vr_c))
    pluto_vz_c = np.zeros(np.shape(pluto_vr_c))

    # Set up SIROCCO r and theta arrays, extended to include ghost cells
    sirocco_r = list(pluto_r_inner)
    sirocco_r.append(pluto_r_inner[-1] + (pluto_r_inner[-1] - pluto_r_inner[-2]))

    sirocco_theta = list(pluto_theta_inner)
    sirocco_theta.append(pluto_theta_inner[-1] + (pluto_theta_inner[-1] - pluto_theta_inner[-2]))

    # Calculate x and z velocities on the PLUTO grid from r and theta components
    for i in range(len(pluto_r_c)):
        for j in range(len(pluto_theta_c)):
            pluto_vx_c[i][j] = (pluto_vr_c[i][j] * np.sin(pluto_theta_c[j]) +
                                pluto_vt_c[i][j] * np.cos(pluto_theta_c[j]))
            pluto_vz_c[i][j] = (pluto_vr_c[i][j] * np.cos(pluto_theta_c[j]) -
                                pluto_vt_c[i][j] * np.sin(pluto_theta_c[j]))

    # Interpolate velocities onto the SIROCCO grid
    # First interpolate in r direction
    vx_temp = np.zeros([len(sirocco_r), len(pluto_theta_c)])
    vy_temp = np.zeros([len(sirocco_r), len(pluto_theta_c)])
    vz_temp = np.zeros([len(sirocco_r), len(pluto_theta_c)])

    for i in range(len(pluto_theta_c)):
        vx = interpolate.interp1d(pluto_r_c, pluto_vx_c[:, i], fill_value='extrapolate')
        vy = interpolate.interp1d(pluto_r_c, pluto_vy_c[:, i], fill_value='extrapolate')
        vz = interpolate.interp1d(pluto_r_c, pluto_vz_c[:, i], fill_value='extrapolate')
        vx_temp[:, i] = vx(sirocco_r)
        vy_temp[:, i] = vy(sirocco_r)
        vz_temp[:, i] = vz(sirocco_r)

    # Now interpolate in theta direction
    sirocco_vx = np.zeros([len(sirocco_r), len(sirocco_theta)])
    sirocco_vy = np.zeros([len(sirocco_r), len(sirocco_theta)])
    sirocco_vz = np.zeros([len(sirocco_r), len(sirocco_theta)])

    for i in range(len(sirocco_r)):
        vx = interpolate.interp1d(pluto_theta_c, vx_temp[i], fill_value='extrapolate')
        vy = interpolate.interp1d(pluto_theta_c, vy_temp[i], fill_value='extrapolate')
        vz = interpolate.interp1d(pluto_theta_c, vz_temp[i], fill_value='extrapolate')
        sirocco_vx[i] = vx(sirocco_theta)
        sirocco_vy[i] = vy(sirocco_theta)
        sirocco_vz[i] = vz(sirocco_theta)

    # Deal with cell-centered values, density and temperature
    sirocco_density = np.zeros([len(sirocco_r), len(sirocco_theta)])  # Includes ghost cells
    sirocco_T_e = np.zeros([len(sirocco_r), len(sirocco_theta)])

    # Assign PLUTO density and temperature to SIROCCO arrays, excluding ghost cells
    sirocco_density[0:-2, 0:-2] = pluto_density
    sirocco_T_e[0:-2, 0:-2] = pluto_temperature


    # Flatten arrays into lists for output
    ir = []
    r = []
    itheta = []
    theta = []
    inwind = []
    v_x = []
    v_y = []
    v_z = []
    density = []
    T = []

    for i in range(len(sirocco_r)):
        for j in range(len(sirocco_theta)):
            ir.append(i)
            r.append(sirocco_r[i])
            itheta.append(j)
            theta.append(np.degrees(sirocco_theta[j]))
            if sirocco_density[i][j] == 0.0:
                inwind.append(-1)  # Ghost cells
            else:
                inwind.append(0)  # Wind cells
            v_x.append(sirocco_vx[i][j])
            v_y.append(sirocco_vy[i][j])
            v_z.append(sirocco_vz[i][j])
            density.append(sirocco_density[i][j])
            T.append(sirocco_T_e[i][j])

    titles = ["ir", "itheta", "inwind", "r", "theta", "v_x", "v_y", "v_z", "density", "T"]

    # Define formats for the output variables
    fmt = '%13.6e'
    fmts = {'ir': '%03i',
            'itheta': '%03i',
            'inwind': '%01i',
            'r': fmt,
            'theta': fmt,
            'v_x': fmt,
            'v_y': fmt,
            'v_z': fmt,
            'density': fmt,
            'T': fmt}

    # Set the filename
    fname = "{:08d}.pluto".format(ifile)
    out = open(fname, 'w')

    # Output the data to file
    out_dat = Table([ir, itheta, inwind, r, theta, v_x, v_y, v_z, density, T], names=titles)
    ascii.write(out_dat, out, formats=fmts)
    out.close()
    return fname



# This makes a SIROCCO input file - if the version of SIROCCO changes such that new parameters are required, this must be edited.
def sirocco_input_file(fname, data, cycles):
    """
    Creates a SIROCCO input file (.pf) for the simulation based on provided parameters.

    Parameters:
    fname (str): Base filename to use for the input file.
    data (dict): Dictionary containing simulation parameters.
    cycles (int): Number of ionization cycles.
    """
    output = open(fname + ".pf", 'w')
    output.write("System_type(star,binary,agn,previous)          " + data["system_type"] + "  \n")
    output.write("\n")
    output.write("### Parameters for the Central Object\n")
    output.write("Central_object.mass(msol)                  " + str(data["CENT_MASS"] / c.M_sun.cgs.value) + "\n")
    output.write("Central_object.radius(cm)                  " + str(data["CENT_RADIUS"]) + "\n")
    output.write("\n")
    output.write("### Parameters for the Disk (if there is one)\n")
    output.write("\n")
    if data["disk_radiation"] == "yes":
        output.write("Disk.type(none,flat,vertically.extended)       flat\n")
        output.write("Disk.radiation(yes,no)      yes\n")
        output.write("Disk.temperature.profile(standard,readin,yso) standard\n")
        # output.write("Disk.T_profile_file() max_40K_disk.dat \n")  # Uncomment if using a user-defined disk temperature file
        output.write("Disk.rad_type_to_make_wind(bb,models) bb\n")
        output.write("Disk.mdot(msol/yr) " + str(data["SIR_DISK_MDOT"]) + "\n")
        output.write("Disk.radmax(cm) " + str(data["DISK_TRUNC_RAD"]) + "\n")
    else:
        output.write("Disk.radiation(yes,no)      no\n")
    output.write("\n")
    output.write("### Parameters for BL or AGN\n")
    output.write("\n")
    if data["boundary_layer"] == "no":
        output.write("Boundary_layer.radiation(yes,no)                   no\n")
    else:
        output.write("Boundary_layer.radiation(yes,no)                   yes\n")
        output.write("Boundary_layer.rad_type_to_make_wind(bb,models,power) bb\n")
        output.write("Boundary_layer.temp(K) " + str(data["T_BL"]) + "\n")
        output.write("Boundary_layer.luminosity(ergs/s)  " + str(data["L_BL"]) + "\n")
    if data["cent_spectype"] == "brem":
        output.write("Central_object.radiation(yes,no)     yes\n")
        output.write("Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems)   brems\n")
        output.write("AGN.bremsstrahlung_temp(K) " + str(data["T_x"]) + "\n")
        output.write("Central_object.luminosity(ergs/s) " + str(data["L_2_10"]) + "\n")
        output.write("Central_object.geometry_for_source(sphere,lamp_post) sphere\n")
    elif data["cent_spectype"] == "bb":
        output.write("Central_object.radiation(yes,no)     yes\n")
        output.write("Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems)   bb\n")
        output.write("Central_object.blackbody_temp(K)                        " + str(data["T_star"]) + "\n")
        output.write("Central_object.geometry_for_source(sphere,lamp_post) sphere\n")
    elif data["cent_spectype"] == "models":
        output.write("Central_object.radiation(yes,no)     yes\n")
        output.write("Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems)   models\n")
        output.write("Input_spectra.model_file          model.ls\n")
        output.write("Central_object.luminosity(ergs/s) " + str(data["L_2_10"]) + "\n")
        output.write("Central_object.geometry_for_source(sphere,lamp_post) sphere\n")
    elif data["cent_spectype"] == "none":
        output.write("Central_object.radiation(yes,no)     no\n")
    output.write("\n")
    output.write("### Parameters describing the various winds or coronae in the system\n")
    output.write("\n")
    if data["wind_radiation"] == "yes":
        output.write("Wind.radiation(yes,no) yes\n")
    else:
        output.write("Wind.radiation(yes,no) no\n")
    output.write("Wind.number_of_components  1\n")
    output.write("Wind.type(SV,star,hydro,corona,kwd,homologous,yso,shell,imported)  imported \n")
    output.write("Wind.coord_system(spherical,cylindrical,polar,cyl_var)  polar\n")
    output.write("Wind.dim.in.x_or_r.direction               30\n")
    output.write("Wind.dim.in.z_or_theta.direction           30\n")
    output.write("\n")
    output.write("### Parameters associated with photon number, cycles, ionization and radiative transfer options\n")
    output.write("\n")
    output.write("Photons_per_cycle        " + str(data["NPHOT"]) + "\n")
    output.write("Ionization_cycles        " + str(cycles) + "\n")
    output.write("Spectrum_cycles          0\n")
    output.write("Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow)  matrix_pow\n")
    if data["line_trans"] == "macro":
        output.write("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms,macro_atoms_thermal_trapping)   macro_atoms_thermal_trapping\n")
        output.write("Atomic_data  zdata/h10_hetop_standard80.dat\n")
        output.write("Matom_transition_mode(mc_jumps,matrix) matrix\n")
    elif data["line_trans"] == "simple":
        output.write("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms,macro_atoms_thermal_trapping)   escape_prob\n")
        output.write("Atomic_data  data/standard80.dat\n")
    output.write("Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)    thermalized.rerad\n")
    output.write("Wind_heating.extra_processes(none,adiabatic,nonthermal,both)   none\n")
    output.write("\n")
    output.write("### Parameters for Domain 0\n")
    output.write("\n")
    output.write("Wind.model2import " + fname + ".pluto\n")
    output.write("Wind.t.init                                10000\n")
    output.write("Wind.filling_factor(1=smooth,<1=clumped)   1\n")
    output.write("\n")
    output.write("### Parameters for Reverberation Modeling (if needed)\n")
    output.write("\n")
    output.write("Reverb.type(none,photon,wind,matom)   none\n")
    output.write("\n")
    output.write("### Other parameters\n")
    output.write("\n")
    output.write("Photon_sampling.approach(T_star,cv,yso,AGN,min_max_freq,user_bands,cloudy_test,wide,logarithmic)  logarithmic\n")
    output.write("Photon_sampling.nbands                     10\n")
    output.write("Photon_sampling.low_energy_limit(eV)       0.13333\n")
    output.write("Photon_sampling.high_energy_limit(eV)      500\n")
    output.close()
    return



# This is the routine that is used to compute prefectors for the blondin heating and cooling rates
# This is a little legacy - I think it still works OK
# The idea is to take the last set of heating and cooling rates from the dbl file
# Compare them to the heating and cooling rates in SIROCCO, and export updated 'prefactors'.

def pre_calc(ifile, radforce = 0):
    """
    Computes prefactors for the Blondin heating and cooling rates.

    Parameters:
    ifile (int): The file index number.
    radforce (int): Flag to indicate if radiation force is included.

    Returns:
    odd (float): A counter for any discrepancies found.
    """
    max_change = 0.9
    max_accel_change = 0.9

    heatcool = ascii.read("{:08d}_py_heatcool.dat".format(ifile))
    D = pp.pload(ifile)

    # We need the definitions file - so we know the conversion factors.
    UNIT_DENSITY, UNIT_LENGTH, UNIT_VELOCITY = get_units('definitions.h')
    UNIT_ACCELERATION = UNIT_VELOCITY * UNIT_VELOCITY / UNIT_LENGTH

    comp_h_pre = []
    comp_c_pre = []
    xray_h_pre = []
    brem_c_pre = []
    line_c_pre = []
    xi_ion_pre = []

    odd = 0.0
    itest = 0

    for i in range(len(heatcool["rho"])):
        if (heatcool["rho"][i] / (D.rho[heatcool["i"][i]][heatcool["j"][i]] * UNIT_DENSITY)) - 1.0 > 1e-6:
            odd += 1
        nenh = D.ne[heatcool["i"][i]][heatcool["j"][i]] * D.nh[heatcool["i"][i]][heatcool["j"][i]]
        nhnh = D.nh[heatcool["i"][i]][heatcool["j"][i]] * D.nh[heatcool["i"][i]][heatcool["j"][i]]

        # Compute comp_h_pre
        ideal_prefactor = (heatcool["heat_comp"][i] / (D.ch[heatcool["i"][i]][heatcool["j"][i]] * nenh))
        change = ideal_prefactor / D.ch_pre[heatcool["i"][i]][heatcool["j"][i]]
        if change < max_change:
            change = max_change
        elif change > (1.0 / max_change):
            change = (1.0 / max_change)
        comp_h_pre.append(change * D.ch_pre[heatcool["i"][i]][heatcool["j"][i]])

        # Compute comp_c_pre
        ideal_prefactor = (heatcool["cool_comp"][i] / (D.cc[heatcool["i"][i]][heatcool["j"][i]] * nenh))
        change = ideal_prefactor / D.cc_pre[heatcool["i"][i]][heatcool["j"][i]]
        if change < max_change:
            change = max_change
        elif change > (1.0 / max_change):
            change = (1.0 / max_change)
        comp_c_pre.append(change * D.cc_pre[heatcool["i"][i]][heatcool["j"][i]])

        # Compute line_c_pre
        ideal_prefactor = (heatcool["cool_lines"][i] / (D.lc[heatcool["i"][i]][heatcool["j"][i]] * nenh))
        change = ideal_prefactor / D.lc_pre[heatcool["i"][i]][heatcool["j"][i]]
        if change < max_change:
            change = max_change
        elif change > (1.0 / max_change):
            change = (1.0 / max_change)
        line_c_pre.append(change * D.lc_pre[heatcool["i"][i]][heatcool["j"][i]])

        # Compute brem_c_pre
        ideal_prefactor = (heatcool["cool_ff"][i] / (D.bc[heatcool["i"][i]][heatcool["j"][i]] * nenh))
        change = ideal_prefactor / D.bc_pre[heatcool["i"][i]][heatcool["j"][i]]
        if change < max_change:
            change = max_change
        elif change > (1.0 / max_change):
            change = (1.0 / max_change)
        brem_c_pre.append(change * D.bc_pre[heatcool["i"][i]][heatcool["j"][i]])

        # Compute xray_h_pre
        ideal_prefactor = (heatcool["heat_xray"][i] / (D.xh[heatcool["i"][i]][heatcool["j"][i]] * nhnh))
        change = ideal_prefactor / D.xh_pre[heatcool["i"][i]][heatcool["j"][i]]
        if change < max_change:
            change = max_change
        elif change > (1.0 / max_change):
            change = (1.0 / max_change)
        xray_h_pre.append(change * D.xh_pre[heatcool["i"][i]][heatcool["j"][i]])

        # Compute xi_ion_pre
        ideal_prefactor = (heatcool["xi"][i] / (D.XI[heatcool["i"][i]][heatcool["j"][i]]))
        change = ideal_prefactor / D.XI_pre[heatcool["i"][i]][heatcool["j"][i]]
        if change < max_change:
            change = max_change
        elif change > (1.0 / max_change):
            change = (1.0 / max_change)
        xi_ion_pre.append(change * D.XI_pre[heatcool["i"][i]][heatcool["j"][i]])

    fmt = '%013.6e'

    # Define formats for the output variables
    fmts = {'ir': '%03i',
            'rcent': fmt,
            'itheta': '%03i',
            'thetacent': fmt,
            'rho': fmt,
            'comp_h_pre': fmt,
            'comp_c_pre': fmt,
            'xray_h_pre': fmt,
            'brem_c_pre': fmt,
            'line_c_pre': fmt,
            'xi_ion_pre': fmt,
            }

    titles = ["ir", "rcent", "itheta", "thetacent", "rho",
              "comp_h_pre", "comp_c_pre", "xray_h_pre", "brem_c_pre", "line_c_pre", "xi_ion_pre"]

    col0 = heatcool["i"]
    col1 = heatcool["rcen"]
    col2 = heatcool["j"]
    col3 = heatcool["thetacen"]
    col4 = heatcool["rho"]
    col5 = comp_h_pre
    col6 = comp_c_pre
    col7 = xray_h_pre
    col8 = brem_c_pre
    col9 = line_c_pre
    col10 = xi_ion_pre

    out = open("prefactors.dat", 'w')

    out_dat = Table([col0, col1, col2, col3, col4, col5, col6, col7, col8, col9, col10], names = titles)
    ascii.write(out_dat, out, formats = fmts)
    out.close()

    return odd
