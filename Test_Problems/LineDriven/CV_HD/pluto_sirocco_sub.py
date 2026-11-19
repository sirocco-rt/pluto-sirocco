#!/usr/bin/env python

import subprocess
import glob
from astropy.io import ascii
from astropy.table import Table
from astropy import constants as c
from astropy import units as u
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy import interpolate

import numpy as np

# =============================================================================
# pload: self-contained PLUTO .dbl file reader (replaces pyPLUTO.pload).
# Reads grid.out, dbl.out (or flt.out), and data.NNNN.dbl (single_file mode).
# Attributes set on the returned object match the pyPLUTO pload interface used
# by pluto2sir_rtheta(), pre_calc(), and the restart logic in the driver:
#   x1, x2, x3, dx1, dx2, dx3   -- cell centres and widths (code units)
#   x1r, x2r, x3r               -- cell-edge positions
#   SimTime, Dt, NStep           -- time information
#   <varname>                    -- 2-D array D.<varname>[i][j] for each
#                                   variable listed in dbl.out  (shape n1 x n2)
# =============================================================================

class pload:
    """
    Minimal PLUTO double-precision binary (.dbl) file reader.
    Replaces pyPLUTO.pload with no external dependencies beyond numpy.

    Usage (identical to pyPLUTO):
        D = pload(ifile)              # reads data.{ifile:04d}.dbl
        D = pload(ifile, w_dir='./') # explicit directory
    """

    def __init__(self, ns, w_dir=None, datatype=None):
        import os, struct, array as _array

        self.NStep = ns
        self.SimTime = 0.0
        self.Dt = 0.0

        nstr = '{:04d}'.format(ns)

        if w_dir is None:
            w_dir = os.getcwd() + '/'
        if not w_dir.endswith('/'):
            w_dir += '/'

        if datatype is None:
            datatype = 'dbl'
        self.datatype = datatype

        if datatype == 'flt':
            dtype_char = 'f'
            varfile  = w_dir + 'flt.out'
            datafile = w_dir + 'data.' + nstr + '.flt'
        else:
            dtype_char = 'd'
            varfile  = w_dir + 'dbl.out'
            datafile = w_dir + 'data.' + nstr + '.dbl'

        # ── read dbl.out / flt.out ──────────────────────────────────────────
        # Format: output_index  time  dt  hydro_step  filetype  endian  var1 var2 ...
        with open(varfile, 'r') as f:
            lines = [l for l in f if l.strip()]

        # find the line for this step number
        varinfo = None
        for line in lines:
            parts = line.split()
            if len(parts) >= 6:
                try:
                    if int(parts[0]) == ns:
                        varinfo = parts
                        break
                except ValueError:
                    pass
        if varinfo is None:
            raise ValueError("Output {} is missing from {}".format(ns, varfile))

        self.SimTime   = float(varinfo[1])
        self.Dt        = float(varinfo[2])
        self.NStep    = int(varinfo[3])
        filetype       = varinfo[4]   # 'single_file' or 'multiple_files'
        endian_token   = varinfo[5]   # 'little' or 'big'
        var_names      = varinfo[6:]

        if filetype != 'single_file':
            raise ValueError("pload supports single_file output, got " + filetype)
        if endian_token not in ('little', 'big'):
            raise ValueError("Unknown output endianness: " + endian_token)
        endian = '<' if endian_token == 'little' else '>' 

        # ── read grid.out ───────────────────────────────────────────────────
        self._read_grid(w_dir + 'grid.out')

        n1, n2, n3 = self.n1, self.n2, self.n3

        # ── read binary data ────────────────────────────────────────────────
        nptype = endian + ('f8' if dtype_char == 'd' else 'f4')
        # CV_HD output contains only cell-centred variables.
        shapes = [(n1, n2, n3)] * len(var_names)

        with open(datafile, 'rb') as fp:
            raw = fp.read()
        itemsize = np.dtype(nptype).itemsize
        expected_bytes = sum(int(np.prod(shape)) for shape in shapes) * itemsize
        if len(raw) != expected_bytes:
            raise ValueError(
                "{}: expected {} bytes from grid/output metadata, got {}".format(
                    datafile, expected_bytes, len(raw)))

        offset = 0
        for vname, shape in zip(var_names, shapes):
            count = int(np.prod(shape))
            arr_flat = np.frombuffer(raw, dtype=np.dtype(nptype),
                                     count=count, offset=offset)
            # Binary storage is [k][j][i]; expose [i][j][k], or [i][j]
            # when this variable has only one k plane.
            arr = arr_flat.reshape(shape[::-1]).transpose(2, 1, 0)
            if shape[2] == 1:
                arr = arr[:, :, 0]
            object.__setattr__(self, vname, arr.copy())
            offset += count * itemsize

    def _read_grid(self, gridfile):
        """Parse grid.out and set x1/x2/x3, dx, xr arrays."""
        xL = []
        xR = []
        nmax = []
        geometry = 'CARTESIAN'

        with open(gridfile, 'r') as gfp:
            for line in gfp:
                parts = line.split()
                if len(parts) == 1:
                    try:
                        nmax.append(int(parts[0]))
                    except ValueError:
                        pass
                elif len(parts) == 3:
                    try:
                        int(parts[0])   # index — proves it is a cell line
                        xL.append(float(parts[1]))
                        xR.append(float(parts[2]))
                    except ValueError:
                        if len(parts) >= 3 and parts[1] == 'GEOMETRY:':
                            geometry = parts[2]

        self.geometry = geometry
        self.n1, self.n2, self.n3 = nmax[0], nmax[1], nmax[2]
        n1 = self.n1
        n1p2 = n1 + self.n2
        n1p2p3 = n1p2 + self.n3

        self.x1  = np.array([0.5*(xL[i]+xR[i]) for i in range(n1)])
        self.dx1 = np.array([xR[i]-xL[i]        for i in range(n1)])
        self.x2  = np.array([0.5*(xL[i]+xR[i]) for i in range(n1, n1p2)])
        self.dx2 = np.array([xR[i]-xL[i]        for i in range(n1, n1p2)])
        self.x3  = np.array([0.5*(xL[i]+xR[i]) for i in range(n1p2, n1p2p3)])
        self.dx3 = np.array([xR[i]-xL[i]        for i in range(n1p2, n1p2p3)])

        # cell-edge arrays (same convention as pyPLUTO)
        self.x1r = np.empty(n1 + 1)
        self.x1r[1:] = self.x1 + self.dx1 / 2.0
        self.x1r[0]  = self.x1r[1] - self.dx1[0]

        self.x2r = np.empty(self.n2 + 1)
        self.x2r[1:] = self.x2 + self.dx2 / 2.0
        self.x2r[0]  = self.x2r[1] - self.dx2[0]

        self.x3r = np.empty(self.n3 + 1)
        self.x3r[1:] = self.x3 + self.dx3 / 2.0
        self.x3r[0]  = self.x3r[1] - self.dx3[0]



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
    output.write("X1-grid 1 {} {} r {}  {}\n".format(data["R_MIN"], data["N_R"], data["s_r"], data["R_MAX"]))
    output.write("X2-grid 1 {} {} r {}  {}\n".format(data["T_MIN"], data["N_T"], data["s_theta"], data["T_MAX"]))
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
    output.write("CFL              0.3\n")
    output.write("CFL_max_var      1.1\n")
    output.write("tstop            {}\n".format(tlim))
    output.write("first_dt         1e-8\n\n")
    output.write("[Solver]\n\n")
    output.write("Solver         hll\n\n")
    output.write("[Boundary]\n\n")
    output.write("X1-beg        userdef\n")
    output.write("X1-end        userdef\n")
    output.write("X2-beg        userdef\n")
    output.write("X2-end        eqtsymmetric\n")
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
    output.write("DFLOOR                      {:4.2e}\n".format(data["DFLOOR"]))
    output.write("GAMMA                       {:0.5f}\n".format(data["GAMMA"]))
    output.write("L_STAR                      {:4.2e}\n".format(data["L_STAR"]))
    output.write("F_X                         {:4.2e}\n".format(data["F_X"]))
    output.write("T_X                         {:4.2e}\n".format(data["T_X"]))
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
    D = pload(ifile)
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

    # Hand-craft header and separator to exactly match the SIROCCO import format:
    #   i(3) j(3) inwind(6) r(12) theta(12) v_x(13) v_y(12) v_z(13) rho(12) t_e(12)
    # Each label is right-justified within its column width; fields are separated
    # by a single space, matching the reference Py_2d.txt layout (without t_r).
    cols = [
        ("i",      3),
        ("j",      3),
        ("inwind", 6),
        ("r",      12),
        ("theta",  12),
        ("v_x",    13),
        ("v_y",    12),
        ("v_z",    13),
        ("rho",    12),
        ("t_e",    12),
    ]
    header    = " ".join(name.rjust(width) for name, width in cols)
    separator = " ".join("-" * width       for _,    width in cols)

    # Integer format widths: i->3, j->3, inwind->6
    # Float format: 12-char scientific notation (%12.6e) for all float columns.
    # v_x and v_z are 13-char to accommodate the leading minus sign.
    fmt_i      = "%3d"
    fmt_j      = "%3d"
    fmt_inwind = "%6d"
    fmt_float  = "%12.6e"
    fmt_signed = "%13.6e"   # v_x and v_z can be negative, need one extra char

    fname = "{:08d}.pluto".format(ifile)
    with open(fname, "w") as out:
        out.write(header + "\n")
        out.write(separator + "\n")
        for k in range(len(ir)):
            out.write(
                (fmt_i      % ir[k])      + " " +
                (fmt_j      % itheta[k])  + " " +
                (fmt_inwind % inwind[k])  + " " +
                (fmt_float  % r[k])       + " " +
                (fmt_float  % theta[k])   + " " +
                (fmt_signed % v_x[k])     + " " +
                (fmt_float  % v_y[k])     + " " +
                (fmt_signed % v_z[k])     + " " +
                (fmt_float  % density[k]) + " " +
                (fmt_float  % T[k])       + "\n"
            )
    return fname



# -----------------------------------------------------------------------------
# Subroutine: atomic_data_file
# Purpose  : Single source of truth for the Atomic_data path written into the
#            .pf files.  Both sirocco_input_file() (ionisation+spectral runs)
#            and sirocco_spec_input_file() (spectral-only runs) call this, so
#            the two can never drift apart.
#
#            This matters: a wind_save is written with a particular atomic
#            dataset, and a spectral-only run restarting from that wind_save
#            must use the SAME dataset or the ion/level structure will not
#            correspond.
#
# Returns  : the relative path, e.g. "data/standard80.dat", exactly as it must
#            appear in the .pf.  The 'data' directory is the one created in the
#            run directory by Setup_Sirocco_Dir on Iridis.
# -----------------------------------------------------------------------------
def atomic_data_file(data):
    """Return the Atomic_data path implied by data['line_trans']."""
    if data["line_trans"] == "macro":
        return "data/h20_hetop_standard80.dat"
    elif data["line_trans"] == "simple":
        return "data/standard80.dat"
    else:
        raise ValueError(
            "atomic_data_file: unrecognised line_trans '{}' "
            "(expected 'macro' or 'simple')".format(data["line_trans"]))


# This makes a SIROCCO input file - if the version of SIROCCO changes such that new parameters are required, this must be edited.
def wind_filling_factor(data):
    """Read the configured volume filling factor (1=smooth, <1=clumped)."""
    import math
    value = float(data["FILLING_FACTOR"])
    if not math.isfinite(value) or not 0.0 < value <= 1.0:
        raise ValueError("FILLING_FACTOR must be finite and satisfy 0 < f <= 1")
    return value


def sirocco_input_file(fname, data, cycles, model_file=None):
    """
    Creates a SIROCCO input file (.pf) for the simulation based on provided parameters.

    Parameters:
    fname (str): Base filename to use for the input file (root name for .pf).
    data (dict): Dictionary containing simulation parameters.
    cycles (int): Number of ionization cycles.
    model_file (str): The .pluto model file for Wind.model2import.
                      If None, defaults to fname (legacy behaviour).
    """
    filling_factor = wind_filling_factor(data)
    if model_file is None:
        model_file = fname
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
        output.write("AGN.bremsstrahlung_temp(K) " + str(data["T_X"]) + "\n")
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
    # output.write("Wind.dim.in.x_or_r.direction               30\n")
    # output.write("Wind.dim.in.z_or_theta.direction           30\n")
    output.write("\n")
    output.write("### Parameters associated with photon number, cycles, ionization and radiative transfer options\n")
    output.write("\n")
    output.write("Photons_per_cycle        " + str(data["NPHOT"]) + "\n")
    output.write("Ionization_cycles        " + str(cycles) + "\n")
    output.write("Spectrum_cycles          0\n")
    output.write("Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow)  matrix_pow\n")
    if data["line_trans"] == "macro":
        output.write("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms,macro_atoms_thermal_trapping)   macro_atoms_thermal_trapping\n")
        output.write("Atomic_data  " + atomic_data_file(data) + "\n")
        output.write("Matom_transition_mode(mc_jumps,matrix) matrix\n")
    elif data["line_trans"] == "simple":
        output.write("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms,macro_atoms_thermal_trapping)   escape_prob\n")
        output.write("Atomic_data  " + atomic_data_file(data) + "\n")
    output.write("Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)    thermalized.rerad\n")
    output.write("Wind_heating.extra_processes(none,adiabatic,nonthermal,both)   none\n")
    output.write("\n")
    output.write("### Parameters for Domain 0\n")
    output.write("\n")
    output.write("Wind.model2import " + model_file + "\n")
    output.write("Wind.t.init                                10000\n")
    output.write("Wind.filling_factor(1=smooth,<1=clumped)    " + str(filling_factor) + "\n")
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


# -----------------------------------------------------------------------------
# Subroutine: sirocco_spec_input_file
# Purpose  : Write a SIROCCO parameter file that runs SPECTRAL CYCLES ONLY,
#            starting from an existing wind_save file rather than re-converging
#            the ionisation state.
#
#            This is the "Knigge test": the ionisation state in the wind_save
#            is the one SIROCCO converged for the RHD (PLUTO) snapshot, so the
#            spectra produced here can be compared directly against those from
#            the normal ionisation+spectral run to check for gross differences.
#
# Key differences from sirocco_input_file():
#   - System_type is 'previous', with Wind.old_windfile pointing at the saved
#     wind_save root.  All wind geometry, densities, velocities and ion
#     populations come from that file.
#   - Ionization_cycles is 0.  Spectrum_cycles is data["SPEC_CYCLES"].
#   - No Wind.type / Wind.model2import / Wind.coord_system / Wind.t.init /
#     Wind.filling_factor lines: those describe how to BUILD a wind, and in
#     'previous' mode the wind already exists. Its filling factor was set
#     from data["FILLING_FACTOR"] when that wind was built; changing the
#     config cannot change an existing wind_save without reconvergence.
#   - A full Spectrum.* block is written (absent from the live-loop .pf).
#   - Wind_heating.extra_processes defaults to 'adiabatic' here, whereas the
#     live loop uses 'none'.
#
# Notes:
#   - fname is the ROOT of the parameter file (e.g. "spec0100" -> spec0100.pf).
#     SIROCCO names all its outputs after this root, so it MUST differ from
#     windsave_root, otherwise the run would overwrite the very wind_save it
#     is reading.
#   - Atomic_data is omitted by default (data["SPEC_ATOMIC_DATA"] = None) so
#     that SIROCCO takes the dataset from the wind_save and cannot mismatch it.
#     Set it to "auto" to write the dataset implied by data["line_trans"].
# -----------------------------------------------------------------------------
def sirocco_spec_input_file(fname, data, windsave_root, cycles=None):
    """
    Creates a spectral-cycles-only SIROCCO input file (.pf).

    Parameters:
    fname (str): Root name of the parameter file to write (no .pf extension).
                 Must NOT equal windsave_root.
    data (dict): Dictionary containing simulation parameters.
    windsave_root (str): Root of the existing wind_save file to restart from,
                 i.e. <windsave_root>.wind_save must exist in the run directory.
    cycles (int): Number of spectral cycles.  Defaults to data["SPEC_CYCLES"].
    """
    filling_factor = wind_filling_factor(data)
    if cycles is None:
        cycles = data.get("SPEC_CYCLES", 30)

    if fname == windsave_root:
        raise ValueError(
            "sirocco_spec_input_file: parameter-file root ({}) must differ "
            "from the wind_save root, or SIROCCO will overwrite the input "
            "wind_save with its own output.".format(fname))

    nphot       = data.get("SPEC_NPHOT", data["NPHOT"])
    extra_proc  = data.get("SPEC_EXTRA_PROC", "adiabatic")
    nwave       = data.get("SPEC_NWAVE", 10000)
    wavemin     = data.get("SPEC_WAVEMIN", 1450)
    wavemax     = data.get("SPEC_WAVEMAX", 1650)
    angles      = data.get("SPEC_ANGLES",
                           [5, 10, 15, 20, 25, 30, 35, 40, 45,
                            50, 55, 60, 65, 70, 75, 80, 85])
    # Atomic_data path for the spectral .pf.
    #   "auto"  (default) -> same path the live .pf uses, via atomic_data_file()
    #   None / "none"     -> omit the line; SIROCCO takes it from the wind_save
    #   any other string  -> written verbatim (explicit override)
    atomic_mode = data.get("SPEC_ATOMIC_DATA", "auto")
    if atomic_mode is None or str(atomic_mode).lower() == "none":
        atomic_path = None
    elif str(atomic_mode).lower() == "auto":
        atomic_path = atomic_data_file(data)
    else:
        atomic_path = str(atomic_mode)

    output = open(fname + ".pf", 'w')

    output.write("### Spectral-cycles-only run restarted from an existing wind_save\n")
    output.write("### Ionisation state is taken from {}.wind_save (RHD snapshot)\n".format(windsave_root))
    output.write("### and is NOT re-converged: Ionization_cycles is 0.\n")
    output.write("### Configured FILLING_FACTOR = {}; actual value is inherited from wind_save.\n".format(filling_factor))
    output.write("### Rebuild/reconverge the wind after changing FILLING_FACTOR.\n")
    output.write("\n")
    output.write("System_type(star,binary,agn,previous)          previous\n")
    output.write("Wind.old_windfile(root_only)                  " + windsave_root + "\n")
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
        output.write("AGN.bremsstrahlung_temp(K) " + str(data["T_X"]) + "\n")
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
    output.write("###\n")
    output.write("### NOTE: no Wind.type / Wind.model2import / Wind.coord_system here.\n")
    output.write("### In 'previous' mode the wind is read from the old wind_save.\n")
    output.write("\n")
    if data["wind_radiation"] == "yes":
        output.write("Wind.radiation(yes,no) yes\n")
    else:
        output.write("Wind.radiation(yes,no) no\n")
    output.write("\n")
    output.write("### Parameters associated with photon number, cycles, ionization and radiative transfer options\n")
    output.write("\n")
    output.write("Photons_per_cycle        " + str(nphot) + "\n")
    output.write("Ionization_cycles        0\n")
    output.write("Spectrum_cycles          " + str(cycles) + "\n")
    output.write("Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow)  matrix_pow\n")
    if data["line_trans"] == "macro":
        output.write("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms,macro_atoms_thermal_trapping)   macro_atoms_thermal_trapping\n")
        if atomic_path is not None:
            output.write("Atomic_data  " + atomic_path + "\n")
        output.write("Matom_transition_mode(mc_jumps,matrix) matrix\n")
    elif data["line_trans"] == "simple":
        output.write("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms,macro_atoms_thermal_trapping)   escape_prob\n")
        if atomic_path is not None:
            output.write("Atomic_data  " + atomic_path + "\n")
    output.write("Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)    thermalized.rerad\n")
    output.write("Wind_heating.extra_processes(none,adiabatic,nonthermal,both)   " + extra_proc + "\n")
    output.write("\n")
    output.write("### Parameters defining the spectra seen by observers\n")
    output.write("\n")
    output.write("Disk.rad_type_in_final_spectrum(bb,models,uniform)                   bb\n")
    output.write("\n")
    output.write("### The minimum and maximum wavelengths in the final spectra and the number of wavelength bins\n")
    output.write("Spectrum.nwave                                " + str(nwave) + "\n")
    output.write("Spectrum.wavemin(Angstroms)                   " + str(wavemin) + "\n")
    output.write("Spectrum.wavemax(Angstroms)                   " + str(wavemax) + "\n")
    output.write("\n")
    output.write("### The observers and their location relative to the system\n")
    output.write("Spectrum.no_observers                         " + str(len(angles)) + "\n")
    for ang in angles:
        output.write("Spectrum.angle(0=pole)                        " + str(ang) + "\n")
    output.write("Spectrum.live_or_die(live.or.die,extract)              extract\n")
    output.write("Spectrum.type(flambda,fnu,basic)              flambda\n")
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


# -----------------------------------------------------------------------------
# Subroutine: spec_cycle_selected
# Purpose  : Decide whether cycle i should get a spectral-only run, according
#            to SPEC_START / SPEC_END / SPEC_STEP in the config.
#
# Example  : SPEC_START=100, SPEC_END=500, SPEC_STEP=5
#            -> True for 100, 105, 110, ... , 495, 500 and False elsewhere.
# -----------------------------------------------------------------------------
def spec_cycle_selected(i, data):
    """Return True if cycle i is one of the selected spectral-test snapshots."""
    if str(data.get("SPEC_RUN", "no")).lower() != "yes":
        return False

    start = data.get("SPEC_START", 0)
    end   = data.get("SPEC_END", -1)
    step  = data.get("SPEC_STEP", 1)

    if step <= 0:
        return False
    if i < start or i > end:
        return False

    return ((i - start) % step) == 0



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

    # AM - 2609 - Safety cap on the absolute size of a prefactor.
    #
    # In cells where SIROCCO reports xi (near) zero -- e.g. shadowed pockets
    # behind dense clumps -- the "bare" analytic Blondin rate (linear in xi,
    # see userdef_output.c: comp_heat = comp_h_pre * 8.9e-36 * xi * tx)
    # collapses toward zero while SIROCCO's own heat_comp/cool_comp/etc. stay
    # finite (diffuse/scattered heating the simple xi-dilution formula can't
    # represent). The "ideal" target ratio there is enormous, and even after
    # the double-division fix below, max_change=0.9 rate-limiting means the
    # prefactor still climbs by the same ~11% ceiling every single cycle this
    # persists -- 130+ consecutive cycles is enough to reach ~1e6. This cap
    # bounds that growth to something physically sane regardless of how the
    # analytic denominator behaves, without changing anything for the
    # well-behaved cells (which sit close to 1 and never approach the cap).
    max_prefactor = 100.0
    min_prefactor = 1.0 / max_prefactor
    n_capped = 0

    def _update_prefactor(sirocco_value, D_value, denom_extra, D_pre_value, prefactor_baked_in=True):
        """Return the new, rate-limited, capped prefactor for one cell.

        sirocco_value : SIROCCO's own heat/cool/xi value for this cell (target numerator)
        D_value       : the corresponding PLUTO "ideal" field (e.g. D.ch[i][j])
        denom_extra   : nenh or nhnh, the extra density factor in the analytic rate
                         (pass 1.0 where there is none, e.g. for XI)
        D_pre_value   : the previous cycle's prefactor for this cell
        prefactor_baked_in : True for ch/cc/lc/bc/xh -- userdef_output.c writes
                         these as D_pre_OLD * bare_analytic (see comp_heat[k][j][i]
                         = comp_heat_pre[k][j][i] * (8.9e-36*xi*tx) etc.), so
                         sirocco_value / (D_value * denom_extra) is ALREADY the
                         ratio (target / D_pre_OLD) -- dividing by D_pre_value
                         again double-counts it. False for XI specifically:
                         userdef_output.c writes xi_out[k][j][i] = xi directly,
                         with no prefactor multiplication anywhere (and
                         xi_ion_pre is in fact never used again in cooling.c --
                         confirmed dead downstream either way), so the original
                         two-step division is the mathematically correct one
                         for this one field and must NOT be "simplified" away.
        """
        nonlocal n_capped
        if not np.isfinite(D_pre_value) or D_pre_value <= 0.0:
            raise ValueError("Previous cooling prefactor must be finite and positive")
        with np.errstate(divide='ignore', invalid='ignore'):
            ideal_prefactor = sirocco_value / (D_value * denom_extra)
            change = ideal_prefactor / D_pre_value if not prefactor_baked_in else ideal_prefactor
        if not np.isfinite(change):
            # 0/0 or x/0 (e.g. an empty cell): keep last cycle's prefactor
            # instead of writing NaN/inf into prefactors.dat.
            change = 1.0
        if change < max_change:
            change = max_change
        elif change > (1.0 / max_change):
            change = (1.0 / max_change)
        new_pre = change * D_pre_value
        if new_pre > max_prefactor:
            new_pre = max_prefactor
            n_capped += 1
        elif new_pre < min_prefactor:
            new_pre = min_prefactor
            n_capped += 1
        return new_pre

    # Use np.loadtxt instead of astropy ascii.read for speed.
    # ascii.read auto-detects format and uses a Python parser — on a
    # 40,000-row sirocco_heatcool.dat it can hang for several minutes.
    # np.loadtxt reads the same file in under a second.
    _hc_fname = "{:08d}_sirocco_heatcool.dat".format(ifile)
    with open(_hc_fname) as _fh:
        _hc_header = _fh.readline().split()
    _hc_arr = np.loadtxt(_hc_fname, skiprows=1)
    # Map column names to column indices from the header
    _hc_col = {name: idx for idx, name in enumerate(_hc_header)}
    class _HeatcoolTable:
        """Minimal dict-like wrapper so pre_calc column access is unchanged."""
        _INT_COLS = {"i", "j"}   # columns used as array indices — must be int
        def __init__(self, arr, col):
            self._arr = arr
            self._col = col
        def __getitem__(self, key):
            col_data = self._arr[:, self._col[key]]
            if key in self._INT_COLS:
                return col_data.astype(int)
            return col_data
        def __len__(self):
            return self._arr.shape[0]
    heatcool = _HeatcoolTable(_hc_arr, _hc_col)
    D = pload(ifile)

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

    # Pre-extract all columns used in the loop to avoid repeated __getitem__ calls
    _hc_i         = heatcool["i"]           # int array, used as PLUTO i-index
    _hc_j         = heatcool["j"]           # int array, used as PLUTO j-index
    _hc_rho       = heatcool["rho"]
    _hc_heat_comp = heatcool["heat_comp"]
    _hc_cool_comp = heatcool["cool_comp"]
    _hc_cool_lines= heatcool["cool_lines"]
    _hc_cool_ff   = heatcool["cool_ff"]
    _hc_heat_xray = heatcool["heat_xray"]
    _hc_xi        = heatcool["xi"]

    for i in range(len(_hc_rho)):
        if (_hc_rho[i] / (D.rho[_hc_i[i]][_hc_j[i]] * UNIT_DENSITY)) - 1.0 > 1e-6:
            odd += 1
        nenh = D.ne[_hc_i[i]][_hc_j[i]] * D.nh[_hc_i[i]][_hc_j[i]]
        nhnh = D.nh[_hc_i[i]][_hc_j[i]] * D.nh[_hc_i[i]][_hc_j[i]]

        # Compute comp_h_pre
        comp_h_pre.append(_update_prefactor(
            _hc_heat_comp[i], D.ch[_hc_i[i]][_hc_j[i]], nenh, D.ch_pre[_hc_i[i]][_hc_j[i]]))

        # Compute comp_c_pre
        comp_c_pre.append(_update_prefactor(
            _hc_cool_comp[i], D.cc[_hc_i[i]][_hc_j[i]], nenh, D.cc_pre[_hc_i[i]][_hc_j[i]]))

        # Compute line_c_pre
        line_c_pre.append(_update_prefactor(
            _hc_cool_lines[i], D.lc[_hc_i[i]][_hc_j[i]], nenh, D.lc_pre[_hc_i[i]][_hc_j[i]]))

        # Compute brem_c_pre
        brem_c_pre.append(_update_prefactor(
            _hc_cool_ff[i], D.bc[_hc_i[i]][_hc_j[i]], nenh, D.bc_pre[_hc_i[i]][_hc_j[i]]))

        # Compute xray_h_pre
        xray_h_pre.append(_update_prefactor(
            _hc_heat_xray[i], D.xh[_hc_i[i]][_hc_j[i]], nhnh, D.xh_pre[_hc_i[i]][_hc_j[i]]))

        # Compute xi_ion_pre -- XI is written unprefactored (xi_out = xi, no
        # xi_ion_pre multiplication anywhere in userdef_output.c), unlike the
        # five rates above, so this needs prefactor_baked_in=False. Also note
        # xi_ion_pre is read into cooling.c but never referenced again there
        # -- it has no live effect on the simulation either way.
        xi_ion_pre.append(_update_prefactor(
            _hc_xi[i], D.XI[_hc_i[i]][_hc_j[i]], 1.0, D.XI_pre[_hc_i[i]][_hc_j[i]],
            prefactor_baked_in=False))

    fmt = '%013.6e'

    # Define formats for the output variables
    fmts = {'i':          '%03i',
            'rcent':      fmt,
            'j':          '%03i',
            'thetacent':  fmt,
            'rho':        fmt,
            'comp_h_pre': fmt,
            'comp_c_pre': fmt,
            'xray_h_pre': fmt,
            'brem_c_pre': fmt,
            'line_c_pre': fmt,
            'xi_ion_pre': fmt,
            }

    titles = ["i", "rcent", "j", "thetacent", "rho",
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

    if n_capped > 0:
        print("pre_calc: {} prefactor update(s) hit the [{:g}, {:g}] safety cap "
              "this cycle (likely xi~0 shadowed/clumped cells)".format(
                  n_capped, min_prefactor, max_prefactor))

    return odd
