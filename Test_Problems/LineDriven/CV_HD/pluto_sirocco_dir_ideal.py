#!/usr/bin/env python3
# =============================================================================
# Script : pluto_sirocco_dir_ideal.py
# Purpose: Drive the iterative PLUTO <-> SIROCCO <-> CAK loop for the cv_ideal
#          (ideal EOS, BLONDIN cooling) disc-wind problem.
#
# Cycle 0 — initialisation (Run A only):
#
#   Run A  (t = 0, no flux files present)
#          If data.0000.dbl does not exist, PLUTO is run once at t=0
#          to generate the initial conditions snapshot automatically.
#          Model file: 00000000.pluto  (from pluto2sir_rtheta(0))
#          SIROCCO: init_sir_cycles = 30, with -p 2 (photon ramp).
#          rad_hydro_files produces directional_flux_{r,theta,phi}.dat.
#          Saved as directional_flux_*_init.dat (restart safety copies).
#          CAK is run.  Blondin prefactors are NOT computed here (see below).
#
# Main loop (i = 1, 2, …):
#   Cycle i=1:
#          PLUTO restarts from file 0 (-restart 0), producing data.0001.dbl
#          at t = dt_rad, with directional flux files already present so
#          lc_NFLUX_ANGLES > 0 from the first timestep.
#          hydro2sirocco updates the Run-A wind_save to the new snapshot.
#          SIROCCO uses -r (restart from wind_save, NO -p 2) with
#          total cycles = init_sir_cycles + sir_cycles (e.g. 30+2=32).
#   Cycle i>1:
#          PLUTO restarts from file i-1.
#          SIROCCO uses -r with total cycles = init_sir_cycles + i*sir_cycles.
#
# KRAD and ALPHARAD have been removed from definitions.h.
# The C code handles missing/zero flux files safely:
#   lc_NFLUX_ANGLES == 0  ->  LineForce returns zero, cooling.c uses
#   analytical xi estimate.
#
# Usage:
#   python3 pluto_sirocco_dir_ideal.py          # fresh start
#   python3 pluto_sirocco_dir_ideal.py <N>      # restore archived cycle N, continue with N+1
#
# Author : A. Mosallanezhad  (a.mosallanezhad@soton.ac.uk)
#          N. Higginbottom   (original pipeline structure)
# Date   : April 2026
# =============================================================================

import os
import sys
import glob
import time
import shutil
import warnings
import subprocess
from datetime import datetime

from pluto_sirocco_config import data
import pluto_sirocco_sub as pss



# ---------------------------------------------------------------------------
# Simulation timing
# ---------------------------------------------------------------------------

dt_rad = 2.0           # radiation time step: PLUTO advances by this many code-time units per cycle

init_sir_cycles = 25   # SIROCCO ionisation cycles in Run A (with -p 2)
sir_cycles      = 2    # additional cycles per main-loop iteration

tot_num_cycle   = 751 # totla number of cycles for main-loop iteration


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logfile = open("run_logfile", "w", buffering=1)


def log(msg: str) -> None:
    out = msg if msg.endswith("\n") else msg + "\n"
    logfile.write(out); logfile.flush()
    sys.stdout.write(out); sys.stdout.flush()


def die(msg: str) -> None:
    log("FATAL: " + msg)
    logfile.close()
    sys.exit(1)


def _short(cmd):
    """Return command as a display string with only basenames (no full paths)."""
    return " ".join(os.path.basename(str(c)) if os.sep in str(c) else str(c)
                    for c in cmd)


def tic():
    """Return current wall-clock time (seconds)."""
    return time.time()


def toc(t_start, label=""):
    """Log elapsed time since t_start.  Returns elapsed seconds."""
    elapsed = time.time() - t_start
    if label:
        log("  [timer] {} : {:.1f} s  ({:.2f} min)".format(
            label, elapsed, elapsed / 60.0))
    return elapsed


def run(cmd, logpath=None, check=True, ignore_fail=False, label="cmd"):
    """Run a subprocess, optionally redirecting stdout to a log file.
    Returns (result, elapsed_seconds).
    """
    log("  [{}] {}".format(label, _short(cmd)))
    t_start = tic()
    try:
        if logpath:
            with open(logpath, "w") as fh:
                result = subprocess.run(cmd, stdout=fh,
                                        stderr=subprocess.STDOUT)
        else:
            result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True)
            if result.stdout:
                log(result.stdout[-4000:])

        elapsed = toc(t_start, label)

        if result.returncode != 0:
            msg = "{} exited {}".format(label, result.returncode)
            if ignore_fail:
                log("WARNING: " + msg + " (ignored)")
            elif check:
                die(msg)
        return result, elapsed

    except FileNotFoundError as e:
        if ignore_fail:
            log("WARNING: {} not found ({}) — skipped".format(label, e))
            return None, 0.0
        die("{} not found: {}".format(label, e))


def sirocco_run(cmd, label="SIROCCO"):
    """Run SIROCCO, appending output to the single 'sirocco_log' file.

    A separator header is written before each run so individual cycle
    outputs remain identifiable within the combined log.
    Returns (result, elapsed_seconds).
    """
    log("  [{}] {}".format(label, _short(cmd)))
    t_start = tic()
    header = (
        "\n" + "=" * 70 + "\n"
        "  {} — {}\n".format(label, _short(cmd)) +
        "  Started: {}\n".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) +
        "=" * 70 + "\n\n"
    )
    try:
        with open("sirocco_log", "a") as fh:
            fh.write(header)
            fh.flush()
            result = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)

        elapsed = toc(t_start, label)

        if result.returncode != 0:
            msg = "{} exited {}".format(label, result.returncode)
            die(msg)
        return result, elapsed

    except FileNotFoundError as e:
        die("{} not found: {}".format(label, e))


def run_pluto(nproc, restart_from=None, logpath="pluto_log", label="PLUTO"):
    """Run PLUTO and verify success by checking the output .dbl file exists.

    On macOS ARM, MPI often returns exit code 133 even after a perfectly clean
    run.  We therefore treat any non-zero exit as a warning rather than a
    fatal error, and instead confirm success by checking that PLUTO wrote the
    expected data.NNNN.dbl file.  If the file is missing we abort.

    Returns (result, elapsed_seconds).
    """
    cmd = pluto_cmd(nproc, restart_from=restart_from)
    log("  [{}] {}".format(label, _short(cmd)))
    t_start = tic()
    try:
        with open(logpath, "w") as fh:
            result = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    except FileNotFoundError as e:
        die("{} not found: {}".format(label, e))

    elapsed = toc(t_start, label)

    if result.returncode != 0:
        log("WARNING: {} exited {} — checking output file to confirm "
            "success (macOS MPI may return 133 on clean exit)".format(
            label, result.returncode))

    # Determine which output file number PLUTO should have produced.
    # With -restart N it writes file N+1; with no restart it writes file 1.
    expected_n = 1 if restart_from is None else restart_from + 1
    expected_dbl = "data.{:04d}.dbl".format(expected_n)

    if not os.path.exists(expected_dbl):
        die("{} failed: {} not produced. "
            "Check {} for details.".format(label, expected_dbl, logpath))

    log("  [{}] confirmed: {} exists".format(label, expected_dbl))
    return result, elapsed


def resolve_executable(*names):
    """Return the first executable found in the current directory or PATH."""
    for name in names:
        candidates = [name]
        if not name.startswith("./"):
            candidates.append("./" + name)
        for candidate in candidates:
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                return candidate
            found = shutil.which(candidate)
            if found:
                return found
    return names[0]


def check_run_dependencies():
    """Check CAK and its local datasets before any expensive cycle starts."""
    cak = os.path.abspath(data.get("CAK_EXEC", "cak_v4"))
    if not os.path.isfile(cak) or not os.access(cak, os.X_OK):
        die("CAK executable is missing or not executable: {}".format(cak))
    for filename in ("atomic_models.txt", "transitiondata.txt"):
        if not os.path.isfile(filename):
            die("Required CAK data file is missing from the run directory: " + filename)
    if shutil.which("mpirun") is None:
        die("mpirun is not on PATH; load the MPI environment before starting the driver")
    log("CAK executable, atomic datasets and MPI launcher checked")


def pluto_cmd(nproc, restart_from=None):
    """Build the mpirun ./pluto command line."""
    base = (["./pluto"] if nproc == 1
            else ["mpirun", "-np", str(nproc), "./pluto"])
    if restart_from is not None:
        base += ["-restart", str(restart_from)]
    return base


def sirocco_cmd(nproc, restart=False, photon_ramp=False):
    """Build the mpirun sirocco command line.

    Parameters
    ----------
    nproc        : int  – number of MPI ranks
    restart      : bool – pass -r to restart from an existing wind_save
    photon_ramp  : bool – pass -p 2 to ramp photons over 2 orders of magnitude.
                          ONLY used for Run A (cycle 0).  The main loop
                          always restarts (-r) with a fixed photon count.
    """
    sir_exe = resolve_executable(
        data.get("RAD_CODE", "sirocco") + data.get("RAD_CODE_VER", ""),
        "sirocco",
        "sirocco-1.1",
    )
    cmd = ["mpirun", "-np", str(nproc), sir_exe, "-f"]
    if restart:
        cmd += ["-r"]
    if photon_ramp:
        cmd += ["-p", "2"]
    cmd += ["input.pf"]
    return cmd


def rad_hydro_cmd(nproc):
    rhf_exe = resolve_executable(
        "rad_hydro_files" + data.get("RAD_CODE_VER", ""),
        "rad_hydro_files",
        "rad_hydro_files-1.1",
    )
    return ["mpirun", "-np", str(nproc), rhf_exe, "input"]


# These are regenerated every cycle by the revised SIROCCO exporter. Remove
# prior-cycle copies first so a failed exporter cannot leave stale data for CAK
# or the next PLUTO restart to consume.
RAD_HYDRO_OUTPUTS = (
    "directional_flux_r.dat", "directional_flux_theta.dat",
    "directional_flux_phi.dat", "sirocco_heatcool.dat",
    "sirocco_driving.dat", "sirocco_ion_data.dat",
    "sirocco_pcon_data.dat", "sirocco_spec_data.dat",
    "sirocco_fluxes.dat", "sirocco_debug_data.dat",
)


def run_rad_hydro(nproc, label, logpath):
    for filename in RAD_HYDRO_OUTPUTS:
        try:
            os.remove(filename)
        except FileNotFoundError:
            pass
    run(rad_hydro_cmd(nproc), logpath=logpath, check=True, label=label)
    required = ("directional_flux_r.dat", "directional_flux_theta.dat",
                "directional_flux_phi.dat", "sirocco_heatcool.dat",
                "sirocco_ion_data.dat", "sirocco_pcon_data.dat",
                "sirocco_spec_data.dat")
    missing = [name for name in required
               if not os.path.isfile(name) or os.path.getsize(name) == 0]
    if missing:
        die("rad_hydro_files did not produce required current-cycle outputs: "
            + ", ".join(missing))


# ---------------------------------------------------------------------------
# PER-CYCLE IONISATION-STATE ARCHIVE  (independent of SPEC_RUN)
#
# `directory` (cycleNNNNNNNN/) has just been fully populated by the normal
# main-loop archive step, and still contains this cycle's input.wind_save.
# That folder itself is deleted wholesale once it is more than 5 cycles old
# (see "old_dir" pruning below), so without this step the ionisation state
# of any earlier cycle would be lost for good.
#
# This function runs windsave2table on that wind_save (producing lightweight
# ascii ion tables), then copies everything now sitting in `directory` --
# those new ion tables plus the same static/sirocco_* files already archived
# there -- into ION_TABLE_DIR/cycleNNNNNNNN/, EXCLUDING the wind_save binary
# itself. Unlike cycleNNNNNNNN/, folders under ION_TABLE_DIR are never
# pruned -- they accumulate for the whole run, independent of whether this
# cycle was also selected for a full spectral-only test run.
#
# Controlled by ION_TABLE_ALL_CYCLES ("yes"/"no", default "yes"). Failures
# are logged as warnings and never abort the campaign.
# ---------------------------------------------------------------------------

def archive_ion_table(i, directory):
    """Run windsave2table on directory/input.wind_save, then copy everything
    in `directory` except the wind_save binary and the files that are
    already kept, unpruned, in the main run directory (data.NNNN.dbl and
    grid.out) into ION_TABLE_DIR/<basename of directory>/.
    """
    if str(data.get("ION_TABLE_ALL_CYCLES", "yes")).lower() != "yes":
        return

    ws_path = os.path.join(directory, "input.wind_save")
    if not os.path.exists(ws_path):
        log("  WARNING: ION_TABLE_ALL_CYCLES set but {} not found -- "
            "skipping ion-table archive for cycle {}".format(ws_path, i))
        return

    w2t_exe = resolve_executable(
        "windsave2table" + data.get("RAD_CODE_VER", ""),
        "windsave2table",
        "windsave2table-1.1",
    )
    w2t_exe = os.path.abspath(w2t_exe) if os.path.exists(w2t_exe) else w2t_exe

    cmd_w2t = [w2t_exe, "input"]
    log("  [ION-TABLE-{}] {}".format(i, _short(cmd_w2t)))
    logpath_w2t = os.path.join(directory, "input_windsave2table_log")
    try:
        with open(logpath_w2t, "w") as fh:
            res_w2t = subprocess.run(cmd_w2t, stdout=fh,
                                     stderr=subprocess.STDOUT,
                                     cwd=directory)
        if res_w2t.returncode != 0:
            log("  WARNING: windsave2table on cycle {} exited {} -- see {}".format(
                i, res_w2t.returncode, logpath_w2t))
            return
        log("  windsave2table -> input.master.txt (+ per-element ion tables)")
    except FileNotFoundError as e:
        log("  WARNING: windsave2table could not start ({}) -- skipping ion "
            "table for cycle {}".format(e, i))
        return
    except Exception as e:
        log("  WARNING: windsave2table failed on cycle {}: {}".format(i, e))
        return

    ion_root_dir = data.get("ION_TABLE_DIR", "ions_data")
    ion_dir = os.path.join(ion_root_dir, os.path.basename(directory))
    os.makedirs(ion_dir, exist_ok=True)

    # AM - 2609 - data.NNNN.dbl is PLUTO's own per-cycle snapshot and is
    # never pruned from the main run directory (every cycle's .dbl stays
    # there for the life of the run), and grid.out is identical, unchanging
    # geometry metadata written once and copied unchanged into every
    # cycleNNNNNNNN/ folder. Copying either into ions_data/ as well would
    # just duplicate data already available elsewhere for no benefit --
    # a later plotting script can simply read data.NNNN.dbl and grid.out
    # straight from the main run directory by cycle number.
    dbl_name = "data.{:04d}.dbl".format(i)
    exclude = {"input.wind_save", "grid.out", dbl_name}

    n_copied = 0
    for fn in os.listdir(directory):
        if fn in exclude:
            continue
        src = os.path.join(directory, fn)
        if not os.path.isfile(src):
            continue
        shutil.copy2(src, os.path.join(ion_dir, fn))
        n_copied += 1

    log("  Ion-state archive: {} -> {} ({} files; {} excluded)".format(
        directory, ion_dir, n_copied, ", ".join(sorted(exclude))))


# ---------------------------------------------------------------------------
# SPECTRAL-ONLY TEST RUN
#
# Called from the main loop immediately after SIROCCO has finished, while
# input.wind_save still holds the ionisation state converged for THIS RHD
# snapshot.  That wind_save is copied aside and a spectral-cycles-only run
# (Ionization_cycles 0) is launched against it, so the resulting line profiles
# can be compared with those from the normal ionisation+spectral approach.
#
# The run happens inline: the main loop waits for it, then continues to the
# next cycle as usual.
#
# CRITICAL — cwd isolation:
#   SIROCCO writes sirocco_heatcool.dat, sirocco_driving.dat, sirocco_ion_data.dat,
#   sirocco_spec_data.dat and a .wind_save into its WORKING DIRECTORY.  If the
#   spectral run were launched in the main run directory it would overwrite
#   the files the live cycle has just produced, and pre_calc(ifile) would then
#   build prefactors.dat from the spectral run instead of the ionisation run —
#   silently corrupting the Blondin prefactors fed back into PLUTO for every
#   subsequent cycle.  The run is therefore executed with cwd=spec_dir.
#
# Failures are logged as warnings and never abort the campaign.
# ---------------------------------------------------------------------------

def run_spectral_cycles(i, sir_model, nproc):
    """Build spec_data/cycleNNNN/ and run spectral cycles from the wind_save.

    Parameters
    ----------
    i         : int  – current main-loop cycle number
    sir_model : str  – the .pluto model file for this cycle (provenance copy)
    nproc     : int  – MPI ranks for the spectral run (data["nproc_spec"])

    Returns True if the spectral run completed (or was only prepared),
    False if anything went wrong.
    """
    spec_root_dir = data.get("SPEC_DIR", "spec_data")
    mode          = str(data.get("SPEC_MODE", "run")).lower()
    n_spec_cycles = data.get("SPEC_CYCLES", 30)

    cyc_name      = "cycle{:04d}".format(i)          # wind_save root
    pf_root       = "spec{:04d}".format(i)           # parameter-file root
    spec_dir      = os.path.join(spec_root_dir, cyc_name)

    log("-" * 70)
    log("SPECTRAL TEST: cycle {} -> {}".format(i, spec_dir))
    log("-" * 70)

    t_spec = tic()

    try:
        os.makedirs(spec_dir, exist_ok=True)
    except Exception as e:
        log("WARNING: cannot create {}: {} — spectral run skipped".format(spec_dir, e))
        return False

    # ── atomic data: give the spec folder the same data/ and zdata/ that ─────
    # Setup_Sirocco_Dir creates in the main run directory.
    #
    # SIROCCO resolves the Atomic_data path (data/standard80.dat) relatively,
    # against its working directory — and since the spectral run executes with
    # cwd=spec_dir, that folder needs its own data/ and zdata/ just as the run
    # directory has.  Preferred route is to run Setup_Sirocco_Dir there, which
    # is exactly what the main directory did; if it is not on PATH we fall back
    # to symlinking whatever the run directory already resolved them to.
    setup_exe = shutil.which("Setup_Sirocco_Dir")
    setup_ok  = False

    if setup_exe:
        try:
            res = subprocess.run([setup_exe], cwd=spec_dir,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True)
            if res.returncode == 0:
                setup_ok = True
                log("  Setup_Sirocco_Dir run in {}".format(spec_dir))
            else:
                log("  WARNING: Setup_Sirocco_Dir exited {} — falling back to "
                    "symlinks".format(res.returncode))
        except Exception as e:
            log("  WARNING: Setup_Sirocco_Dir failed ({}) — falling back to "
                "symlinks".format(e))

    # Fall back to (or top up) symlinks for anything still missing.
    for sub in ("data", "zdata"):
        dst = os.path.join(spec_dir, sub)
        if os.path.exists(dst):
            continue
        try:
            if os.path.exists(sub):
                # mirror whatever the run directory points at
                os.symlink(os.path.realpath(sub), dst)
                log("  linked {} -> {}".format(sub, os.path.realpath(sub)))
            elif os.environ.get("SIROCCO"):
                target = os.path.join(os.environ["SIROCCO"], sub)
                if os.path.exists(target):
                    os.symlink(target, dst)
                    log("  linked {} -> {}".format(sub, target))
                else:
                    log("  WARNING: {} not found in $SIROCCO".format(sub))
            elif not setup_ok:
                log("  WARNING: no {} in run dir and $SIROCCO unset — "
                    "spectral run may not find atomic data".format(sub))
        except Exception as e:
            log("  WARNING: cannot create {} link: {}".format(sub, e))

    # ── copy the converged wind_save (the reference for this snapshot) ────────
    ws_src = "input.wind_save"
    ws_dst = os.path.join(spec_dir, "{}.wind_save".format(cyc_name))
    if not os.path.exists(ws_src):
        log("WARNING: {} not found — spectral run skipped".format(ws_src))
        return False
    try:
        shutil.copy(ws_src, ws_dst)
        log("  wind_save -> {}".format(os.path.basename(ws_dst)))
    except Exception as e:
        log("WARNING: cannot copy wind_save: {} — spectral run skipped".format(e))
        return False

    # ── copy the .pluto model file (provenance only; not read in 'previous') ──
    try:
        if sir_model and os.path.exists(sir_model):
            shutil.copy(sir_model, os.path.join(spec_dir, os.path.basename(sir_model)))
            log("  model file -> {} (provenance only)".format(os.path.basename(sir_model)))
    except Exception as e:
        log("  WARNING: cannot copy {}: {}".format(sir_model, e))

    # ── write the spectral .pf ────────────────────────────────────────────────
    # NOTE: pf_root ('specNNNN') deliberately differs from cyc_name
    # ('cycleNNNN'), because SIROCCO names its outputs after the .pf root.
    # Identical names would make the run overwrite its own input wind_save.
    cwd_save = os.getcwd()
    try:
        os.chdir(spec_dir)
        pss.sirocco_spec_input_file(pf_root, data,
                                    windsave_root=cyc_name,
                                    cycles=n_spec_cycles)
        log("  wrote {}.pf  (Ionization_cycles 0, Spectrum_cycles {})".format(
            pf_root, n_spec_cycles))
    except Exception as e:
        os.chdir(cwd_save)
        log("WARNING: cannot write spectral .pf: {} — spectral run skipped".format(e))
        return False
    finally:
        os.chdir(cwd_save)

    if mode != "run":
        toc(t_spec, "spectral-prepare-{}".format(i))
        log("  SPEC_MODE = '{}' — folder prepared, run not launched".format(mode))
        return True

    # ── run SIROCCO, spectral cycles only, inside the spec folder ────────────
    spec_exe = resolve_executable(
        data.get("RAD_CODE", "sirocco") + data.get("RAD_CODE_VER", ""),
        "sirocco",
        "sirocco-1.1",
    )
    spec_exe = os.path.abspath(spec_exe) if os.path.exists(spec_exe) else spec_exe

    cmd = ["mpirun", "-np", str(nproc), spec_exe, "-f", "{}.pf".format(pf_root)]
    log("  [SPEC-{}] {}".format(i, _short(cmd)))

    logpath = os.path.join(spec_dir, "sirocco_spec_log")
    try:
        with open(logpath, "w") as fh:
            fh.write("=" * 70 + "\n")
            fh.write("  Spectral-only run, cycle {}\n".format(i))
            fh.write("  wind_save : {}.wind_save\n".format(cyc_name))
            fh.write("  Started   : {}\n".format(
                datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            fh.write("=" * 70 + "\n\n")
            fh.flush()
            result = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                                    cwd=spec_dir)
    except FileNotFoundError as e:
        log("WARNING: spectral run could not start ({}) — continuing".format(e))
        return False
    except Exception as e:
        log("WARNING: spectral run failed ({}) — continuing".format(e))
        return False

    elapsed = toc(t_spec, "SPEC-{}".format(i))

    if result.returncode != 0:
        log("WARNING: spectral run for cycle {} exited {} — see {}".format(
            i, result.returncode, logpath))
        return False

    spec_out = os.path.join(spec_dir, "{}.spec".format(pf_root))
    if os.path.exists(spec_out):
        log("  spectral run complete: {} ({:.1f} min)".format(
            os.path.basename(spec_out), elapsed / 60.0))
    else:
        log("  WARNING: spectral run exited cleanly but {}.spec not found".format(pf_root))

    # ── windsave2table: ionisation-state tables from the wind_save ──────
    #
    # wind_save binaries are large; windsave2table converts one into a set
    # of lightweight ascii tables — one master table plus one ion-fraction
    # table per element — that are cheap to keep long-term as the durable
    # record of the ionisation state, once the wind_save itself is pruned
    # (see the retention step below). Controlled by SPEC_WINDSAVE2TABLE
    # ("yes"/"no", default "yes").
    #
    # Run on both wind_save files sitting in spec_dir:
    #   cycleNNNN  — the RHD-snapshot ionisation state that seeded this run
    #   specNNNN   — SIROCCO's own output for this run (same ionisation
    #                state, since Ionization_cycles=0, but written under the
    #                .pf root SIROCCO itself considers canonical)
    if str(data.get("SPEC_WINDSAVE2TABLE", "yes")).lower() == "yes":
        w2t_exe = resolve_executable(
            "windsave2table" + data.get("RAD_CODE_VER", ""),
            "windsave2table",
            "windsave2table-1.1",
        )
        w2t_exe = os.path.abspath(w2t_exe) if os.path.exists(w2t_exe) else w2t_exe

        for root in (cyc_name, pf_root):
            ws_path = os.path.join(spec_dir, "{}.wind_save".format(root))
            if not os.path.exists(ws_path):
                continue
            cmd_w2t = [w2t_exe, root]
            log("  [W2T-{}] {}".format(i, _short(cmd_w2t)))
            logpath_w2t = os.path.join(spec_dir, "{}_windsave2table_log".format(root))
            try:
                with open(logpath_w2t, "w") as fh:
                    res_w2t = subprocess.run(cmd_w2t, stdout=fh,
                                             stderr=subprocess.STDOUT,
                                             cwd=spec_dir)
                if res_w2t.returncode != 0:
                    log("  WARNING: windsave2table on {} exited {} — see {}".format(
                        root, res_w2t.returncode, logpath_w2t))
                else:
                    log("  windsave2table -> {}.master.txt (+ per-element ion tables)".format(root))
            except FileNotFoundError as e:
                log("  WARNING: windsave2table could not start ({}) — skipping ion "
                    "tables for this cycle".format(e))
                break   # binary genuinely missing/unresolved; retrying for the
                        # other root would just fail the same way
            except Exception as e:
                log("  WARNING: windsave2table failed on {}: {}".format(root, e))

    # ── Rolling wind_save retention for spec_data/ ───────────────────────
    #
    # SPEC_KEEP_WINDSAVE controls how many of the most-recent spectral
    # snapshots keep their *.wind_save binaries:
    #   "all"   -> never delete any wind_save (disk-hungry, full fidelity)
    #   N (int) -> keep the wind_save files for the N most recently
    #              selected spectral cycles (default 5); older
    #              spec_data/cycleNNNN folders keep their .spec, .pf, log
    #              and windsave2table ion tables, but have their
    #              *.wind_save removed, since those tables are now the
    #              durable record of the ionisation state.
    #
    # Selected spectral cycles are equally spaced by SPEC_STEP starting at
    # SPEC_START, so "N cycles ago" is simply i - N*SPEC_STEP.
    keep = data.get("SPEC_KEEP_WINDSAVE", 5)
    if str(keep).lower() != "all":
        try:
            n_keep = int(keep)
        except (TypeError, ValueError):
            log("  WARNING: SPEC_KEEP_WINDSAVE={!r} not 'all' or an int — "
                "defaulting to 5".format(keep))
            n_keep = 5

        step = data.get("SPEC_STEP", 1)
        if n_keep > 0 and step > 0:
            old_i = i - n_keep * step
            if old_i >= data.get("SPEC_START", 0):
                old_dir = os.path.join(spec_root_dir, "cycle{:04d}".format(old_i))
                removed = []
                for root in ("cycle{:04d}".format(old_i), "spec{:04d}".format(old_i)):
                    old_ws = os.path.join(old_dir, "{}.wind_save".format(root))
                    if os.path.exists(old_ws):
                        try:
                            os.remove(old_ws)
                            removed.append(os.path.basename(old_ws))
                        except Exception as e:
                            log("  WARNING: cannot remove {}: {}".format(old_ws, e))
                if removed:
                    log("  Pruned wind_save from {} (kept .spec/.pf/ion tables): {}"
                        .format(old_dir, ", ".join(removed)))

    return True


check_run_dependencies()

# ---------------------------------------------------------------------------
# Code units
# ---------------------------------------------------------------------------

try:
    UNIT_DENSITY, UNIT_LENGTH, UNIT_VELOCITY = pss.get_units()
    log("Code units read from definitions.h")
except Exception as e:
    die("Cannot read definitions.h: {}".format(e))

# # ---------------------------------------------------------------------------

nproc_pluto    = data.get("nproc_pluto",    1)
nproc_sir      = data.get("nproc_sir",      1)
nproc_cak      = data.get("nproc_cak",      1)
nproc_radhydro = data.get("nproc_radhydro", nproc_sir)
nproc_spec     = data.get("nproc_spec",     nproc_sir)

sir_exe = data.get("RAD_CODE", "sirocco") + data.get("RAD_CODE_VER", "")

# ---------------------------------------------------------------------------
# Restart / fresh-start logic
# ---------------------------------------------------------------------------

if len(sys.argv) < 2:
    istart   = 0
    sim_time = dt_rad  # first main-loop PLUTO run advances to dt_rad
    log("Fresh start: cycle 0")
else:
    try:
        istart = int(sys.argv[1])
    except ValueError:
        die("Invalid restart argument: {}".format(sys.argv[1]))

    if istart > 0:
        cyc_dir = "cycle{:08d}".format(istart)
        log("Restarting from {}".format(cyc_dir))
        try:
            for fn in os.listdir(cyc_dir):
                src = os.path.join(cyc_dir, fn)
                if os.path.isfile(src):
                    shutil.copy(src, ".")
        except Exception as e:
            die("Cannot restore from {}: {}".format(cyc_dir, e))
        try:
            last     = pss.pload(istart)
            sim_time = last.SimTime + dt_rad
            log("Restarting at t = {:.6f}".format(sim_time))
        except Exception as e:
            die("Cannot load PLUTO restart data: {}".format(e))
        # cycle N is complete and archived: continue with cycle N+1
        # (PLUTO -restart N -> data.{N+1}.dbl), do not redo cycle N.
        istart += 1
        log("Continuing with cycle {}".format(istart))
    else:
        istart   = 0
        sim_time = dt_rad
        log("Restart index <= 0 -> treating as fresh start")

# ---------------------------------------------------------------------------
# CYCLE 0 — INITIALISATION (Run A)
#
# Only executed on a fresh start (istart == 0).
#
#  1. Generate data.0000.dbl if absent (short PLUTO run at t=0).
#  2. Convert data.0000.dbl -> 00000000.pluto (SIROCCO model).
#  3. Run SIROCCO 30 cycles with -p 2  (fresh start, photon ramp).
#  4. Run rad_hydro_files -> directional_flux_*.dat.
#  5. Save flux files as *_init.dat (restart safety copies).
#  6. Run CAK from data.0000.dbl.
#
# Blondin prefactors are intentionally NOT computed in cycle 0.
# pre_calc(ifile) divides SIROCCO heating/cooling rates by the Blondin
# rate coefficients stored as uservars in data.{ifile}.dbl.  For the
# result to be physically meaningful those coefficients must be positive.
# data.0000.dbl is the raw PLUTO initial snapshot; its ch, cc, lc, bc, xh
# uservars may be zero or NaN for shadowed inner-disc cells (zero X-ray
# heating) or if the file was produced with stale flux/prefactor files from
# a prior run.  Dividing by zero or NaN produces NaN prefactors that then
# propagate into every subsequent cycle through the prefactors.dat -> C read
# -> userdef_output write -> pload chain.
#
# Instead, the main loop cycle i=1 is the first pre_calc call.  At that
# point pload(1) reads data.0001.dbl, which was written by PLUTO restarting
# from data.0000.dbl with:
#   - NFLUX_ANGLES > 0   (flux files from Run A are present)
#   - comp_h_pre = 1.0   (no prior prefactors.dat -> read_sirocco_heatcool
#                          initialises all *_pre arrays to unity)
#   - sirocco_xi > 0     (read from sirocco_heatcool.dat produced by Run A)
# So D.ch = 8.9e-36 * xi_sirocco * t_r > 0 everywhere, and ideal_prefactor
# is a well-defined physical ratio.  This is exactly what Version 1 did.
#
# After cycle 0 lc_NFLUX_ANGLES will be > 0 in all subsequent PLUTO runs.
# The main loop begins at i=1 with sim_time = dt_rad.
# ---------------------------------------------------------------------------

if istart == 0:

    t_cyc0 = tic()

    log("=" * 70)
    log("CYCLE 0 INITIALISATION")
    log("=" * 70)

    # ------------------------------------------------------------------
    # 1. Ensure data.0000.dbl exists
    # ------------------------------------------------------------------

    log("")
    log("--- Cycle 0: building SIROCCO model from data.0000.dbl ---")

    if not os.path.exists("data.0000.dbl"):
        log("data.0000.dbl not found — running PLUTO once to generate initial conditions")
        pss.pluto_input_file(dt_rad, data)
        # Initialize and write snapshot zero without evolving unilluminated gas.
        # A failed MPI/PLUTO launch must stop the driver even if it wrote data.0000.
        run(pluto_cmd(nproc_pluto) + ["-maxsteps", "0"],
            logpath="pluto_log_init", check=True, label="PLUTO-init")
        if not os.path.exists("data.0000.dbl"):
            die("data.0000.dbl still not produced — check pluto_log_init.")
        log("  data.0000.dbl generated successfully")

    # ------------------------------------------------------------------
    # 2. Convert data.0000.dbl -> SIROCCO model
    # ------------------------------------------------------------------

    log("Converting data.0000.dbl -> SIROCCO model (ifile=0 -> 00000000.pluto)")
    try:
        sir_model_init = pss.pluto2sir_rtheta(0)
    except Exception as e:
        die("pluto2sir_rtheta(0) failed: {}".format(e))
    log("SIROCCO model: {}".format(sir_model_init))

    root_init = "00000000"
    try:
        pss.sirocco_input_file(root_init, data, cycles=init_sir_cycles,
                               model_file=sir_model_init)
    except Exception as e:
        die("sirocco_input_file failed (init): {}".format(e))
    shutil.copy("{}.pf".format(root_init), "input.pf")

    # ------------------------------------------------------------------
    # 3. SIROCCO Run A — 30 cycles, -p 2 (photon ramp, fresh start)
    # ------------------------------------------------------------------
    #

    log("--- Cycle 0: SIROCCO ({} cycles) ---".format(init_sir_cycles))
    sirocco_run(sirocco_cmd(nproc_sir, restart=False, photon_ramp=False),
                label="SIROCCO-0A")

    # ------------------------------------------------------------------
    # 4. rad_hydro_files -> directional flux files
    # ------------------------------------------------------------------

    log("--- Cycle 0: rad_hydro_files ---")
    run_rad_hydro(nproc_radhydro, "rad_hydro_files-0A",
                  "rad_hydro_files_output_0A")

    # ------------------------------------------------------------------
    # 5. Save flux files as *_init.dat
    # ------------------------------------------------------------------

    flux_bases = ("directional_flux_r",
                  "directional_flux_theta",
                  "directional_flux_phi")

    for base in flux_bases:
        src = "{}.dat".format(base)
        dst = "{}_init.dat".format(base)
        if os.path.exists(src):
            shutil.copy(src, dst)
            log("Saved {}".format(dst))
        else:
            die("{} not produced by rad_hydro_files — cannot continue.\n"
                "    Check rad_hydro_files_output_0A for errors.".format(src))

    log("Directional flux files ready.")

    # ------------------------------------------------------------------
    # 6. CAK from data.0000.dbl
    # ------------------------------------------------------------------

    log("--- Cycle 0: CAK ---")
    run(["mpirun", "-np", str(nproc_cak),
         os.path.abspath(data.get("CAK_EXEC", "cak_v4"))],
        logpath="cak_output_0", check=True, label="CAK-0")

    # ------------------------------------------------------------------
    # 7. Archive sirocco_* outputs from SIROCCO Run A with 00000000_ prefix.
    #    (pre_calc is not called here — see cycle-0 comment block above.)
    # ------------------------------------------------------------------

    for src, dst in [
        ("sirocco_heatcool.dat",  "00000000_sirocco_heatcool.dat"),
        ("sirocco_driving.dat",   "00000000_sirocco_driving.dat"),
        ("sirocco_pcon_data.dat", "00000000_sirocco_pcon_data.dat"),
        ("sirocco_ion_data.dat",  "00000000_sirocco_ion_data.dat"),
        ("sirocco_spec_data.dat", "00000000_sirocco_spec_data.dat"),
        ("sirocco_fluxes.dat",    "00000000_sirocco_fluxes.dat"),
    ]:
        if os.path.exists(src):
            shutil.copy(src, dst)
        else:
            log("  WARNING: {} not found — skipping".format(src))

    # ------------------------------------------------------------------
    # 8. Archive cycle 0 -> cycle00000000/
    # ------------------------------------------------------------------

    dir0 = "cycle00000000"
    if os.path.exists(dir0):
        backup = dir0 + "_old"
        if os.path.exists(backup):
            shutil.rmtree(backup)
        shutil.move(dir0, backup)
    os.mkdir(dir0)

    cycle0_files = [
        "data.0000.dbl", "dbl.out", "grid.out", "restart.out",
        "vtk.out", "flt.out", "dbl.h5.out", "flt.h5.out",
        "pluto.ini", "definitions.h",
        "00000000.pluto", "00000000.pf",
        "directional_flux_r.dat",
        "directional_flux_theta.dat",
        "directional_flux_phi.dat",
        "directional_flux_r_init.dat",
        "directional_flux_theta_init.dat",
        "directional_flux_phi_init.dat",
        "M_UV_data.dat", "M_UV_data.meta",
        "input.wind_save", "sirocco_log",
        "00000000_sirocco_heatcool.dat",
        "00000000_sirocco_driving.dat",
        "00000000_sirocco_pcon_data.dat",
        "00000000_sirocco_ion_data.dat",
        "00000000_sirocco_spec_data.dat",
        "00000000_sirocco_fluxes.dat",
    ]
    for fn in cycle0_files:
        if os.path.exists(fn):
            dst = os.path.join(dir0, os.path.basename(fn))
            if os.path.exists(dst):
                os.remove(dst)
            shutil.copy(fn, dir0)
    for fn in glob.glob("sirocco_*.dat"):
        dst = os.path.join(dir0, os.path.basename(fn))
        if os.path.exists(dst):
            os.remove(dst)
        shutil.copy(fn, dir0)

    archive_ion_table(0, dir0)

    toc(t_cyc0, "Cycle 0 total")
    log("Cycle 0 archived to {}".format(dir0))
    log("Initialisation complete — main loop will start PLUTO from data.0000.dbl")
    log("")

    # Main loop starts at i=1; sim_time is already set to dt_rad above
    istart = 1

# ---------------------------------------------------------------------------
# Ensure flux files are present before the main loop
# (needed on restart; on fresh start they were created in cycle 0 above)
# ---------------------------------------------------------------------------

for base in ("directional_flux_r",
             "directional_flux_theta",
             "directional_flux_phi"):
    if not os.path.exists("{}.dat".format(base)):
        init_copy = "{}_init.dat".format(base)
        if os.path.exists(init_copy):
            shutil.copy(init_copy, "{}.dat".format(base))
            log("Restored {}.dat from {}".format(base, init_copy))
        else:
            die("Neither {base}.dat nor {base}_init.dat found.\n"
                "    Re-run from cycle 0.".format(base=base))

# ---------------------------------------------------------------------------
# Main PLUTO <-> SIROCCO <-> CAK loop
#
# Loop variable  i : cycle number (starts at istart, which is 1 on fresh run)
# restart_file     : PLUTO restarts from data.{restart_file:04d}.dbl  = i-1
# ifile            : PLUTO produces data.{ifile:04d}.dbl               = i
# archive dir      : cycle{ifile:08d}
#
# SIROCCO total cycles = init_sir_cycles + i * sir_cycles
#   i=1 -> 30+2=32   (continues from Run-A wind_save via -r, no -p 2)
#   i=2 -> 30+4=34
#   …
# ---------------------------------------------------------------------------

for i in range(istart, tot_num_cycle):

    t_cyc = tic()

    log("=" * 70)
    log("MAIN LOOP: cycle {}   t = {:.6f}".format(i, sim_time))
    log("=" * 70)

    restart_file     = i - 1
    ifile            = i
    root             = "{:08d}".format(ifile)
    dbl              = "data.{:04d}.dbl".format(ifile)
    directory        = "cycle" + root

    # SIROCCO: always restart from existing wind_save (no -p 2)
    total_sir_cycles = init_sir_cycles + i * sir_cycles

    # ── PLUTO ──────────────────────────────────────────────────────────

    pss.pluto_input_file(sim_time, data)

    run_pluto(nproc_pluto, restart_from=restart_file,
              logpath="pluto_log", label="PLUTO-{}".format(i))

    sim_time += dt_rad

    # ── Convert PLUTO -> SIROCCO model ────────────────────────────────

    log("Converting {} -> SIROCCO model".format(dbl))
    try:
        sir_model = pss.pluto2sir_rtheta(ifile)
    except Exception as e:
        die("pluto2sir_rtheta failed at cycle {}: {}".format(i, e))

    # ── hydro2sirocco then SIROCCO (-r, no -p 2) ───────────────────────

    pss.sirocco_input_file(root, data,
                           cycles=total_sir_cycles,
                           model_file=sir_model)
    shutil.copy("{}.pf".format(root), "input.pf")

    run([resolve_executable("hydro2sirocco" + data.get("RAD_CODE_VER", ""),
                            "hydro2sirocco", "hydro2sirocco-1.1"),
         "-model_file", sir_model, "input"],
        logpath="mod_wind_log", check=True,
        label="hydro2sirocco-{}".format(i))

    shutil.copy("new.wind_save", "input.wind_save")

    log("--- Cycle {}: SIROCCO ({} cycles, -r) ---".format(i, total_sir_cycles))
    sirocco_run(sirocco_cmd(nproc_sir, restart=True, photon_ramp=False),
                label="SIROCCO-{}".format(i))

    # ── Spectral-only test run (if this cycle is selected) ────────────
    #
    # Runs here, immediately after SIROCCO, because input.wind_save now
    # holds the ionisation state converged for THIS RHD snapshot.  Must
    # come before rad_hydro_files / pre_calc, and is executed with
    # cwd=spec_data/cycleNNNN so it cannot overwrite the live sirocco_*.dat
    # files that pre_calc depends on.

    if pss.spec_cycle_selected(i, data):
        run_spectral_cycles(i, sir_model, nproc_spec)

    # ── rad_hydro_files ───────────────────────────────────────────────

    run_rad_hydro(nproc_radhydro, "rad_hydro_files-{}".format(i),
                  "rad_hydro_files_output")

    # ── Archive sirocco_* outputs ──────────────────────────────────────────

    for src, dst in [
        ("sirocco_heatcool.dat",  "{}_sirocco_heatcool.dat".format(root)),
        ("sirocco_driving.dat",   "{}_sirocco_driving.dat".format(root)),
        ("sirocco_pcon_data.dat", "{}_sirocco_pcon_data.dat".format(root)),
        ("sirocco_ion_data.dat",  "{}_sirocco_ion_data.dat".format(root)),
        ("sirocco_spec_data.dat", "{}_sirocco_spec_data.dat".format(root)),
        ("sirocco_fluxes.dat",    "{}_sirocco_fluxes.dat".format(root)),
    ]:
        if os.path.exists(src):
            if os.path.exists(dst): os.remove(dst)
            shutil.copy(src, dst)

    # ── CAK ───────────────────────────────────────────────────────────

    run(["mpirun", "-np", str(nproc_cak),
         os.path.abspath(data.get("CAK_EXEC", "cak_v4"))],
        logpath="cak_output", check=True,
        label="CAK-{}".format(i))

    # ── Blondin prefactors ────────────────────────────────────────────

    log("Computing Blondin prefactors")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pss.pre_calc(ifile)
    except Exception as e:
        die("pre_calc failed at cycle {}: {}".format(i, e))

    # ── Archive cycle directory ───────────────────────────────────────

    old_dir = "cycle{:08d}".format(ifile - 5)
    if ifile > 5 and os.path.exists(old_dir):
        try:
            shutil.rmtree(old_dir)
            log("Removed old archive {}".format(old_dir))
        except Exception as e:
            log("WARNING: cannot remove {}: {}".format(old_dir, e))

    if os.path.exists(directory):
        backup = directory + "_old"
        if os.path.exists(backup): shutil.rmtree(backup)
        shutil.move(directory, backup)
    os.mkdir(directory)

    static_files = [
        "dbl.out", "pluto.ini", "restart.out", "grid.out",
        "vtk.out", "flt.out", "dbl.h5.out", "flt.h5.out",
        "definitions.h", dbl,
        "M_UV_data.dat", "M_UV_data.meta", "prefactors.dat",
        "directional_flux_r.dat",
        "directional_flux_theta.dat",
        "directional_flux_phi.dat",
        "input.wind_save", "sirocco_log",
    ]
    for fn in static_files:
        if os.path.exists(fn):
            dst = os.path.join(directory, os.path.basename(fn))
            if os.path.exists(dst): os.remove(dst)
            shutil.copy(fn, directory)

    for fn in glob.glob("sirocco_*.dat"):
        dst = os.path.join(directory, os.path.basename(fn))
        if os.path.exists(dst): os.remove(dst)
        shutil.copy(fn, directory)

    for fn in [sir_model, "{}.pf".format(root),
               "{}_sirocco_heatcool.dat".format(root),
               "{}_sirocco_driving.dat".format(root),
               "{}_sirocco_pcon_data.dat".format(root),
               "{}_sirocco_ion_data.dat".format(root),
               "{}_sirocco_spec_data.dat".format(root),
               "{}_sirocco_fluxes.dat".format(root)]:
        if os.path.exists(fn):
            dst = os.path.join(directory, os.path.basename(fn))
            if os.path.exists(dst): os.remove(dst)
            shutil.move(fn, directory)

    archive_ion_table(i, directory)

    elapsed = toc(t_cyc, "")
    log("Cycle {} archived to {}  [total: {:.1f} s / {:.2f} min]".format(
        i, directory, elapsed, elapsed / 60.0))

logfile.close()
print("Fin")
