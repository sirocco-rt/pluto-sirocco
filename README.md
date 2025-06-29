# PLUTO-SIROCCO v4.4.3

**Monte Carlo Radiation Hydrodynamics (MC-RHD) of Line-Driven Disk Winds**

This repository provides a working version of PLUTO v4.4.3 modified to couple with the SIROCCO radiative transfer code and an external CAK solver to simulate line-driven winds from accretion disks.

---

## 🧩 Code Components

- **PLUTO v4.4.3** — Handles hydrodynamics  
- **SIROCCO** — Performs Monte Carlo radiative transfer  
- **CAK** — Computes force multipliers  
- **Python Scripts** — Coordinate sequential PLUTO-SIROCCO-CAK execution  

---

## 🔧 Installation Guide

### 1. Install SIROCCO

```bash
git clone https://github.com/sirocco-rt/sirocco.git
cd sirocco/source
# Replace the rad_hydro_files.c with the custom version
```

**Update your environment:**
```bash
export SIROCCO=/path/to/sirocco
export PATH=$SIROCCO/bin:$PATH
source ~/.bashrc  # or ~/.zshrc
```

**Build SIROCCO:**
```bash
cd $SIROCCO
./configure
make install 2>&1 | tee install_log.txt
make clean
```

---

### 2. Install PLUTO-SIROCCO v4.4.3

```bash
git clone https://github.com/aminmosallanezhad/PLUTO-SIROCCO.git
# Unzip if downloaded manually, and rename the folder to PLUTO
```

**Update your environment:**
```bash
export PLUTO_DIR=/path/to/PLUTO
source ~/.bashrc  # or ~/.zshrc
```

---

### 3. Install pyPLUTO (custom version)

```bash
cd $PLUTO_DIR/Tools/pyPLUTO
python3 setup.py install  # Or with sudo if needed
```

> ⚠️ Do **not** use pyPLUTO from other versions. This one includes required changes to `pload.py`.  
[More info](https://groups.google.com/g/pluto_users/c/7OziFlHLdU0/m/j4robn-rAwAJ)

---

### 4. Install CAK

```bash
cd $PLUTO_DIR/Test_Problems/LineDrivenWind/cak_v3
make
```

> Copy the `cak_v3` executable (only) to your working directory. Do **not** run the makefile outside the `cak_v3` folder.

---

## 🧪 Running the Code

```bash
cd $PLUTO_DIR/Test_Problems/LineDrivenWind
cp -r cv_idl ~/my_sim_dir
cd ~/my_sim_dir
python3 $PLUTO_DIR/setup.py
```

**Configuration highlights:**
- EOS: IDEAL or ISOTHERMAL  
- COOLING: BLONDIN (for IDEAL)  
- BODY_FORCE: VECTOR  
- LINE_DRIVEN_WIND: SIROCCO_MODE

Then:
```bash
make
make clean
```

**Review config file:**
```bash
nano pluto_sirocco_config.py
```

**Run simulation:**
```bash
python3 ./pluto_sirocco_dir_iso.py
```

Monitor output:
```bash
tail -f sirocco_log
```

If config is changed:
```bash
python3 ./pluto_sirocco_init.py
```

---

## ⚙️ Code Structure Summary

- `cak_v3/` — `atomic_models.txt`, `transitiondata.txt`, `cak_v3` executable  
- `pluto_sirocco_*.py` — Python drivers and control scripts  
- `init.c`, `line_connect.c`, `cooling.c`, etc. — Custom PLUTO source code

---

## 📡 Running on Iridis HPC

**Example Slurm Script (`job.slurm`):**
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=192
#SBATCH --time=48:00:00
#SBATCH --partition=highmem
#SBATCH --mem=2500000M
#SBATCH --job-name=idl_macro
#SBATCH --output=idl_output.out
#SBATCH --error=idl_error.err

module load intel-mpi/2021.14
export PLUTO_DIR=/home/am1f24/PLUTO
export SIROCCO=/home/am1f24/sirocco
export PATH=$SIROCCO/bin:$PATH
source /iridisfs/i6software/conda/miniconda-py3/etc/profile.d/conda.sh
conda activate myenv39

python3 ./pluto_sirocco_dir_iso.py
```

**Submit job:**
```bash
sbatch job.slurm
```

**Monitor:**
```bash
myqueue
```

**Cancel:**
```bash
scancel JOBID
```

---

## 👩‍💻 Contributors

- Amin Mosallanezhad  
- SIROCCO & PLUTO user community

## 📄 License

This repository is distributed under the **GPL-2.0 license**, in line with PLUTO and SIROCCO.

---

## 📬 Contact

For questions, please contact:  
**Amin Mosallanezhad** — a.mosallanezhad@soton.ac.uk
