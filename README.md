# PLUTO-SIROCCO v1.0

**Monte Carlo Radiation Hydrodynamics (MC-RHD) of Line-Driven Disc Winds**

This repository contains a version of the PLUTO v4.4.3 astrophysical fluid dynamics code that has been custom-modified to interface with the SIROCCO radiative transfer code and a separate CAK force multiplier solver**. Together, this software should allow the user to perform radiation–hydrodynamics (RHD) simulations of line-driven disc winds (or in principle, radiative or thermal type of disc wind). 

We implement a modular, operator-split approach to simulate disc winds:
- **PLUTO** handles the hydrodynamic evolution of the gas (density, velocity, pressure, temperature).
- **SIROCCO** performs frequency-dependent **Monte Carlo radiative transfer**, calculating the local radiation field and ionization states of the gas.
- **CAK (Castor-Abbott-Klein)** solver calculates the **force multiplier** — a function that quantifies how much the radiation pressure is boosted due to line absorption.

These three modules are called in sequence at each radiation timestep, forming a coupled loop. For further detail, please consult the following papers:

_Monte-Carlo radiation hydrodynamic simulations of line-driven disc winds: relaxing the isothermal approximation_, Mosallanezhad et al., MNRAS in press

_State-of-the-art simulations of line-driven accretion disc winds: realistic radiation hydrodynamics leads to weaker outflows_, Higginbottom, Scepi et al., MNRAS, Volume 527, Issue 3, pp.9236-9249, [arXiv:2312.06042](https://arxiv.org/abs/2312.06042)

**Please cite these papers if you use the software, together with the [SIROCCO release paper](https://ui.adsabs.harvard.edu/abs/2025MNRAS.536..879M/abstract) and PLUTO code paper ([Mignone et al. 2007](https://ui.adsabs.harvard.edu/abs/2007ApJS..170..228M/abstract)).** 

## ⚙️ Code Structure Summary

- `cak_v3/` — `atomic_models.txt`, `transitiondata.txt`, `cak_v3` executable  
- `pluto_sirocco_*.py` — Python drivers and control scripts  
- `init.c`, `line_connect.c`, `cooling.c`, etc. — Custom PLUTO source code

Setups for running models from Mosallanezhad et al. (including further parameter descriptions) can be found in `Test_Problems/LineDrivenWind`

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

### 2. Install PLUTO-SIROCCO v1.0

```bash
git clone https://github.com/sirocco-rt/pluto-sirocco.git
```

> If you downloaded the `.zip` file manually, unzip it and **rename the extracted folder to `PLUTO`** for consistency with environment variables used in the rest of this guide.


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

## 👩‍💻 Contributors

- Amin Mosallanezhad
- Nicolas Scepi
- Nick Higginbottom 
- SIROCCO collaboration
- PLUTO developers (see main PLUTO repository here) 

## 📄 License

This repository is distributed under the **GPL-2.0 license**, in line with PLUTO and SIROCCO.

---

## 📬 Contact

For questions, please contact:  
**Amin Mosallanezhad** — a.mosallanezhad@soton.ac.uk
